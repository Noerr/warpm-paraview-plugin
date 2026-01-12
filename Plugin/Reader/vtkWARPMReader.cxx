/**
 * @file vtkWARPMReader.cxx
 * @brief Implementation of VTK reader for WARPM HDF5 simulation output
 *
 * This is the base reader for physical-space (x, y, z) WARPM data.
 * For phase space (Vlasov-Maxwell) data, use vtkWARPMPhaseSpaceReader.
 */

#include "vtkWARPMReader.h"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArraySelection.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>

#include "vtk_hdf5.h"

// For VTK's embedded Python interpreter (string-based communication only)
#include <vtkPythonInterpreter.h>

#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <cmath>

vtkStandardNewMacro(vtkWARPMReader);

//----------------------------------------------------------------------------
// Helper to read a scalar integer attribute
static bool ReadIntAttribute(hid_t loc, const char* name, int& value)
{
  hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
  if (attr < 0) return false;
  herr_t status = H5Aread(attr, H5T_NATIVE_INT, &value);
  H5Aclose(attr);
  return status >= 0;
}

// Helper to read an integer array attribute
static bool ReadIntArrayAttribute(hid_t loc, const char* name, std::vector<int>& values)
{
  hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
  if (attr < 0) return false;

  hid_t space = H5Aget_space(attr);
  int ndims = H5Sget_simple_extent_ndims(space);
  if (ndims != 1) {
    H5Sclose(space);
    H5Aclose(attr);
    return false;
  }

  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, nullptr);
  values.resize(dims[0]);

  herr_t status = H5Aread(attr, H5T_NATIVE_INT, values.data());
  H5Sclose(space);
  H5Aclose(attr);
  return status >= 0;
}

// Helper to read a scalar double attribute
static bool ReadDoubleAttribute(hid_t loc, const char* name, double& value)
{
  hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
  if (attr < 0) return false;
  herr_t status = H5Aread(attr, H5T_NATIVE_DOUBLE, &value);
  H5Aclose(attr);
  return status >= 0;
}

// Helper to read a string attribute
static bool ReadStringAttribute(hid_t loc, const char* name, std::string& value)
{
  hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
  if (attr < 0) return false;

  hid_t atype = H5Aget_type(attr);
  size_t size = H5Tget_size(atype);

  if (H5Tis_variable_str(atype)) {
    char* str = nullptr;
    hid_t memtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(memtype, H5T_VARIABLE);
    herr_t status = H5Aread(attr, memtype, &str);
    if (status >= 0 && str) {
      value = str;
      H5free_memory(str);
    }
    H5Tclose(memtype);
  } else {
    std::vector<char> buf(size + 1, '\0');
    herr_t status = H5Aread(attr, atype, buf.data());
    if (status >= 0) {
      value = buf.data();
    }
  }

  H5Tclose(atype);
  H5Aclose(attr);
  return true;
}

// Helper to read a string array attribute
static bool ReadStringArrayAttribute(hid_t loc, const char* name, std::vector<std::string>& values)
{
  hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
  if (attr < 0) return false;

  hid_t space = H5Aget_space(attr);
  hid_t atype = H5Aget_type(attr);

  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, nullptr);

  values.clear();

  if (H5Tis_variable_str(atype)) {
    std::vector<char*> strs(dims[0], nullptr);
    hid_t memtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(memtype, H5T_VARIABLE);
    H5Aread(attr, memtype, strs.data());
    for (size_t i = 0; i < dims[0]; ++i) {
      if (strs[i]) {
        values.push_back(strs[i]);
        H5free_memory(strs[i]);
      }
    }
    H5Tclose(memtype);
  }

  H5Tclose(atype);
  H5Sclose(space);
  H5Aclose(attr);
  return true;
}

// Check if a domain is physical-only (no velocity dimensions)
// Velocity dimensions have coordinate names starting with 'v' or 'w'
static bool IsPhysicalOnlyDomain(hid_t domainGroup)
{
  std::vector<std::string> coordNames;
  if (!ReadStringArrayAttribute(domainGroup, "CoordinateNames", coordNames))
  {
    // If no CoordinateNames attribute, assume physical (legacy compatibility)
    return true;
  }

  for (const auto& name : coordNames)
  {
    if (!name.empty() && (name[0] == 'v' || name[0] == 'w'))
    {
      return false; // Has velocity dimension
    }
  }
  return true; // All coordinates are physical (x, y, z, etc.)
}

// Check if a variable's domain is physical-only
static bool IsVariableOnPhysicalDomain(hid_t file, hid_t varGroup)
{
  std::string domainName;
  if (!ReadStringAttribute(varGroup, "OnDomain", domainName))
  {
    return false;
  }

  std::string domainPath = "/domains/" + domainName;
  hid_t domGroup = H5Gopen(file, domainPath.c_str(), H5P_DEFAULT);
  if (domGroup < 0) return false;

  bool isPhysical = IsPhysicalOnlyDomain(domGroup);
  H5Gclose(domGroup);
  return isPhysical;
}

//----------------------------------------------------------------------------
// Python-based coordinate expression evaluator (vectorized with numpy)
//
// WARPM VertexCoordinateExpressions are functions x(k) where k is the vertex
// index. These can be arbitrary C-like math expressions evaluated by WARPM's
// HIP language, including polynomials and math functions.
//
// The vertex index k ranges from startIndex to startIndex + numCells (inclusive).
// startIndex is NOT always 0 - it can be negative (e.g., velocity dimensions
// often use symmetric ranges like -3 to +3).
//
// Examples:
//   "x=-0.000542 + 0.001115*k;"                    (linear)
//   "x=-0.000542 + 0.001115*k - 0.0001*k*k;"       (quadratic)
//   "x= k*(0.0148*k*k + 0.2);"                     (cubic, factored)
//   "x=sin(k * 0.1) * 2.0;"                        (trigonometric)
//
// We use numpy for vectorized evaluation - single Python call for all vertices.
//----------------------------------------------------------------------------

// Compute all vertex positions for a dimension by evaluating the coordinate
// expression using numpy vectorization. Returns positions for vertices at
// k = startIndex, startIndex+1, ..., startIndex+numCells (numCells+1 values).
//
// Uses string-based communication with Python to avoid linking against Python directly.
// The result is passed back as a comma-separated string in a global variable.
static std::vector<double> ComputeVertexPositions(
  const std::string& expr, int startIndex, int numCells)
{
  int numVertices = numCells + 1;
  std::vector<double> positions(numVertices);

  // Ensure Python is initialized
  if (!vtkPythonInterpreter::IsInitialized())
  {
    vtkPythonInterpreter::Initialize();
  }

  // Extract RHS of expression (after "x=" or similar)
  std::string rhs = expr;
  size_t eqPos = expr.find('=');
  if (eqPos != std::string::npos)
  {
    rhs = expr.substr(eqPos + 1);
  }

  // Remove trailing semicolon and whitespace
  while (!rhs.empty() && (rhs.back() == ';' || rhs.back() == ' ' || rhs.back() == '\t'))
  {
    rhs.pop_back();
  }

  // Build Python code for vectorized evaluation
  // Store result as comma-separated string in a global variable
  std::ostringstream pyCode;
  pyCode << "import numpy as np\n";
  pyCode << "from math import *\n";
  pyCode << "k = np.arange(" << startIndex << ", " << (startIndex + numCells + 1) << ", dtype=np.float64)\n";
  pyCode << "_warpm_result = " << rhs << "\n";
  pyCode << "if np.isscalar(_warpm_result):\n";
  pyCode << "    _warpm_result = np.full(" << numVertices << ", _warpm_result)\n";
  pyCode << "_warpm_result_str = ','.join(f'{v:.17g}' for v in _warpm_result)\n";

  // Execute the Python code
  vtkPythonInterpreter::RunSimpleString(pyCode.str().c_str());

  // Retrieve the result string using vtkPythonInterpreter
  // We need to write it to a file and read it back (no direct API access)
  // Alternative: use a temporary file or pipe
  // Simpler approach: use vtkPythonInterpreter to write to a known location

  // Write result to temp file and read it back
  std::string tempFile = "/tmp/_warpm_vertex_positions.txt";
  std::ostringstream writeCode;
  writeCode << "with open('" << tempFile << "', 'w') as f:\n";
  writeCode << "    f.write(_warpm_result_str)\n";
  writeCode << "del _warpm_result, _warpm_result_str\n";

  vtkPythonInterpreter::RunSimpleString(writeCode.str().c_str());

  // Read the result from the temp file
  std::ifstream resultFile(tempFile);
  if (resultFile.is_open())
  {
    std::string resultStr;
    std::getline(resultFile, resultStr);
    resultFile.close();

    // Parse comma-separated values
    std::istringstream iss(resultStr);
    std::string token;
    size_t idx = 0;
    while (std::getline(iss, token, ',') && idx < static_cast<size_t>(numVertices))
    {
      try
      {
        positions[idx++] = std::stod(token);
      }
      catch (...)
      {
        // Keep default value of 0.0
      }
    }
  }

  return positions;
}

// Gauss-Lobatto-Legendre nodes on [-1, 1] for orders 1-8
static std::vector<double> GetGLLNodes(int order)
{
  std::vector<double> nodes;
  switch (order) {
    case 1:
      nodes = {-1.0, 1.0};
      break;
    case 2:
      nodes = {-1.0, 0.0, 1.0};
      break;
    case 3:
      nodes = {-1.0, -0.4472135954999579, 0.4472135954999579, 1.0};
      break;
    case 4:
      nodes = {-1.0, -0.6546536707079771, 0.0, 0.6546536707079771, 1.0};
      break;
    case 5:
      nodes = {-1.0, -0.7650553239294647, -0.2852315164806451,
               0.2852315164806451, 0.7650553239294647, 1.0};
      break;
    case 6:
      nodes = {-1.0, -0.8302238962785669, -0.4688487934707142, 0.0,
               0.4688487934707142, 0.8302238962785669, 1.0};
      break;
    case 7:
      nodes = {-1.0, -0.8717401485096066, -0.5917001814331423,
               -0.2092992179024789, 0.2092992179024789,
               0.5917001814331423, 0.8717401485096066, 1.0};
      break;
    case 8:
      nodes = {-1.0, -0.8997579954114602, -0.6771862795107378,
               -0.3631174638261782, 0.0, 0.3631174638261782,
               0.6771862795107378, 0.8997579954114602, 1.0};
      break;
    default:
      // For higher orders, use uniform spacing as fallback
      nodes.resize(order + 1);
      for (int i = 0; i <= order; ++i) {
        nodes[i] = -1.0 + 2.0 * i / order;
      }
  }
  return nodes;
}

// Build mapping from WARPM row-major node index to VTK Lagrange node index
// WARPM: row-major tensor product (y varies fastest, then x)
//   linear index = x * nodesY + y
// VTK Lagrange: corners, then edges (ccw), then interior (row-major)
static std::vector<int> BuildWarpmToVTKLagrangeMapping(int orderX, int orderY)
{
  int nodesX = orderX + 1;
  int nodesY = orderY + 1;
  int totalNodes = nodesX * nodesY;

  // Map from (x, y) grid position to VTK Lagrange index
  std::vector<int> xyToVTK(totalNodes, -1);

  int vtkIdx = 0;

  // 1. Four corners (counter-clockwise from bottom-left)
  // Corner 0: (0, 0) - bottom-left
  // Corner 1: (orderX, 0) - bottom-right
  // Corner 2: (orderX, orderY) - top-right
  // Corner 3: (0, orderY) - top-left
  xyToVTK[0 * nodesY + 0] = vtkIdx++;                    // (0, 0)
  xyToVTK[orderX * nodesY + 0] = vtkIdx++;              // (orderX, 0)
  xyToVTK[orderX * nodesY + orderY] = vtkIdx++;         // (orderX, orderY)
  xyToVTK[0 * nodesY + orderY] = vtkIdx++;              // (0, orderY)

  // 2. Edge nodes (interior nodes on each edge, corners excluded)
  // VTK Lagrange edges all go in INCREASING parameter direction
  // (See vtkHigherOrderQuadrilateral::PointIndexFromIJK)
  // Edge 0: bottom (y=0), x from 1 to orderX-1
  for (int x = 1; x < orderX; ++x) {
    xyToVTK[x * nodesY + 0] = vtkIdx++;
  }
  // Edge 1: right (x=orderX), y from 1 to orderY-1
  for (int y = 1; y < orderY; ++y) {
    xyToVTK[orderX * nodesY + y] = vtkIdx++;
  }
  // Edge 2: top (y=orderY), x from 1 to orderX-1 (NOT reversed!)
  for (int x = 1; x < orderX; ++x) {
    xyToVTK[x * nodesY + orderY] = vtkIdx++;
  }
  // Edge 3: left (x=0), y from 1 to orderY-1 (NOT reversed!)
  for (int y = 1; y < orderY; ++y) {
    xyToVTK[0 * nodesY + y] = vtkIdx++;
  }

  // 3. Interior nodes (VTK convention: x/first parametric direction varies fastest)
  for (int y = 1; y < orderY; ++y) {
    for (int x = 1; x < orderX; ++x) {
      xyToVTK[x * nodesY + y] = vtkIdx++;
    }
  }

  // Now invert to get WARPM-to-VTK mapping
  // WARPM linear index = x * nodesY + y (same as our xyToVTK indexing)
  // So xyToVTK IS the WARPM-to-VTK mapping!
  return xyToVTK;
}

//----------------------------------------------------------------------------
vtkWARPMReader::vtkWARPMReader()
  : FileName(nullptr)
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1); // Single port for physical space
  this->PointDataArraySelection = vtkSmartPointer<vtkDataArraySelection>::New();
}

//----------------------------------------------------------------------------
vtkWARPMReader::~vtkWARPMReader()
{
  delete[] this->FileName;
}

//----------------------------------------------------------------------------
int vtkWARPMReader::FillOutputPortInformation(int port, vtkInformation* info)
{
  // Single output port for physical space mesh
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}

//----------------------------------------------------------------------------
// Static helper: Check if an HDF5 file has any phase space variables
// (variables on domains with velocity coordinates starting with 'v' or 'w')
static bool HasPhaseSpaceVariables(const char* fname)
{
  hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) return false;

  bool hasPhaseSpace = false;

  hid_t domainsGroup = H5Gopen(file, "/domains", H5P_DEFAULT);
  hid_t varsGroup = H5Gopen(file, "/variables", H5P_DEFAULT);

  if (domainsGroup >= 0 && varsGroup >= 0)
  {
    // Collect variable names first
    std::vector<std::string> varNames;
    H5Literate(varsGroup, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr,
      [](hid_t, const char* name, const H5L_info_t*, void* opdata) -> herr_t {
        static_cast<std::vector<std::string>*>(opdata)->push_back(name);
        return 0;
      }, &varNames);

    // Check each variable's domain for velocity coordinates
    for (const auto& varName : varNames)
    {
      if (hasPhaseSpace) break;

      hid_t varGroup = H5Gopen(varsGroup, varName.c_str(), H5P_DEFAULT);
      if (varGroup >= 0)
      {
        std::string domainName;
        ReadStringAttribute(varGroup, "OnDomain", domainName);
        H5Gclose(varGroup);

        if (!domainName.empty())
        {
          hid_t domGroup = H5Gopen(domainsGroup, domainName.c_str(), H5P_DEFAULT);
          if (domGroup >= 0)
          {
            std::vector<std::string> coordNames;
            ReadStringArrayAttribute(domGroup, "CoordinateNames", coordNames);
            H5Gclose(domGroup);

            for (const auto& coord : coordNames)
            {
              if (!coord.empty() && (coord[0] == 'v' || coord[0] == 'w'))
              {
                hasPhaseSpace = true;
                break;
              }
            }
          }
        }
      }
    }
  }

  if (domainsGroup >= 0) H5Gclose(domainsGroup);
  if (varsGroup >= 0) H5Gclose(varsGroup);
  H5Fclose(file);

  return hasPhaseSpace;
}

//----------------------------------------------------------------------------
// Static helper: Check if an HDF5 file has WARPM structure
static bool IsWarpmHDF5File(const char* fname)
{
  hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0)
  {
    return false;
  }

  // Check for required WARPM groups: /domains, /variables, /timeData
  bool hasRequired = true;

  hid_t domainsGroup = H5Gopen(file, "/domains", H5P_DEFAULT);
  if (domainsGroup < 0)
  {
    hasRequired = false;
  }
  else
  {
    H5Gclose(domainsGroup);
  }

  hid_t varsGroup = H5Gopen(file, "/variables", H5P_DEFAULT);
  if (varsGroup < 0)
  {
    hasRequired = false;
  }
  else
  {
    H5Gclose(varsGroup);
  }

  hid_t timeGroup = H5Gopen(file, "/timeData", H5P_DEFAULT);
  if (timeGroup < 0)
  {
    hasRequired = false;
  }
  else
  {
    H5Gclose(timeGroup);
  }

  H5Fclose(file);
  return hasRequired;
}

//----------------------------------------------------------------------------
int vtkWARPMReader::CanReadFile(const char* fname)
{
  if (!fname || strlen(fname) == 0)
  {
    return 0;
  }

  std::filesystem::path filePath(fname);
  std::string ext = filePath.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  if (ext == ".h5" || ext == ".hdf5")
  {
    // Direct HDF5 file - check for WARPM structure
    // Return false if file has phase space variables (use WARPMPhaseSpaceReader instead)
    if (!IsWarpmHDF5File(fname))
    {
      return 0;
    }
    return HasPhaseSpaceVariables(fname) ? 0 : 1;
  }
  else if (ext == ".warpm")
  {
    // Metadata file - parse to get first H5 and validate it
    std::filesystem::path baseDir = filePath.parent_path();
    if (baseDir.empty())
    {
      baseDir = ".";
    }

    std::ifstream file(fname);
    if (!file.is_open())
    {
      return 0;
    }

    std::string line;
    while (std::getline(file, line))
    {
      // Trim whitespace
      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);

      // Skip empty lines and comments
      if (line.empty() || line[0] == '#')
      {
        continue;
      }

      // Handle glob pattern - find first matching file
      if (line.find('*') != std::string::npos || line.find('?') != std::string::npos)
      {
        std::filesystem::path pattern = baseDir / line;
        std::filesystem::path parentDir = pattern.parent_path();
        std::string globPattern = pattern.filename().string();

        if (std::filesystem::exists(parentDir) && std::filesystem::is_directory(parentDir))
        {
          for (const auto& entry : std::filesystem::directory_iterator(parentDir))
          {
            if (entry.is_regular_file())
            {
              std::string filename = entry.path().filename().string();
              size_t starPos = globPattern.find('*');
              if (starPos != std::string::npos)
              {
                std::string prefix = globPattern.substr(0, starPos);
                std::string suffix = globPattern.substr(starPos + 1);
                if (filename.size() >= prefix.size() + suffix.size() &&
                    filename.substr(0, prefix.size()) == prefix &&
                    filename.substr(filename.size() - suffix.size()) == suffix)
                {
                  // Found a matching file - validate it
                  std::string h5path = entry.path().string();
                  if (!IsWarpmHDF5File(h5path.c_str()))
                  {
                    return 0;
                  }
                  // Return false if file has phase space variables
                  return HasPhaseSpaceVariables(h5path.c_str()) ? 0 : 1;
                }
              }
            }
          }
        }
      }
      else
      {
        // Explicit filename
        std::filesystem::path framePath = baseDir / line;
        if (std::filesystem::exists(framePath))
        {
          std::string h5path = framePath.string();
          if (!IsWarpmHDF5File(h5path.c_str()))
          {
            return 0;
          }
          // Return false if file has phase space variables
          return HasPhaseSpaceVariables(h5path.c_str()) ? 0 : 1;
        }
      }
    }

    return 0; // No valid H5 file found
  }

  return 0; // Unknown extension
}

//----------------------------------------------------------------------------
void vtkWARPMReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: " << (this->FileName ? this->FileName : "(none)") << "\n";
  os << indent << "Number of frames: " << this->FrameFiles.size() << "\n";
}

//----------------------------------------------------------------------------
vtkDataArraySelection* vtkWARPMReader::GetPointDataArraySelection()
{
  return this->PointDataArraySelection;
}

//----------------------------------------------------------------------------
int vtkWARPMReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkWARPMReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkWARPMReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkWARPMReader::SetPointArrayStatus(const char* name, int status)
{
  if (status)
  {
    this->PointDataArraySelection->EnableArray(name);
  }
  else
  {
    this->PointDataArraySelection->DisableArray(name);
  }
  this->Modified();
}

//----------------------------------------------------------------------------
bool vtkWARPMReader::ParseWarpmFile()
{
  if (!this->FileName || strlen(this->FileName) == 0)
  {
    vtkErrorMacro("No filename specified");
    return false;
  }

  this->FrameFiles.clear();

  // Get base directory for resolving relative paths
  std::filesystem::path filePath(this->FileName);
  this->BaseDirectory = filePath.parent_path().string();
  if (this->BaseDirectory.empty())
  {
    this->BaseDirectory = ".";
  }

  std::ifstream file(this->FileName);
  if (!file.is_open())
  {
    vtkErrorMacro("Could not open file: " << this->FileName);
    return false;
  }

  std::string line;
  while (std::getline(file, line))
  {
    // Trim whitespace
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);

    // Skip empty lines and comments
    if (line.empty() || line[0] == '#')
    {
      continue;
    }

    // Check if it's a glob pattern (contains * or ?)
    if (line.find('*') != std::string::npos || line.find('?') != std::string::npos)
    {
      // Handle glob pattern
      std::filesystem::path pattern = std::filesystem::path(this->BaseDirectory) / line;
      std::filesystem::path parentDir = pattern.parent_path();
      std::string globPattern = pattern.filename().string();

      if (std::filesystem::exists(parentDir) && std::filesystem::is_directory(parentDir))
      {
        std::vector<std::string> matches;
        for (const auto& entry : std::filesystem::directory_iterator(parentDir))
        {
          if (entry.is_regular_file())
          {
            std::string filename = entry.path().filename().string();
            size_t starPos = globPattern.find('*');
            if (starPos != std::string::npos)
            {
              std::string prefix = globPattern.substr(0, starPos);
              std::string suffix = globPattern.substr(starPos + 1);
              if (filename.size() >= prefix.size() + suffix.size() &&
                  filename.substr(0, prefix.size()) == prefix &&
                  filename.substr(filename.size() - suffix.size()) == suffix)
              {
                matches.push_back(entry.path().string());
              }
            }
          }
        }
        std::sort(matches.begin(), matches.end());
        for (const auto& match : matches)
        {
          this->FrameFiles.push_back(match);
        }
      }
      else
      {
        vtkWarningMacro("Directory does not exist: " << parentDir.string());
      }
    }
    else
    {
      std::filesystem::path framePath = std::filesystem::path(this->BaseDirectory) / line;
      if (std::filesystem::exists(framePath))
      {
        this->FrameFiles.push_back(framePath.string());
      }
      else
      {
        vtkWarningMacro("Frame file not found: " << framePath.string());
      }
    }
  }

  return !this->FrameFiles.empty();
}

//----------------------------------------------------------------------------
bool vtkWARPMReader::ReadTimeValues()
{
  this->TimeValues.clear();
  this->TimeValues.reserve(this->FrameFiles.size());

  for (const auto& frameFile : this->FrameFiles)
  {
    hid_t file = H5Fopen(frameFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0)
    {
      vtkWarningMacro("Could not open " << frameFile << " to read time");
      this->TimeValues.push_back(static_cast<double>(this->TimeValues.size()));
      continue;
    }

    // Open /timeData group
    hid_t timeGroup = H5Gopen(file, "/timeData", H5P_DEFAULT);
    if (timeGroup < 0)
    {
      vtkWarningMacro("No /timeData group in " << frameFile);
      this->TimeValues.push_back(static_cast<double>(this->TimeValues.size()));
      H5Fclose(file);
      continue;
    }

    // Read time attribute
    double timeValue = static_cast<double>(this->TimeValues.size()); // fallback
    ReadDoubleAttribute(timeGroup, "time", timeValue);

    this->TimeValues.push_back(timeValue);

    H5Gclose(timeGroup);
    H5Fclose(file);
  }

  return !this->TimeValues.empty();
}

//----------------------------------------------------------------------------
int vtkWARPMReader::FindClosestFrame(double requestedTime) const
{
  if (this->TimeValues.empty())
  {
    return 0;
  }

  int closestFrame = 0;
  double minDiff = std::abs(this->TimeValues[0] - requestedTime);

  for (size_t i = 1; i < this->TimeValues.size(); ++i)
  {
    double diff = std::abs(this->TimeValues[i] - requestedTime);
    if (diff < minDiff)
    {
      minDiff = diff;
      closestFrame = static_cast<int>(i);
    }
  }

  return closestFrame;
}

//----------------------------------------------------------------------------
int vtkWARPMReader::RequestInformation(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  if (!this->FileName || strlen(this->FileName) == 0)
  {
    vtkErrorMacro("No filename specified");
    return 0;
  }

  // Determine file type from extension
  std::filesystem::path filePath(this->FileName);
  std::string ext = filePath.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  if (ext == ".h5" || ext == ".hdf5")
  {
    // Direct HDF5 file - use it as single frame
    this->FrameFiles.clear();
    this->FrameFiles.push_back(this->FileName);
    this->BaseDirectory = filePath.parent_path().string();
    if (this->BaseDirectory.empty())
    {
      this->BaseDirectory = ".";
    }
  }
  else if (ext == ".warpm")
  {
    // Metadata file - parse for frame list
    if (!this->ParseWarpmFile())
    {
      return 0;
    }
  }
  else
  {
    vtkErrorMacro("Unsupported file extension: " << ext);
    return 0;
  }

  // Read actual time values from HDF5 files
  this->ReadTimeValues();

  // Enumerate available variables from first HDF5 file
  // Only add variables on physical-only domains (skip phase space variables)
  if (!this->FrameFiles.empty())
  {
    hid_t file = H5Fopen(this->FrameFiles[0].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file >= 0)
    {
      hid_t varsGroup = H5Gopen(file, "/variables", H5P_DEFAULT);
      if (varsGroup >= 0)
      {
        // Collect all variable names
        std::vector<std::string> varNames;
        H5Literate(varsGroup, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr,
          [](hid_t, const char* name, const H5L_info_t*, void* opdata) -> herr_t {
            static_cast<std::vector<std::string>*>(opdata)->push_back(name);
            return 0; // Continue iterating
          }, &varNames);

        // Add only variables on physical-only domains (no velocity dimensions)
        for (const auto& name : varNames)
        {
          // Check if variable is on a physical-only domain (no velocity dimensions)
          hid_t varGroup = H5Gopen(varsGroup, name.c_str(), H5P_DEFAULT);
          if (varGroup >= 0)
          {
            if (IsVariableOnPhysicalDomain(file, varGroup))
            {
              if (!this->PointDataArraySelection->ArrayExists(name.c_str()))
              {
                this->PointDataArraySelection->AddArray(name.c_str());
                this->PointDataArraySelection->EnableArray(name.c_str());
              }
            }
            H5Gclose(varGroup);
          }
        }

        H5Gclose(varsGroup);
      }
      H5Fclose(file);
    }
  }

  // Report available time steps
  if (!this->TimeValues.empty())
  {
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    if (outInfo)
    {
      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                   this->TimeValues.data(),
                   static_cast<int>(this->TimeValues.size()));

      double timeRange[2] = { this->TimeValues.front(), this->TimeValues.back() };
      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    }
  }
  else
  {
    vtkWarningMacro("No time values available");
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkWARPMReader::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid* output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!output)
  {
    vtkErrorMacro("Could not get output");
    return 0;
  }

  if (this->FrameFiles.empty())
  {
    vtkErrorMacro("No frame files available");
    return 0;
  }

  // Determine which frame to read based on requested time
  int frameIndex = 0;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
  {
    double requestedTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    frameIndex = this->FindClosestFrame(requestedTime);
  }

  const std::string& frameFile = this->FrameFiles[frameIndex];

  // Open HDF5 file
  hid_t file = H5Fopen(frameFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0)
  {
    vtkErrorMacro("Could not open HDF5 file: " << frameFile);
    return 0;
  }

  // Open groups
  hid_t varsGroup = H5Gopen(file, "/variables", H5P_DEFAULT);
  if (varsGroup < 0)
  {
    vtkErrorMacro("Could not open variables group");
    H5Fclose(file);
    return 0;
  }

  // Collect enabled variable names that are on physical-only domains
  std::vector<std::string> allVarNames;
  H5Literate(varsGroup, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr,
    [](hid_t, const char* name, const H5L_info_t*, void* opdata) -> herr_t {
      static_cast<std::vector<std::string>*>(opdata)->push_back(name);
      return 0;
    }, &allVarNames);

  std::vector<std::string> enabledVars;
  std::string physicalDomainName; // Will be set from the first physical variable found
  for (const auto& name : allVarNames)
  {
    if (this->PointDataArraySelection->ArrayIsEnabled(name.c_str()))
    {
      // Safety check: verify variable is on a physical-only domain
      hid_t varGroup = H5Gopen(varsGroup, name.c_str(), H5P_DEFAULT);
      if (varGroup >= 0)
      {
        if (IsVariableOnPhysicalDomain(file, varGroup))
        {
          std::string varDomainName;
          ReadStringAttribute(varGroup, "OnDomain", varDomainName);

          if (physicalDomainName.empty())
          {
            // First physical variable - use its domain
            physicalDomainName = varDomainName;
            enabledVars.push_back(name);
          }
          else if (varDomainName == physicalDomainName)
          {
            // Same domain as first variable
            enabledVars.push_back(name);
          }
          else
          {
            // Different domain - warn and skip
            vtkWarningMacro("Skipping variable '" << name << "' on domain '"
              << varDomainName << "' (differs from '" << physicalDomainName << "')");
          }
        }
        H5Gclose(varGroup);
      }
    }
  }

  if (enabledVars.empty())
  {
    vtkWarningMacro("No physical-space variables enabled (use WARPM Phase Space Reader for phase space data)");
    H5Gclose(varsGroup);
    H5Fclose(file);
    return 1;
  }

  // Read domain info using the domain name from the first enabled variable
  std::string domainPath = "/domains/" + physicalDomainName;
  hid_t domain = H5Gopen(file, domainPath.c_str(), H5P_DEFAULT);
  if (domain < 0)
  {
    vtkErrorMacro("Could not open domain group: " << physicalDomainName);
    H5Gclose(varsGroup);
    H5Fclose(file);
    return 0;
  }

  int ndims = 0;
  std::vector<int> numCells;
  std::vector<int> startIndices;
  std::vector<std::string> coordExprs;

  ReadIntAttribute(domain, "ndims", ndims);
  ReadIntArrayAttribute(domain, "numCells", numCells);
  ReadIntArrayAttribute(domain, "startIndices", startIndices);
  ReadStringArrayAttribute(domain, "VertexCoordinateExpressions", coordExprs);

  H5Gclose(domain);

  if (ndims < 2 || numCells.size() < 2 || coordExprs.size() < 2)
  {
    vtkErrorMacro("Invalid domain dimensions");
    H5Gclose(varsGroup);
    H5Fclose(file);
    return 0;
  }

  // Default startIndices to 0 if not present
  if (startIndices.size() < 2)
  {
    startIndices.resize(ndims, 0);
  }

  // Evaluate coordinate expressions to get vertex positions
  // vertexX[i] is the x-coordinate of vertex i (where i is relative to cell index)
  std::vector<double> vertexX = ComputeVertexPositions(coordExprs[0], startIndices[0], numCells[0]);
  std::vector<double> vertexY = ComputeVertexPositions(coordExprs[1], startIndices[1], numCells[1]);

  // Read metadata from first enabled variable (for geometry)
  hid_t firstVarGroup = H5Gopen(varsGroup, enabledVars[0].c_str(), H5P_DEFAULT);
  std::vector<int> elementOrder;
  int entriesPerElement = 0;

  ReadIntArrayAttribute(firstVarGroup, "ElementOrder", elementOrder);
  ReadIntAttribute(firstVarGroup, "EntriesPerElement", entriesPerElement);

  // Get data dimensions from first variable
  hid_t firstDataset = H5Dopen(firstVarGroup, "data", H5P_DEFAULT);
  hid_t firstDataspace = H5Dget_space(firstDataset);
  int dataNdims = H5Sget_simple_extent_ndims(firstDataspace);
  std::vector<hsize_t> dataDims(dataNdims);
  H5Sget_simple_extent_dims(firstDataspace, dataDims.data(), nullptr);
  H5Sclose(firstDataspace);
  H5Dclose(firstDataset);
  H5Gclose(firstVarGroup);

  int dataNx = static_cast<int>(dataDims[0]);
  int dataNy = static_cast<int>(dataDims[1]);
  int nodesPerElement = static_cast<int>(dataDims[2]);

  // Determine element order from nodes per element
  int orderX = (elementOrder.size() > 0) ? elementOrder[0] : static_cast<int>(std::sqrt(nodesPerElement)) - 1;
  int orderY = (elementOrder.size() > 1) ? elementOrder[1] : orderX;
  int nodesX = orderX + 1;
  int nodesY = orderY + 1;

  // Check if we can use cached geometry
  bool canUseCache = this->GeometryCached &&
                     this->CachedNx == dataNx &&
                     this->CachedNy == dataNy &&
                     this->CachedNodesPerElement == nodesPerElement;

  if (!canUseCache)
  {
    // Build geometry from scratch
    auto gllX = GetGLLNodes(orderX);
    auto gllY = GetGLLNodes(orderY);
    this->CachedWarpmToVTK = BuildWarpmToVTKLagrangeMapping(orderX, orderY);

    // Create points
    this->CachedPoints = vtkSmartPointer<vtkPoints>::New();
    int totalNodes = dataNx * dataNy * nodesPerElement;
    this->CachedPoints->SetNumberOfPoints(totalNodes);

    // Fill point coordinates using evaluated vertex positions
    for (int ey = 0; ey < dataNy; ++ey)
    {
      for (int ex = 0; ex < dataNx; ++ex)
      {
        // Cell bounds from pre-computed vertex positions
        double elemXMin = vertexX[ex];
        double elemXMax = vertexX[ex + 1];
        double elemYMin = vertexY[ey];
        double elemYMax = vertexY[ey + 1];

        int elemStartIdx = (ey * dataNx + ex) * nodesPerElement;

        for (int ix = 0; ix < nodesX; ++ix)
        {
          for (int iy = 0; iy < nodesY; ++iy)
          {
            int warpmLocalIdx = ix * nodesY + iy;
            int vtkLocalIdx = this->CachedWarpmToVTK[warpmLocalIdx];
            int vtkPointIdx = elemStartIdx + vtkLocalIdx;

            double xi = gllX[ix];
            double eta = gllY[iy];
            double x = elemXMin + (xi + 1.0) * 0.5 * (elemXMax - elemXMin);
            double y = elemYMin + (eta + 1.0) * 0.5 * (elemYMax - elemYMin);

            this->CachedPoints->SetPoint(vtkPointIdx, x, y, 0.0);
          }
        }
      }
    }

    // Create cells
    this->CachedCells = vtkSmartPointer<vtkCellArray>::New();
    for (int ey = 0; ey < dataNy; ++ey)
    {
      for (int ex = 0; ex < dataNx; ++ex)
      {
        int elemStartIdx = (ey * dataNx + ex) * nodesPerElement;
        std::vector<vtkIdType> cellPts(nodesPerElement);
        for (int i = 0; i < nodesPerElement; ++i)
        {
          cellPts[i] = elemStartIdx + i;
        }
        this->CachedCells->InsertNextCell(nodesPerElement, cellPts.data());
      }
    }

    // Cache metadata
    this->CachedNx = dataNx;
    this->CachedNy = dataNy;
    this->CachedNodesPerElement = nodesPerElement;
    this->GeometryCached = true;
  }

  // Set cached geometry on output
  output->SetPoints(this->CachedPoints);
  output->SetCells(VTK_LAGRANGE_QUADRILATERAL, this->CachedCells);

  int totalNodes = dataNx * dataNy * nodesPerElement;

  // Load data for each enabled variable
  for (const auto& varName : enabledVars)
  {
    hid_t varGroup = H5Gopen(varsGroup, varName.c_str(), H5P_DEFAULT);
    if (varGroup < 0)
    {
      vtkWarningMacro("Could not open variable: " << varName);
      continue;
    }

    // Read component names for this variable
    std::vector<std::string> componentNames;
    ReadStringArrayAttribute(varGroup, "ComponentNames", componentNames);

    // Read data
    hid_t dataset = H5Dopen(varGroup, "data", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    int varNdims = H5Sget_simple_extent_ndims(dataspace);
    std::vector<hsize_t> varDims(varNdims);
    H5Sget_simple_extent_dims(dataspace, varDims.data(), nullptr);

    int varNumComponents = static_cast<int>(varDims[3]);

    std::vector<double> data(dataNx * dataNy * nodesPerElement * varNumComponents);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Gclose(varGroup);

    // Create field arrays for each component
    std::vector<vtkSmartPointer<vtkDoubleArray>> fieldArrays;
    for (int c = 0; c < varNumComponents; ++c)
    {
      vtkNew<vtkDoubleArray> arr;
      std::string name = (c < static_cast<int>(componentNames.size()))
                         ? componentNames[c]
                         : varName + "_" + std::to_string(c);
      arr->SetName(name.c_str());
      arr->SetNumberOfTuples(totalNodes);
      fieldArrays.push_back(arr);
    }

    // Fill field data (using cached mapping)
    for (int ey = 0; ey < dataNy; ++ey)
    {
      for (int ex = 0; ex < dataNx; ++ex)
      {
        int elemStartIdx = (ey * dataNx + ex) * nodesPerElement;

        for (int ix = 0; ix < nodesX; ++ix)
        {
          for (int iy = 0; iy < nodesY; ++iy)
          {
            int warpmLocalIdx = ix * nodesY + iy;
            int vtkLocalIdx = this->CachedWarpmToVTK[warpmLocalIdx];
            int vtkPointIdx = elemStartIdx + vtkLocalIdx;

            int warpmDataIdx = ((ex * dataNy + ey) * nodesPerElement + warpmLocalIdx) * varNumComponents;
            for (int c = 0; c < varNumComponents; ++c)
            {
              fieldArrays[c]->SetValue(vtkPointIdx, data[warpmDataIdx + c]);
            }
          }
        }
      }
    }

    // Add field arrays to point data
    for (auto& arr : fieldArrays)
    {
      output->GetPointData()->AddArray(arr);
    }
  }

  H5Gclose(varsGroup);
  H5Fclose(file);

  return 1;
}

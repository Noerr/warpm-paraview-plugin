/**
 * @file vtkWARPMPhaseSpaceReader.cxx
 * @brief Implementation of VTK reader for WARPM HDF5 phase space (Vlasov-Maxwell) data
 */

#include "vtkWARPMPhaseSpaceReader.h"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkDataArraySelection.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkVertex.h>

#include "vtk_hdf5.h"

// For VTK's embedded Python interpreter (string-based communication only)
#include <vtkPythonInterpreter.h>

#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <cmath>

vtkStandardNewMacro(vtkWARPMPhaseSpaceReader);

//----------------------------------------------------------------------------
// Helper to read a string attribute (local copy to avoid linkage issues)
static bool PSReadStringAttribute(hid_t loc, const char* name, std::string& value)
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
static bool PSReadStringArrayAttribute(hid_t loc, const char* name, std::vector<std::string>& values)
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

// Helper to read an integer attribute
static bool PSReadIntAttribute(hid_t loc, const char* name, int& value)
{
  hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
  if (attr < 0) return false;
  herr_t status = H5Aread(attr, H5T_NATIVE_INT, &value);
  H5Aclose(attr);
  return status >= 0;
}

// Helper to read an integer array attribute
static bool PSReadIntArrayAttribute(hid_t loc, const char* name, std::vector<int>& values)
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

// Helper to read a double attribute
static bool PSReadDoubleAttribute(hid_t loc, const char* name, double& value)
{
  hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
  if (attr < 0) return false;
  herr_t status = H5Aread(attr, H5T_NATIVE_DOUBLE, &value);
  H5Aclose(attr);
  return status >= 0;
}

//----------------------------------------------------------------------------
// Python-based coordinate expression evaluator (vectorized with numpy)
//
// WARPM VertexCoordinateExpressions are functions x(k) where k is the vertex
// index. The vertex index k ranges from startIndex to startIndex + numCells.
// startIndex is NOT always 0 - it can be negative (e.g., velocity dimensions
// often use symmetric ranges like -3 to +3).
//----------------------------------------------------------------------------

// Compute all vertex positions for a dimension by evaluating the coordinate
// expression using numpy vectorization. Returns positions for vertices at
// k = startIndex, startIndex+1, ..., startIndex+numCells (numCells+1 values).
//
// Uses string-based communication with Python to avoid linking against Python directly.
static std::vector<double> PSComputeVertexPositions(
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
  // Store result as comma-separated string
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

  // Write result to temp file and read it back
  std::string tempFile = "/tmp/_warpm_ps_vertex_positions.txt";
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
static std::vector<double> PSGetGLLNodes(int order)
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
      nodes.resize(order + 1);
      for (int i = 0; i <= order; ++i) {
        nodes[i] = -1.0 + 2.0 * i / order;
      }
  }
  return nodes;
}

//----------------------------------------------------------------------------
// N-Dimensional Helper Functions for Arbitrary Dimension Support
//----------------------------------------------------------------------------

// Decompose a flat index into N-dimensional indices (row-major, last dim fastest)
// Example: flatIdx=5, sizes={2,3} -> indices={1,2} (5 = 1*3 + 2)
static std::vector<int> PSDecomposeIndex(int flatIdx, const std::vector<int>& sizes)
{
  std::vector<int> indices(sizes.size());
  for (int d = static_cast<int>(sizes.size()) - 1; d >= 0; --d)
  {
    indices[d] = flatIdx % sizes[d];
    flatIdx /= sizes[d];
  }
  return indices;
}

// Compute flat index from N-dimensional indices (row-major, last dim fastest)
// Example: indices={1,2}, sizes={2,3} -> flatIdx=5 (1*3 + 2)
static int PSComputeFlatIndex(const std::vector<int>& indices, const std::vector<int>& sizes)
{
  int flatIdx = 0;
  for (size_t d = 0; d < indices.size(); ++d)
  {
    flatIdx = flatIdx * sizes[d] + indices[d];
  }
  return flatIdx;
}

// Compute total nodes for an N-dimensional element with given orders
// Example: orders={2,3} -> (2+1)*(3+1) = 12 nodes
static int PSComputeTotalNodes(const std::vector<int>& orders)
{
  int total = 1;
  for (int o : orders)
  {
    total *= (o + 1);
  }
  return total;
}

// Get VTK Lagrange cell type for given number of dimensions
static int PSGetLagrangeCellType(int ndims)
{
  switch (ndims)
  {
    case 1: return VTK_LAGRANGE_CURVE;          // Type 68
    case 2: return VTK_LAGRANGE_QUADRILATERAL;  // Type 70
    case 3: return VTK_LAGRANGE_HEXAHEDRON;     // Type 72
    default:
      return -1; // Invalid
  }
}

// Detect if a domain is a phase space domain
static bool PSDetectPhaseSpaceDomain(hid_t domainGroup, int& numPhysical, int& numVelocity,
                                     std::vector<std::string>& coordNames)
{
  numPhysical = 0;
  numVelocity = 0;
  coordNames.clear();

  if (!PSReadStringArrayAttribute(domainGroup, "CoordinateNames", coordNames))
  {
    int ndims = 0;
    if (PSReadIntAttribute(domainGroup, "ndims", ndims))
    {
      numPhysical = ndims;
    }
    return false;
  }

  for (const auto& name : coordNames)
  {
    if (!name.empty() && (name[0] == 'v' || name[0] == 'w'))
    {
      ++numVelocity;
    }
    else
    {
      ++numPhysical;
    }
  }

  return numVelocity > 0;
}

//----------------------------------------------------------------------------
// Node Ordering Mappings: WARPM (row-major) to VTK Lagrange
//----------------------------------------------------------------------------

// 1D mapping: VTK_LAGRANGE_CURVE
// VTK order: endpoints first (0, order), then interior (1 to order-1)
// WARPM order: sequential (0, 1, 2, ..., order)
static std::vector<int> PSBuildWarpmToVTKMapping1D(int order)
{
  int numNodes = order + 1;
  std::vector<int> mapping(numNodes, -1);
  int vtkIdx = 0;

  // 1. Two endpoints
  mapping[0] = vtkIdx++;           // Left endpoint
  mapping[order] = vtkIdx++;       // Right endpoint

  // 2. Interior nodes (1 to order-1)
  for (int i = 1; i < order; ++i)
  {
    mapping[i] = vtkIdx++;
  }

  return mapping;
}

// 2D mapping: VTK_LAGRANGE_QUADRILATERAL (existing function, renamed for consistency)
static std::vector<int> PSBuildWarpmToVTKMapping2D(int orderX, int orderY)
{
  int nodesX = orderX + 1;
  int nodesY = orderY + 1;
  int totalNodes = nodesX * nodesY;

  std::vector<int> xyToVTK(totalNodes, -1);
  int vtkIdx = 0;

  // 1. Four corners (counter-clockwise from bottom-left)
  xyToVTK[0 * nodesY + 0] = vtkIdx++;
  xyToVTK[orderX * nodesY + 0] = vtkIdx++;
  xyToVTK[orderX * nodesY + orderY] = vtkIdx++;
  xyToVTK[0 * nodesY + orderY] = vtkIdx++;

  // 2. Edge nodes (interior nodes on each edge)
  for (int x = 1; x < orderX; ++x) {
    xyToVTK[x * nodesY + 0] = vtkIdx++;
  }
  for (int y = 1; y < orderY; ++y) {
    xyToVTK[orderX * nodesY + y] = vtkIdx++;
  }
  for (int x = 1; x < orderX; ++x) {
    xyToVTK[x * nodesY + orderY] = vtkIdx++;
  }
  for (int y = 1; y < orderY; ++y) {
    xyToVTK[0 * nodesY + y] = vtkIdx++;
  }

  // 3. Interior nodes
  for (int y = 1; y < orderY; ++y) {
    for (int x = 1; x < orderX; ++x) {
      xyToVTK[x * nodesY + y] = vtkIdx++;
    }
  }

  return xyToVTK;
}

// 3D mapping: VTK_LAGRANGE_HEXAHEDRON
// VTK order: 8 corners, 12 edges (interior), 6 faces (interior), volume interior
// WARPM order: row-major x*nodesY*nodesZ + y*nodesZ + z
static std::vector<int> PSBuildWarpmToVTKMapping3D(int orderX, int orderY, int orderZ)
{
  int nodesX = orderX + 1;
  int nodesY = orderY + 1;
  int nodesZ = orderZ + 1;
  int totalNodes = nodesX * nodesY * nodesZ;

  // Helper lambda to compute WARPM flat index from (x,y,z)
  auto warpmIdx = [nodesY, nodesZ](int x, int y, int z) {
    return x * nodesY * nodesZ + y * nodesZ + z;
  };

  std::vector<int> mapping(totalNodes, -1);
  int vtkIdx = 0;

  // 1. Eight corners (VTK ordering)
  mapping[warpmIdx(0, 0, 0)] = vtkIdx++;
  mapping[warpmIdx(orderX, 0, 0)] = vtkIdx++;
  mapping[warpmIdx(orderX, orderY, 0)] = vtkIdx++;
  mapping[warpmIdx(0, orderY, 0)] = vtkIdx++;
  mapping[warpmIdx(0, 0, orderZ)] = vtkIdx++;
  mapping[warpmIdx(orderX, 0, orderZ)] = vtkIdx++;
  mapping[warpmIdx(orderX, orderY, orderZ)] = vtkIdx++;
  mapping[warpmIdx(0, orderY, orderZ)] = vtkIdx++;

  // 2. Twelve edges (interior nodes only, all in increasing parameter direction)
  // Bottom face edges (z=0)
  for (int x = 1; x < orderX; ++x)  // Edge 0-1: y=0, z=0
    mapping[warpmIdx(x, 0, 0)] = vtkIdx++;
  for (int y = 1; y < orderY; ++y)  // Edge 1-2: x=orderX, z=0
    mapping[warpmIdx(orderX, y, 0)] = vtkIdx++;
  for (int x = 1; x < orderX; ++x)  // Edge 3-2: y=orderY, z=0
    mapping[warpmIdx(x, orderY, 0)] = vtkIdx++;
  for (int y = 1; y < orderY; ++y)  // Edge 0-3: x=0, z=0
    mapping[warpmIdx(0, y, 0)] = vtkIdx++;

  // Top face edges (z=orderZ)
  for (int x = 1; x < orderX; ++x)  // Edge 4-5: y=0, z=orderZ
    mapping[warpmIdx(x, 0, orderZ)] = vtkIdx++;
  for (int y = 1; y < orderY; ++y)  // Edge 5-6: x=orderX, z=orderZ
    mapping[warpmIdx(orderX, y, orderZ)] = vtkIdx++;
  for (int x = 1; x < orderX; ++x)  // Edge 7-6: y=orderY, z=orderZ
    mapping[warpmIdx(x, orderY, orderZ)] = vtkIdx++;
  for (int y = 1; y < orderY; ++y)  // Edge 4-7: x=0, z=orderZ
    mapping[warpmIdx(0, y, orderZ)] = vtkIdx++;

  // Vertical edges
  for (int z = 1; z < orderZ; ++z)  // Edge 0-4: x=0, y=0
    mapping[warpmIdx(0, 0, z)] = vtkIdx++;
  for (int z = 1; z < orderZ; ++z)  // Edge 1-5: x=orderX, y=0
    mapping[warpmIdx(orderX, 0, z)] = vtkIdx++;
  for (int z = 1; z < orderZ; ++z)  // Edge 2-6: x=orderX, y=orderY
    mapping[warpmIdx(orderX, orderY, z)] = vtkIdx++;
  for (int z = 1; z < orderZ; ++z)  // Edge 3-7: x=0, y=orderY
    mapping[warpmIdx(0, orderY, z)] = vtkIdx++;

  // 3. Six faces (interior nodes only)
  // Face at z=0 (bottom)
  for (int y = 1; y < orderY; ++y)
    for (int x = 1; x < orderX; ++x)
      mapping[warpmIdx(x, y, 0)] = vtkIdx++;

  // Face at z=orderZ (top)
  for (int y = 1; y < orderY; ++y)
    for (int x = 1; x < orderX; ++x)
      mapping[warpmIdx(x, y, orderZ)] = vtkIdx++;

  // Face at y=0 (front)
  for (int z = 1; z < orderZ; ++z)
    for (int x = 1; x < orderX; ++x)
      mapping[warpmIdx(x, 0, z)] = vtkIdx++;

  // Face at y=orderY (back)
  for (int z = 1; z < orderZ; ++z)
    for (int x = 1; x < orderX; ++x)
      mapping[warpmIdx(x, orderY, z)] = vtkIdx++;

  // Face at x=0 (left)
  for (int z = 1; z < orderZ; ++z)
    for (int y = 1; y < orderY; ++y)
      mapping[warpmIdx(0, y, z)] = vtkIdx++;

  // Face at x=orderX (right)
  for (int z = 1; z < orderZ; ++z)
    for (int y = 1; y < orderY; ++y)
      mapping[warpmIdx(orderX, y, z)] = vtkIdx++;

  // 4. Interior volume nodes
  for (int z = 1; z < orderZ; ++z)
    for (int y = 1; y < orderY; ++y)
      for (int x = 1; x < orderX; ++x)
        mapping[warpmIdx(x, y, z)] = vtkIdx++;

  return mapping;
}

// Dispatcher: select 1D/2D/3D mapping based on orders vector size
static std::vector<int> PSBuildWarpmToVTKMapping(const std::vector<int>& orders)
{
  switch (orders.size())
  {
    case 1:
      return PSBuildWarpmToVTKMapping1D(orders[0]);
    case 2:
      return PSBuildWarpmToVTKMapping2D(orders[0], orders[1]);
    case 3:
      return PSBuildWarpmToVTKMapping3D(orders[0], orders[1], orders[2]);
    default:
      // Return empty mapping for unsupported dimensions
      return std::vector<int>();
  }
}

// Check if an HDF5 file has WARPM structure
static bool PSIsWarpmHDF5File(const char* fname)
{
  hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) return false;

  bool hasRequired = true;

  hid_t domainsGroup = H5Gopen(file, "/domains", H5P_DEFAULT);
  if (domainsGroup < 0) hasRequired = false;
  else H5Gclose(domainsGroup);

  hid_t varsGroup = H5Gopen(file, "/variables", H5P_DEFAULT);
  if (varsGroup < 0) hasRequired = false;
  else H5Gclose(varsGroup);

  hid_t timeGroup = H5Gopen(file, "/timeData", H5P_DEFAULT);
  if (timeGroup < 0) hasRequired = false;
  else H5Gclose(timeGroup);

  H5Fclose(file);
  return hasRequired;
}

// Check if an HDF5 file has phase space variables
static bool PSHasPhaseSpaceVariables(const char* fname)
{
  hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) return false;

  bool hasPhaseSpace = false;

  hid_t domainsGroup = H5Gopen(file, "/domains", H5P_DEFAULT);
  hid_t varsGroup = H5Gopen(file, "/variables", H5P_DEFAULT);

  if (domainsGroup >= 0 && varsGroup >= 0)
  {
    // Iterate over variables to find one on a phase space domain
    std::vector<std::string> varNames;
    H5Literate(varsGroup, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr,
      [](hid_t, const char* name, const H5L_info_t*, void* opdata) -> herr_t {
        static_cast<std::vector<std::string>*>(opdata)->push_back(name);
        return 0;
      }, &varNames);

    for (const auto& varName : varNames)
    {
      if (hasPhaseSpace) break;

      hid_t varGroup = H5Gopen(varsGroup, varName.c_str(), H5P_DEFAULT);
      if (varGroup >= 0)
      {
        std::string domainName;
        if (PSReadStringAttribute(varGroup, "OnDomain", domainName))
        {
          hid_t domGroup = H5Gopen(domainsGroup, domainName.c_str(), H5P_DEFAULT);
          if (domGroup >= 0)
          {
            int numPhys = 0, numVel = 0;
            std::vector<std::string> coordNames;
            if (PSDetectPhaseSpaceDomain(domGroup, numPhys, numVel, coordNames))
            {
              hasPhaseSpace = true;
            }
            H5Gclose(domGroup);
          }
        }
        H5Gclose(varGroup);
      }
    }
  }

  if (varsGroup >= 0) H5Gclose(varsGroup);
  if (domainsGroup >= 0) H5Gclose(domainsGroup);
  H5Fclose(file);

  return hasPhaseSpace;
}

//----------------------------------------------------------------------------
vtkWARPMPhaseSpaceReader::vtkWARPMPhaseSpaceReader()
  : PhysicalNodeIndex(0)
  , HasPhaseSpaceVariables(false)
  , NumPhysicalDims(0)
  , NumVelocityDims(0)
{
  this->PhysicalSliceIndices[0] = 0;
  this->PhysicalSliceIndices[1] = 0;
  this->PhysicalSliceIndices[2] = 0;
  // Override base class: set 3 output ports
  // Port 0: Physical space mesh
  // Port 1: Velocity space mesh
  // Port 2: Probe location (single point)
  this->SetNumberOfOutputPorts(3);
}

//----------------------------------------------------------------------------
vtkWARPMPhaseSpaceReader::~vtkWARPMPhaseSpaceReader()
{
}

//----------------------------------------------------------------------------
int vtkWARPMPhaseSpaceReader::FillOutputPortInformation(int port, vtkInformation* info)
{
  if (port == 2)
  {
    // Port 2: Probe location (single point showing slice position in physical space)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  }
  else
  {
    // Ports 0 and 1: vtkUnstructuredGrid
    // Port 0: Physical space mesh
    // Port 1: Velocity space mesh
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkWARPMPhaseSpaceReader::CanReadFile(const char* fname)
{
  if (!fname || strlen(fname) == 0)
  {
    return 0;
  }

  std::filesystem::path filePath(fname);
  std::string ext = filePath.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  std::string h5FileToCheck;

  if (ext == ".h5" || ext == ".hdf5")
  {
    // Direct HDF5 file
    if (!PSIsWarpmHDF5File(fname)) return 0;
    h5FileToCheck = fname;
  }
  else if (ext == ".warpm")
  {
    // Metadata file - parse to get first H5 and validate it
    std::filesystem::path baseDir = filePath.parent_path();
    if (baseDir.empty()) baseDir = ".";

    std::ifstream file(fname);
    if (!file.is_open()) return 0;

    std::string line;
    while (std::getline(file, line))
    {
      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);
      if (line.empty() || line[0] == '#') continue;

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
                  h5FileToCheck = entry.path().string();
                  break;
                }
              }
            }
          }
        }
      }
      else
      {
        std::filesystem::path framePath = baseDir / line;
        if (std::filesystem::exists(framePath))
        {
          h5FileToCheck = framePath.string();
        }
      }
      if (!h5FileToCheck.empty()) break;
    }
  }
  else
  {
    return 0; // Unknown extension
  }

  if (h5FileToCheck.empty()) return 0;

  // Key difference from base class: only return true if phase space variables exist
  return PSHasPhaseSpaceVariables(h5FileToCheck.c_str()) ? 1 : 0;
}

//----------------------------------------------------------------------------
void vtkWARPMPhaseSpaceReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "PhysicalSliceIndices: " << this->PhysicalSliceIndices[0]
     << ", " << this->PhysicalSliceIndices[1]
     << ", " << this->PhysicalSliceIndices[2] << "\n";
  os << indent << "PhysicalNodeIndex: " << this->PhysicalNodeIndex << "\n";
  os << indent << "HasPhaseSpaceVariables: " << this->HasPhaseSpaceVariables << "\n";
}

//----------------------------------------------------------------------------
int vtkWARPMPhaseSpaceReader::RequestInformation(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // Call base class to parse file, populate FrameFiles, TimeValues, and enumerate variables
  if (!this->Superclass::RequestInformation(request, inputVector, outputVector))
  {
    return 0;
  }

  // Reset phase space detection
  this->HasPhaseSpaceVariables = false;
  this->NumPhysicalDims = 0;
  this->NumVelocityDims = 0;
  this->PhaseSpaceDomainName.clear();

  // Detect phase space domains
  if (!this->FrameFiles.empty())
  {
    hid_t file = H5Fopen(this->FrameFiles[0].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file >= 0)
    {
      hid_t domainsGroup = H5Gopen(file, "/domains", H5P_DEFAULT);
      hid_t varsGroup = H5Gopen(file, "/variables", H5P_DEFAULT);

      if (domainsGroup >= 0 && varsGroup >= 0)
      {
        std::vector<std::string> varNames;
        H5Literate(varsGroup, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr,
          [](hid_t, const char* name, const H5L_info_t*, void* opdata) -> herr_t {
            static_cast<std::vector<std::string>*>(opdata)->push_back(name);
            return 0;
          }, &varNames);

        // First pass: detect phase space domain and add phase space variables to selection
        for (const auto& name : varNames)
        {
          hid_t varGroup = H5Gopen(varsGroup, name.c_str(), H5P_DEFAULT);
          if (varGroup >= 0)
          {
            std::string domainName;
            if (PSReadStringAttribute(varGroup, "OnDomain", domainName))
            {
              hid_t domGroup = H5Gopen(domainsGroup, domainName.c_str(), H5P_DEFAULT);
              if (domGroup >= 0)
              {
                int numPhys = 0, numVel = 0;
                std::vector<std::string> coordNames;
                if (PSDetectPhaseSpaceDomain(domGroup, numPhys, numVel, coordNames))
                {
                  // Found a phase space variable - record the domain info
                  if (!this->HasPhaseSpaceVariables)
                  {
                    this->HasPhaseSpaceVariables = true;
                    this->NumPhysicalDims = numPhys;
                    this->NumVelocityDims = numVel;
                    this->PhaseSpaceDomainName = domainName;
                  }

                  // Add phase space variable to selection (base class only added physical vars)
                  if (!this->PointDataArraySelection->ArrayExists(name.c_str()))
                  {
                    this->PointDataArraySelection->AddArray(name.c_str());
                    this->PointDataArraySelection->EnableArray(name.c_str());
                  }
                }
                H5Gclose(domGroup);
              }
            }
            H5Gclose(varGroup);
          }
        }

        H5Gclose(varsGroup);
      }
      if (domainsGroup >= 0) H5Gclose(domainsGroup);
      H5Fclose(file);
    }
  }

  // Set time information on BOTH output ports (critical for multi-port time animation)
  if (!this->TimeValues.empty())
  {
    for (int port = 0; port < this->GetNumberOfOutputPorts(); ++port)
    {
      vtkInformation* outInfo = outputVector->GetInformationObject(port);
      if (outInfo)
      {
        outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                     this->TimeValues.data(),
                     static_cast<int>(this->TimeValues.size()));

        double timeRange[2] = { this->TimeValues.front(), this->TimeValues.back() };
        outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
      }
    }
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkWARPMPhaseSpaceReader::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // Debug: track RequestData calls
  static int callCount = 0;
  callCount++;
  {
    std::ofstream debugFile("/tmp/_warpm_debug.txt", std::ios::app);
    debugFile << "=== RequestData call #" << callCount << " ===\n";
  }

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid* output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!output)
  {
    vtkErrorMacro("Could not get output port 0");
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
  hid_t domainsGroup = H5Gopen(file, "/domains", H5P_DEFAULT);
  hid_t varsGroup = H5Gopen(file, "/variables", H5P_DEFAULT);
  if (varsGroup < 0)
  {
    vtkErrorMacro("Could not open variables group");
    if (domainsGroup >= 0) H5Gclose(domainsGroup);
    H5Fclose(file);
    return 0;
  }

  // Collect all variable names
  std::vector<std::string> allVarNames;
  H5Literate(varsGroup, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr,
    [](hid_t, const char* name, const H5L_info_t*, void* opdata) -> herr_t {
      static_cast<std::vector<std::string>*>(opdata)->push_back(name);
      return 0;
    }, &allVarNames);

  // Separate enabled variables into physical-space and phase-space
  std::vector<std::string> physicalVars;
  std::vector<std::string> phaseSpaceVars;
  std::string physicalDomainName; // Will be set from the first physical variable found
  std::string phaseSpaceDomainName; // Will be set from the first phase space variable found

  for (const auto& name : allVarNames)
  {
    if (!this->PointDataArraySelection->ArrayIsEnabled(name.c_str()))
    {
      continue;
    }

    // Check which domain this variable uses (based on coordinate names, not domain name)
    bool isPhaseSpace = false;
    std::string varDomainName;
    hid_t varGroup = H5Gopen(varsGroup, name.c_str(), H5P_DEFAULT);
    if (varGroup >= 0)
    {
      if (PSReadStringAttribute(varGroup, "OnDomain", varDomainName) && domainsGroup >= 0)
      {
        hid_t domGroup = H5Gopen(domainsGroup, varDomainName.c_str(), H5P_DEFAULT);
        if (domGroup >= 0)
        {
          int numPhys = 0, numVel = 0;
          std::vector<std::string> coordNames;
          isPhaseSpace = PSDetectPhaseSpaceDomain(domGroup, numPhys, numVel, coordNames);
          H5Gclose(domGroup);
        }
      }
      H5Gclose(varGroup);
    }

    if (isPhaseSpace)
    {
      if (phaseSpaceDomainName.empty())
      {
        // First phase space variable - use its domain
        phaseSpaceDomainName = varDomainName;
        phaseSpaceVars.push_back(name);
      }
      else if (varDomainName == phaseSpaceDomainName)
      {
        // Same domain as first phase space variable
        phaseSpaceVars.push_back(name);
      }
      else
      {
        // Different domain - warn and skip
        vtkWarningMacro("Skipping variable '" << name << "' on domain '"
          << varDomainName << "' (differs from '" << phaseSpaceDomainName << "')");
      }
    }
    else
    {
      if (physicalDomainName.empty())
      {
        // First physical variable - use its domain
        physicalDomainName = varDomainName;
        physicalVars.push_back(name);
      }
      else if (varDomainName == physicalDomainName)
      {
        // Same domain as first physical variable
        physicalVars.push_back(name);
      }
      else
      {
        // Different domain - warn and skip
        vtkWarningMacro("Skipping variable '" << name << "' on domain '"
          << varDomainName << "' (differs from '" << physicalDomainName << "')");
      }
    }
  }

  // Debug output to trace variable detection
  {
    std::ofstream debugFile("/tmp/_warpm_debug.txt");
    debugFile << "Physical vars (" << physicalVars.size() << "): ";
    for (const auto& v : physicalVars) debugFile << v << " ";
    debugFile << "\n";
    debugFile << "Physical domain: " << physicalDomainName << "\n";
    debugFile << "Phase space vars (" << phaseSpaceVars.size() << "): ";
    for (const auto& v : phaseSpaceVars) debugFile << v << " ";
    debugFile << "\n";
    debugFile << "Phase space domain: " << phaseSpaceDomainName << "\n";
  }

  if (physicalVars.empty() && phaseSpaceVars.empty())
  {
    vtkWarningMacro("No variables enabled");
    if (domainsGroup >= 0) H5Gclose(domainsGroup);
    H5Gclose(varsGroup);
    H5Fclose(file);
    return 1;
  }

  // ============================================================
  // PORT 0: Physical space variables
  // ============================================================
  if (!physicalVars.empty())
  {
    std::string domainPath = "/domains/" + physicalDomainName;
    hid_t domain = H5Gopen(file, domainPath.c_str(), H5P_DEFAULT);
    if (domain < 0)
    {
      vtkErrorMacro("Could not open domain group: " << physicalDomainName);
      if (domainsGroup >= 0) H5Gclose(domainsGroup);
      H5Gclose(varsGroup);
      H5Fclose(file);
      return 0;
    }

    int ndims = 0;
    std::vector<int> numCells;
    std::vector<int> startIndices;
    std::vector<std::string> coordExprs;

    PSReadIntAttribute(domain, "ndims", ndims);
    PSReadIntArrayAttribute(domain, "numCells", numCells);
    PSReadIntArrayAttribute(domain, "startIndices", startIndices);
    PSReadStringArrayAttribute(domain, "VertexCoordinateExpressions", coordExprs);

    H5Gclose(domain);

    if (ndims < 1 || numCells.empty() || coordExprs.empty())
    {
      vtkErrorMacro("Invalid domain dimensions");
      if (domainsGroup >= 0) H5Gclose(domainsGroup);
      H5Gclose(varsGroup);
      H5Fclose(file);
      return 0;
    }

    // Limit to 3D max
    int numPhysDims = std::min(ndims, 3);

    // Default startIndices to 0 if not present
    if (startIndices.size() < static_cast<size_t>(numPhysDims))
    {
      startIndices.resize(numPhysDims, 0);
    }

    // Get element order from first variable
    hid_t firstVarGroup = H5Gopen(varsGroup, physicalVars[0].c_str(), H5P_DEFAULT);
    std::vector<int> elementOrder;
    int entriesPerElement = 0;

    PSReadIntArrayAttribute(firstVarGroup, "ElementOrder", elementOrder);
    PSReadIntAttribute(firstVarGroup, "EntriesPerElement", entriesPerElement);

    hid_t firstDataset = H5Dopen(firstVarGroup, "data", H5P_DEFAULT);
    hid_t firstDataspace = H5Dget_space(firstDataset);
    int dataNdims = H5Sget_simple_extent_ndims(firstDataspace);
    std::vector<hsize_t> dataDims(dataNdims);
    H5Sget_simple_extent_dims(firstDataspace, dataDims.data(), nullptr);
    H5Sclose(firstDataspace);
    H5Dclose(firstDataset);
    H5Gclose(firstVarGroup);

    // Build dimension arrays for N-dimensional physical mesh
    std::vector<int> physNumCellsVec(numPhysDims);
    std::vector<int> physOrders(numPhysDims);
    std::vector<std::vector<double>> physVertices(numPhysDims);
    std::vector<std::vector<double>> physGLL(numPhysDims);

    for (int d = 0; d < numPhysDims; ++d)
    {
      physNumCellsVec[d] = numCells[d];
      physOrders[d] = (d < static_cast<int>(elementOrder.size())) ? elementOrder[d] : 2;
      physVertices[d] = PSComputeVertexPositions(coordExprs[d], startIndices[d], numCells[d]);
      physGLL[d] = PSGetGLLNodes(physOrders[d]);
    }

    int nodesPerElement = PSComputeTotalNodes(physOrders);
    int totalCells = 1;
    for (int nc : physNumCellsVec) totalCells *= nc;
    int totalNodes = totalCells * nodesPerElement;

    // Build nodes-per-dimension array for index decomposition
    std::vector<int> physNodesPerDim(numPhysDims);
    for (int d = 0; d < numPhysDims; ++d)
    {
      physNodesPerDim[d] = physOrders[d] + 1;
    }

    // Check if we can use cached geometry
    bool canUseCache = this->GeometryCached &&
                       this->CachedNodesPerElement == nodesPerElement &&
                       this->CachedNdims == numPhysDims &&
                       this->CachedDataDims == physNumCellsVec;

    if (!canUseCache)
    {
      this->CachedWarpmToVTK = PSBuildWarpmToVTKMapping(physOrders);

      this->CachedPoints = vtkSmartPointer<vtkPoints>::New();
      this->CachedPoints->SetNumberOfPoints(totalNodes);

      // Iterate over all physical cells using flat index
      for (int cellIdx = 0; cellIdx < totalCells; ++cellIdx)
      {
        // Decompose cell index into per-dimension indices
        std::vector<int> cellIndices = PSDecomposeIndex(cellIdx, physNumCellsVec);

        // Get cell bounds from vertex positions for each dimension
        std::vector<double> cellMin(numPhysDims), cellMax(numPhysDims);
        for (int d = 0; d < numPhysDims; ++d)
        {
          cellMin[d] = physVertices[d][cellIndices[d]];
          cellMax[d] = physVertices[d][cellIndices[d] + 1];
        }

        int elemStartIdx = cellIdx * nodesPerElement;

        // Iterate over all nodes within this cell
        for (int nodeIdx = 0; nodeIdx < nodesPerElement; ++nodeIdx)
        {
          // Decompose node index into per-dimension indices
          std::vector<int> nodeIndices = PSDecomposeIndex(nodeIdx, physNodesPerDim);

          // Compute physical position (VTK always uses 3D coordinates)
          double pos[3] = {0.0, 0.0, 0.0};
          for (int d = 0; d < numPhysDims; ++d)
          {
            double xi = physGLL[d][nodeIndices[d]];
            pos[d] = cellMin[d] + (xi + 1.0) * 0.5 * (cellMax[d] - cellMin[d]);
          }

          int vtkLocalIdx = this->CachedWarpmToVTK[nodeIdx];
          int vtkPointIdx = elemStartIdx + vtkLocalIdx;
          this->CachedPoints->SetPoint(vtkPointIdx, pos[0], pos[1], pos[2]);
        }
      }

      this->CachedCells = vtkSmartPointer<vtkCellArray>::New();
      for (int cellIdx = 0; cellIdx < totalCells; ++cellIdx)
      {
        int elemStartIdx = cellIdx * nodesPerElement;
        std::vector<vtkIdType> cellPts(nodesPerElement);
        for (int i = 0; i < nodesPerElement; ++i)
        {
          cellPts[i] = elemStartIdx + i;
        }
        this->CachedCells->InsertNextCell(nodesPerElement, cellPts.data());
      }

      this->CachedNdims = numPhysDims;
      this->CachedDataDims = physNumCellsVec;
      this->CachedNodesPerElement = nodesPerElement;
      this->GeometryCached = true;
    }

    output->SetPoints(this->CachedPoints);
    output->SetCells(PSGetLagrangeCellType(numPhysDims), this->CachedCells);

    // Load data for each physical-space variable
    for (const auto& varName : physicalVars)
    {
      hid_t varGroup = H5Gopen(varsGroup, varName.c_str(), H5P_DEFAULT);
      if (varGroup < 0)
      {
        vtkWarningMacro("Could not open variable: " << varName);
        continue;
      }

      std::vector<std::string> componentNames;
      PSReadStringArrayAttribute(varGroup, "ComponentNames", componentNames);

      hid_t dataset = H5Dopen(varGroup, "data", H5P_DEFAULT);
      hid_t dataspace = H5Dget_space(dataset);

      int varNdims = H5Sget_simple_extent_ndims(dataspace);
      std::vector<hsize_t> varDims(varNdims);
      H5Sget_simple_extent_dims(dataspace, varDims.data(), nullptr);

      // Last dimension is components, second-to-last is nodesPerElement
      int varNumComponents = static_cast<int>(varDims[varNdims - 1]);

      // Read all data
      hsize_t totalSize = 1;
      for (int i = 0; i < varNdims; ++i) totalSize *= varDims[i];
      std::vector<double> data(totalSize);
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

      H5Sclose(dataspace);
      H5Dclose(dataset);
      H5Gclose(varGroup);

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

      // Generalized N-dimensional data loading loop
      for (int cellIdx = 0; cellIdx < totalCells; ++cellIdx)
      {
        // Decompose cell index into per-dimension indices
        std::vector<int> cellIndices = PSDecomposeIndex(cellIdx, physNumCellsVec);

        int elemStartIdx = cellIdx * nodesPerElement;

        for (int nodeIdx = 0; nodeIdx < nodesPerElement; ++nodeIdx)
        {
          int vtkLocalIdx = this->CachedWarpmToVTK[nodeIdx];
          int vtkPointIdx = elemStartIdx + vtkLocalIdx;

          // Compute flat data index (row-major order in HDF5 data)
          int dataCellIdx = PSComputeFlatIndex(cellIndices, physNumCellsVec);
          int dataIdx = (dataCellIdx * nodesPerElement + nodeIdx) * varNumComponents;

          for (int c = 0; c < varNumComponents; ++c)
          {
            if (dataIdx + c < static_cast<int>(data.size()))
            {
              fieldArrays[c]->SetValue(vtkPointIdx, data[dataIdx + c]);
            }
          }
        }
      }

      for (auto& arr : fieldArrays)
      {
        output->GetPointData()->AddArray(arr);
      }
    }
  }

  // ============================================================
  // PORT 1: Phase space variables (velocity space mesh)
  // ============================================================
  if (!phaseSpaceVars.empty())
  {
    vtkInformation* phaseOutInfo = outputVector->GetInformationObject(1);
    vtkUnstructuredGrid* phaseOutput = vtkUnstructuredGrid::SafeDownCast(
      phaseOutInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Debug: trace velocity mesh conditions
    {
      std::ofstream debugFile("/tmp/_warpm_debug.txt", std::ios::app);
      debugFile << "--- Velocity mesh building ---\n";
      debugFile << "phaseOutInfo: " << (phaseOutInfo ? "valid" : "NULL") << "\n";
      debugFile << "phaseOutput: " << (phaseOutput ? "valid" : "NULL") << "\n";
      debugFile << "domainsGroup: " << domainsGroup << "\n";
      debugFile << "phaseSpaceDomainName: '" << phaseSpaceDomainName << "'\n";
    }

    if (phaseOutput && domainsGroup >= 0 && !phaseSpaceDomainName.empty())
    {
      hid_t phaseDomGroup = H5Gopen(domainsGroup, phaseSpaceDomainName.c_str(), H5P_DEFAULT);
      if (phaseDomGroup >= 0)
      {
        std::vector<int> phaseNumCells;
        std::vector<int> phaseStartIndices;
        std::vector<std::string> phaseCoordExprs;
        PSReadIntArrayAttribute(phaseDomGroup, "numCells", phaseNumCells);
        PSReadIntArrayAttribute(phaseDomGroup, "startIndices", phaseStartIndices);
        PSReadStringArrayAttribute(phaseDomGroup, "VertexCoordinateExpressions", phaseCoordExprs);

        std::vector<std::string> coordNames;
        PSReadStringArrayAttribute(phaseDomGroup, "CoordinateNames", coordNames);

        // Default startIndices if not present
        if (phaseStartIndices.size() < coordNames.size())
        {
          phaseStartIndices.resize(coordNames.size(), 0);
        }

        std::vector<int> velDimIndices;
        std::vector<int> physDimIndices;
        for (size_t i = 0; i < coordNames.size(); ++i)
        {
          if (!coordNames[i].empty() && (coordNames[i][0] == 'v' || coordNames[i][0] == 'w'))
          {
            velDimIndices.push_back(static_cast<int>(i));
          }
          else
          {
            physDimIndices.push_back(static_cast<int>(i));
          }
        }

        // Debug: trace velocity dimension detection
        {
          std::ofstream debugFile("/tmp/_warpm_debug.txt", std::ios::app);
          debugFile << "--- Phase space domain analysis ---\n";
          debugFile << "Total dimensions: " << coordNames.size() << "\n";
          debugFile << "Physical dimensions: " << physDimIndices.size() << " (indices: ";
          for (auto c : physDimIndices) debugFile << c << " ";
          debugFile << ")\n";
          debugFile << "Velocity dimensions: " << velDimIndices.size() << " (indices: ";
          for (auto c : velDimIndices) debugFile << c << " ";
          debugFile << ")\n";
          debugFile << "coordNames: ";
          for (const auto& c : coordNames) debugFile << c << " ";
          debugFile << "\n";
          debugFile << "phaseNumCells: ";
          for (auto c : phaseNumCells) debugFile << c << " ";
          debugFile << "\n";
          debugFile << "phaseStartIndices: ";
          for (auto c : phaseStartIndices) debugFile << c << " ";
          debugFile << "\n";
          debugFile << "Velocity coordinate expressions:\n";
          for (size_t i = 0; i < velDimIndices.size(); ++i)
          {
            int idx = velDimIndices[i];
            if (idx < static_cast<int>(phaseCoordExprs.size()))
            {
              debugFile << "  dim " << i << " (idx " << idx << ", " << coordNames[idx] << "): "
                        << phaseCoordExprs[idx] << "\n";
            }
          }
          debugFile << "Physical coordinate expressions:\n";
          for (size_t i = 0; i < physDimIndices.size(); ++i)
          {
            int idx = physDimIndices[i];
            if (idx < static_cast<int>(phaseCoordExprs.size()))
            {
              debugFile << "  dim " << i << " (idx " << idx << ", " << coordNames[idx] << "): "
                        << phaseCoordExprs[idx] << "\n";
            }
          }
        }

        int numPhysDims = static_cast<int>(physDimIndices.size());
        int numVelDims = static_cast<int>(velDimIndices.size());

        // Validate PhysicalSliceIndices for dimensions beyond actual physical dimensions
        // For 1D physical data, indices [1] and [2] should be 0
        // For 2D physical data, index [2] should be 0
        for (int d = numPhysDims; d < 3; ++d)
        {
          if (this->PhysicalSliceIndices[d] != 0)
          {
            vtkWarningMacro("PhysicalSliceIndices[" << d << "] = "
              << this->PhysicalSliceIndices[d]
              << " is invalid for " << numPhysDims << "D physical data (only "
              << numPhysDims << " physical dimension" << (numPhysDims > 1 ? "s" : "")
              << "), value ignored (treated as 0)");
          }
        }

        if (numVelDims >= 1 && numVelDims <= 3)
        {
          // Build velocity dimension arrays
          std::vector<int> velNumCellsVec(numVelDims);
          std::vector<int> velOrders(numVelDims);
          std::vector<std::vector<double>> velVertices(numVelDims);
          std::vector<std::vector<double>> velGLL(numVelDims);

          // Debug: trace velocity mesh computation
          {
            std::ofstream debugFile("/tmp/_warpm_debug.txt", std::ios::app);
            debugFile << "--- Computing " << numVelDims << "D velocity mesh (GENERALIZED) ---\n";
            debugFile << "velNumCells: [";
            for (int d = 0; d < numVelDims; ++d)
            {
              int idx = velDimIndices[d];
              debugFile << phaseNumCells[idx] << (d < numVelDims - 1 ? ", " : "");
            }
            debugFile << "]\n";
            debugFile << "Computing vertex positions for each velocity dimension:\n";
            for (int d = 0; d < numVelDims; ++d)
            {
              int idx = velDimIndices[d];
              debugFile << "  vel dim " << d << ": startIndex=" << phaseStartIndices[idx]
                        << ", numCells=" << phaseNumCells[idx] << "\n";
            }
          }

          // Get element orders from first phase space variable
          hid_t firstPhaseVar = H5Gopen(varsGroup, phaseSpaceVars[0].c_str(), H5P_DEFAULT);
          std::vector<int> phaseElemOrder;
          PSReadIntArrayAttribute(firstPhaseVar, "ElementOrder", phaseElemOrder);
          H5Gclose(firstPhaseVar);

          // Compute vertex positions and setup for each velocity dimension
          for (int d = 0; d < numVelDims; ++d)
          {
            int idx = velDimIndices[d];
            velNumCellsVec[d] = phaseNumCells[idx];
            velOrders[d] = (idx < static_cast<int>(phaseElemOrder.size()))
                           ? phaseElemOrder[idx] : 2;
            velVertices[d] = PSComputeVertexPositions(
              phaseCoordExprs[idx], phaseStartIndices[idx], phaseNumCells[idx]);
            velGLL[d] = PSGetGLLNodes(velOrders[d]);
          }

          // Debug: check vertex results
          {
            std::ofstream debugFile("/tmp/_warpm_debug.txt", std::ios::app);
            for (int d = 0; d < numVelDims; ++d)
            {
              debugFile << "velVertices[" << d << "] (" << velVertices[d].size() << " values): ";
              for (double v : velVertices[d]) debugFile << v << " ";
              debugFile << "\n";
            }
          }

          // Compute total cells and nodes
          int velNodesPerElement = PSComputeTotalNodes(velOrders);
          int velTotalCells = 1;
          for (int nc : velNumCellsVec) velTotalCells *= nc;
          int velTotalNodes = velTotalCells * velNodesPerElement;

          // Build node ordering mapping
          auto velWarpmToVTK = PSBuildWarpmToVTKMapping(velOrders);

          // Build nodes-per-dimension array for index decomposition
          std::vector<int> velNodesPerDim(numVelDims);
          for (int d = 0; d < numVelDims; ++d)
          {
            velNodesPerDim[d] = velOrders[d] + 1;
          }

          // Create velocity mesh points
          vtkNew<vtkPoints> velPoints;
          velPoints->SetNumberOfPoints(velTotalNodes);

          // Iterate over all velocity cells using flat index
          for (int cellIdx = 0; cellIdx < velTotalCells; ++cellIdx)
          {
            // Decompose cell index into per-dimension indices
            std::vector<int> cellIndices = PSDecomposeIndex(cellIdx, velNumCellsVec);

            // Get cell bounds from vertex positions for each dimension
            std::vector<double> cellMin(numVelDims), cellMax(numVelDims);
            for (int d = 0; d < numVelDims; ++d)
            {
              cellMin[d] = velVertices[d][cellIndices[d]];
              cellMax[d] = velVertices[d][cellIndices[d] + 1];
            }

            int elemStartIdx = cellIdx * velNodesPerElement;

            // Iterate over all nodes within this cell
            for (int nodeIdx = 0; nodeIdx < velNodesPerElement; ++nodeIdx)
            {
              // Decompose node index into per-dimension indices
              std::vector<int> nodeIndices = PSDecomposeIndex(nodeIdx, velNodesPerDim);

              // Compute physical position (VTK always uses 3D coordinates)
              double pos[3] = {0.0, 0.0, 0.0};
              for (int d = 0; d < numVelDims; ++d)
              {
                double xi = velGLL[d][nodeIndices[d]];
                pos[d] = cellMin[d] + (xi + 1.0) * 0.5 * (cellMax[d] - cellMin[d]);
              }

              int vtkLocalIdx = velWarpmToVTK[nodeIdx];
              int vtkPointIdx = elemStartIdx + vtkLocalIdx;
              velPoints->SetPoint(vtkPointIdx, pos[0], pos[1], pos[2]);
            }
          }

          // Create velocity mesh cells
          vtkNew<vtkCellArray> velCells;
          for (int cellIdx = 0; cellIdx < velTotalCells; ++cellIdx)
          {
            int elemStartIdx = cellIdx * velNodesPerElement;
            std::vector<vtkIdType> cellPts(velNodesPerElement);
            for (int i = 0; i < velNodesPerElement; ++i)
            {
              cellPts[i] = elemStartIdx + i;
            }
            velCells->InsertNextCell(velNodesPerElement, cellPts.data());
          }

          // Set mesh with appropriate cell type
          phaseOutput->SetPoints(velPoints);
          phaseOutput->SetCells(PSGetLagrangeCellType(numVelDims), velCells);

          // Debug: verify velocity mesh was set correctly
          {
            std::ofstream debugFile("/tmp/_warpm_debug.txt", std::ios::app);
            debugFile << "--- After " << numVelDims << "D velocity mesh setup ---\n";
            debugFile << "velTotalNodes=" << velTotalNodes << ", velTotalCells=" << velTotalCells << "\n";
            debugFile << "phaseOutput points: " << phaseOutput->GetNumberOfPoints() << "\n";
            debugFile << "phaseOutput cells: " << phaseOutput->GetNumberOfCells() << "\n";
            debugFile << "Cell type: " << PSGetLagrangeCellType(numVelDims) << "\n";
            double bounds[6];
            phaseOutput->GetBounds(bounds);
            debugFile << "phaseOutput bounds: X=[" << bounds[0] << ", " << bounds[1]
                      << "], Y=[" << bounds[2] << ", " << bounds[3]
                      << "], Z=[" << bounds[4] << ", " << bounds[5] << "]\n";
          }

          // ============================================================
          // PORT 2: Probe location (single point in physical space)
          // ============================================================
          vtkInformation* probeOutInfo = outputVector->GetInformationObject(2);
          vtkPolyData* probeOutput = vtkPolyData::SafeDownCast(
            probeOutInfo->Get(vtkDataObject::DATA_OBJECT()));

          int numPhysDims = static_cast<int>(physDimIndices.size());
          if (probeOutput && numPhysDims >= 1 && numPhysDims <= 3)
          {
            // Get physical element orders
            hid_t probePhaseVar = H5Gopen(varsGroup, phaseSpaceVars[0].c_str(), H5P_DEFAULT);
            std::vector<int> probeElemOrder;
            PSReadIntArrayAttribute(probePhaseVar, "ElementOrder", probeElemOrder);
            H5Gclose(probePhaseVar);

            // Build physical dimension arrays
            std::vector<std::vector<double>> physVertices(numPhysDims);
            std::vector<int> physOrders(numPhysDims);
            std::vector<std::vector<double>> physGLL(numPhysDims);
            std::vector<int> sliceIndices(numPhysDims);

            for (int d = 0; d < numPhysDims; ++d)
            {
              int idx = physDimIndices[d];
              physOrders[d] = (idx < static_cast<int>(probeElemOrder.size()))
                              ? probeElemOrder[idx] : 2;
              physVertices[d] = PSComputeVertexPositions(
                phaseCoordExprs[idx], phaseStartIndices[idx], phaseNumCells[idx]);
              physGLL[d] = PSGetGLLNodes(physOrders[d]);

              // Get slice index with bounds checking
              int maxSlice = phaseNumCells[idx] - 1;
              sliceIndices[d] = std::max(0, std::min(this->PhysicalSliceIndices[d], maxSlice));
            }

            // Compute total nodes per physical cell and decompose node index
            int physNodesPerCell = PSComputeTotalNodes(physOrders);
            std::vector<int> physNodesPerDim(numPhysDims);
            for (int d = 0; d < numPhysDims; ++d)
            {
              physNodesPerDim[d] = physOrders[d] + 1;
            }

            int physNodeIdx = this->PhysicalNodeIndex;
            int maxNodeIdx = physNodesPerCell - 1;
            physNodeIdx = std::max(0, std::min(physNodeIdx, maxNodeIdx));

            // Decompose node index into per-dimension indices
            std::vector<int> nodeIndices = PSDecomposeIndex(physNodeIdx, physNodesPerDim);

            // Compute probe position (VTK always uses 3D coordinates)
            double probePos[3] = {0.0, 0.0, 0.0};
            for (int d = 0; d < numPhysDims; ++d)
            {
              double cellMin = physVertices[d][sliceIndices[d]];
              double cellMax = physVertices[d][sliceIndices[d] + 1];
              double xi = physGLL[d][nodeIndices[d]];
              probePos[d] = cellMin + (xi + 1.0) * 0.5 * (cellMax - cellMin);
            }

            // Create probe point
            vtkNew<vtkPoints> probePoints;
            probePoints->InsertNextPoint(probePos[0], probePos[1], probePos[2]);

            vtkNew<vtkCellArray> probeVerts;
            vtkIdType ptId = 0;
            probeVerts->InsertNextCell(1, &ptId);

            probeOutput->SetPoints(probePoints);
            probeOutput->SetVerts(probeVerts);
          }

          // Load phase space variable data with slicing
          for (const auto& varName : phaseSpaceVars)
          {
            hid_t varGroup = H5Gopen(varsGroup, varName.c_str(), H5P_DEFAULT);
            if (varGroup < 0) continue;

            std::vector<std::string> componentNames;
            PSReadStringArrayAttribute(varGroup, "ComponentNames", componentNames);

            hid_t dataset = H5Dopen(varGroup, "data", H5P_DEFAULT);
            hid_t dataspace = H5Dget_space(dataset);

            int varNdims = H5Sget_simple_extent_ndims(dataspace);
            std::vector<hsize_t> varDims(varNdims);
            H5Sget_simple_extent_dims(dataspace, varDims.data(), nullptr);

            int numComponents = static_cast<int>(varDims[varNdims - 1]);
            int nodesPerElem = static_cast<int>(varDims[varNdims - 2]);

            std::vector<hsize_t> start(varNdims, 0);
            std::vector<hsize_t> count(varNdims);

            // Set hyperslab for physical dimensions (slice to single cell)
            for (size_t i = 0; i < physDimIndices.size() && i < 3; ++i)
            {
              int idx = physDimIndices[i];
              int sliceVal = this->PhysicalSliceIndices[i];
              int maxVal = static_cast<int>(varDims[idx]) - 1;
              if (sliceVal < 0 || sliceVal > maxVal)
              {
                vtkWarningMacro("PhysicalSliceIndices[" << i << "] = " << sliceVal
                  << " is out of range [0, " << maxVal << "], clamping");
                sliceVal = std::min(sliceVal, maxVal);
                sliceVal = std::max(sliceVal, 0);
              }
              start[idx] = sliceVal;
              count[idx] = 1;
            }

            for (size_t i = 0; i < velDimIndices.size(); ++i)
            {
              int idx = velDimIndices[i];
              start[idx] = 0;
              count[idx] = varDims[idx];
            }

            count[varNdims - 2] = nodesPerElem;
            count[varNdims - 1] = numComponents;

            hsize_t memSize = 1;
            for (auto c : count) memSize *= c;
            hid_t memspace = H5Screate_simple(1, &memSize, nullptr);

            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);

            std::vector<double> slicedData(memSize);
            H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, slicedData.data());

            H5Sclose(memspace);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            H5Gclose(varGroup);

            std::vector<int> phaseElemOrderRead;
            hid_t varGroupReopen = H5Gopen(varsGroup, varName.c_str(), H5P_DEFAULT);
            PSReadIntArrayAttribute(varGroupReopen, "ElementOrder", phaseElemOrderRead);
            H5Gclose(varGroupReopen);

            int physNodesPerElem = 1;
            for (size_t i = 0; i < physDimIndices.size() && i < phaseElemOrderRead.size(); ++i)
            {
              physNodesPerElem *= (phaseElemOrderRead[physDimIndices[i]] + 1);
            }

            int physNodeIdx = this->PhysicalNodeIndex;
            int maxNodeIdx = physNodesPerElem - 1;
            if (physNodeIdx < 0 || physNodeIdx > maxNodeIdx)
            {
              vtkWarningMacro("PhysicalNodeIndex = " << physNodeIdx
                << " is out of range [0, " << maxNodeIdx << "], clamping");
              physNodeIdx = std::min(physNodeIdx, maxNodeIdx);
              physNodeIdx = std::max(physNodeIdx, 0);
            }

            std::vector<vtkSmartPointer<vtkDoubleArray>> velFieldArrays;
            for (int c = 0; c < numComponents; ++c)
            {
              vtkNew<vtkDoubleArray> arr;
              std::string name = (c < static_cast<int>(componentNames.size()))
                                 ? componentNames[c]
                                 : varName + "_" + std::to_string(c);
              arr->SetName(name.c_str());
              arr->SetNumberOfTuples(velTotalNodes);
              velFieldArrays.push_back(arr);
            }

            int fullNodesPerElem = nodesPerElem;
            int velNodesActual = fullNodesPerElem / physNodesPerElem;

            // Generalized N-dimensional data loading loop
            for (int cellIdx = 0; cellIdx < velTotalCells; ++cellIdx)
            {
              // Decompose cell index into per-dimension indices
              std::vector<int> cellIndices = PSDecomposeIndex(cellIdx, velNumCellsVec);

              int velElemStartIdx = cellIdx * velNodesPerElement;

              for (int nodeIdx = 0; nodeIdx < velNodesPerElement; ++nodeIdx)
              {
                // Decompose node index into per-dimension indices
                std::vector<int> nodeIndices = PSDecomposeIndex(nodeIdx, velNodesPerDim);

                int velVtkLocalIdx = velWarpmToVTK[nodeIdx];
                int vtkPointIdx = velElemStartIdx + velVtkLocalIdx;

                // Compute flat data index (row-major order in HDF5 data)
                // Data layout: velocity cells are in row-major order (last dim varies fastest)
                int velCellFlatIdx = PSComputeFlatIndex(cellIndices, velNumCellsVec);
                int fullNodeIdx = physNodeIdx * velNodesActual + nodeIdx;
                int dataIdx = (velCellFlatIdx * fullNodesPerElem + fullNodeIdx) * numComponents;

                for (int c = 0; c < numComponents; ++c)
                {
                  if (dataIdx + c < static_cast<int>(slicedData.size()))
                  {
                    velFieldArrays[c]->SetValue(vtkPointIdx, slicedData[dataIdx + c]);
                  }
                }
              }
            }

            for (auto& arr : velFieldArrays)
            {
              phaseOutput->GetPointData()->AddArray(arr);
            }
          }
        }

        H5Gclose(phaseDomGroup);
      }
    }
  }

  if (domainsGroup >= 0) H5Gclose(domainsGroup);
  H5Gclose(varsGroup);
  H5Fclose(file);

  // Debug: final state of all outputs
  {
    std::ofstream debugFile("/tmp/_warpm_debug.txt", std::ios::app);
    debugFile << "=== Final output state ===\n";

    for (int port = 0; port < 3; ++port)
    {
      vtkInformation* info = outputVector->GetInformationObject(port);
      if (info)
      {
        vtkDataObject* obj = info->Get(vtkDataObject::DATA_OBJECT());
        debugFile << "Port " << port << ": ptr=" << obj;
        if (obj)
        {
          vtkDataSet* ds = vtkDataSet::SafeDownCast(obj);
          if (ds)
          {
            double bounds[6];
            ds->GetBounds(bounds);
            debugFile << ", pts=" << ds->GetNumberOfPoints()
                      << ", cells=" << ds->GetNumberOfCells()
                      << ", bounds=[" << bounds[0] << "," << bounds[1] << "]x["
                      << bounds[2] << "," << bounds[3] << "]";
          }
        }
        debugFile << "\n";
      }
    }
  }

  return 1;
}

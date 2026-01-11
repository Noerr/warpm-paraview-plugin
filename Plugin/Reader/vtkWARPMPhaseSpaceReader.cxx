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
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkDataArraySelection.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>

#include "vtk_hdf5.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <cmath>
#include <regex>

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

// Parse coordinate expression like "x=-6.28... + 2.513...*k;"
static bool PSParseCoordinateExpression(const std::string& expr, double& offset, double& scale)
{
  std::regex pattern(R"(x\s*=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*\+\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*\*\s*k)");
  std::smatch match;
  if (std::regex_search(expr, match, pattern)) {
    offset = std::stod(match[1].str());
    scale = std::stod(match[2].str());
    return true;
  }
  return false;
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

// Build mapping from WARPM row-major node index to VTK Lagrange node index
static std::vector<int> PSBuildWarpmToVTKLagrangeMapping(int orderX, int orderY)
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
  // Override base class: set 2 output ports
  this->SetNumberOfOutputPorts(2);
}

//----------------------------------------------------------------------------
vtkWARPMPhaseSpaceReader::~vtkWARPMPhaseSpaceReader()
{
}

//----------------------------------------------------------------------------
int vtkWARPMPhaseSpaceReader::FillOutputPortInformation(int port, vtkInformation* info)
{
  // Both ports output vtkUnstructuredGrid
  // Port 0: Physical space mesh
  // Port 1: Velocity space mesh
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
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
     << ", " << this->PhysicalSliceIndices[1] << "\n";
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

  for (const auto& name : allVarNames)
  {
    if (!this->PointDataArraySelection->ArrayIsEnabled(name.c_str()))
    {
      continue;
    }

    // Check which domain this variable uses
    bool isPhaseSpace = false;
    hid_t varGroup = H5Gopen(varsGroup, name.c_str(), H5P_DEFAULT);
    if (varGroup >= 0)
    {
      std::string domainName;
      if (PSReadStringAttribute(varGroup, "OnDomain", domainName) && domainsGroup >= 0)
      {
        hid_t domGroup = H5Gopen(domainsGroup, domainName.c_str(), H5P_DEFAULT);
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
      phaseSpaceVars.push_back(name);
    }
    else
    {
      physicalVars.push_back(name);
    }
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
    hid_t domain = H5Gopen(file, "/domains/physical_space_domain", H5P_DEFAULT);
    if (domain < 0)
    {
      vtkErrorMacro("Could not open physical_space_domain group");
      if (domainsGroup >= 0) H5Gclose(domainsGroup);
      H5Gclose(varsGroup);
      H5Fclose(file);
      return 0;
    }

    int ndims = 0;
    std::vector<int> numCells;
    std::vector<std::string> coordExprs;

    PSReadIntAttribute(domain, "ndims", ndims);
    PSReadIntArrayAttribute(domain, "numCells", numCells);
    PSReadStringArrayAttribute(domain, "VertexCoordinateExpressions", coordExprs);

    H5Gclose(domain);

    if (ndims < 2 || numCells.size() < 2 || coordExprs.size() < 2)
    {
      vtkErrorMacro("Invalid domain dimensions");
      if (domainsGroup >= 0) H5Gclose(domainsGroup);
      H5Gclose(varsGroup);
      H5Fclose(file);
      return 0;
    }

    double xOffset = 0, xScale = 1, yOffset = 0, yScale = 1;
    PSParseCoordinateExpression(coordExprs[0], xOffset, xScale);
    PSParseCoordinateExpression(coordExprs[1], yOffset, yScale);

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

    int dataNx = static_cast<int>(dataDims[0]);
    int dataNy = static_cast<int>(dataDims[1]);
    int nodesPerElement = static_cast<int>(dataDims[2]);

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
      auto gllX = PSGetGLLNodes(orderX);
      auto gllY = PSGetGLLNodes(orderY);
      this->CachedWarpmToVTK = PSBuildWarpmToVTKLagrangeMapping(orderX, orderY);

      this->CachedPoints = vtkSmartPointer<vtkPoints>::New();
      int totalNodes = dataNx * dataNy * nodesPerElement;
      this->CachedPoints->SetNumberOfPoints(totalNodes);

      for (int ey = 0; ey < dataNy; ++ey)
      {
        for (int ex = 0; ex < dataNx; ++ex)
        {
          double elemXMin = xOffset + xScale * ex;
          double elemXMax = xOffset + xScale * (ex + 1);
          double elemYMin = yOffset + yScale * ey;
          double elemYMax = yOffset + yScale * (ey + 1);

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

      this->CachedNx = dataNx;
      this->CachedNy = dataNy;
      this->CachedNodesPerElement = nodesPerElement;
      this->GeometryCached = true;
    }

    output->SetPoints(this->CachedPoints);
    output->SetCells(VTK_LAGRANGE_QUADRILATERAL, this->CachedCells);

    int totalNodes = dataNx * dataNy * nodesPerElement;

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

      int varNumComponents = static_cast<int>(varDims[3]);

      std::vector<double> data(dataNx * dataNy * nodesPerElement * varNumComponents);
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

    if (phaseOutput && domainsGroup >= 0)
    {
      hid_t phaseDomGroup = H5Gopen(domainsGroup, this->PhaseSpaceDomainName.c_str(), H5P_DEFAULT);
      if (phaseDomGroup >= 0)
      {
        std::vector<int> phaseNumCells;
        std::vector<std::string> phaseCoordExprs;
        PSReadIntArrayAttribute(phaseDomGroup, "numCells", phaseNumCells);
        PSReadStringArrayAttribute(phaseDomGroup, "VertexCoordinateExpressions", phaseCoordExprs);

        std::vector<std::string> coordNames;
        PSReadStringArrayAttribute(phaseDomGroup, "CoordinateNames", coordNames);

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

        if (velDimIndices.size() >= 2)
        {
          int velNx = phaseNumCells[velDimIndices[0]];
          int velNy = phaseNumCells[velDimIndices[1]];

          double velXOffset = 0, velXScale = 1, velYOffset = 0, velYScale = 1;
          if (velDimIndices[0] < static_cast<int>(phaseCoordExprs.size()))
            PSParseCoordinateExpression(phaseCoordExprs[velDimIndices[0]], velXOffset, velXScale);
          if (velDimIndices[1] < static_cast<int>(phaseCoordExprs.size()))
            PSParseCoordinateExpression(phaseCoordExprs[velDimIndices[1]], velYOffset, velYScale);

          hid_t firstPhaseVar = H5Gopen(varsGroup, phaseSpaceVars[0].c_str(), H5P_DEFAULT);
          std::vector<int> phaseElemOrder;
          PSReadIntArrayAttribute(firstPhaseVar, "ElementOrder", phaseElemOrder);

          int velOrderX = phaseElemOrder.size() > static_cast<size_t>(velDimIndices[0])
                          ? phaseElemOrder[velDimIndices[0]] : 2;
          int velOrderY = phaseElemOrder.size() > static_cast<size_t>(velDimIndices[1])
                          ? phaseElemOrder[velDimIndices[1]] : 2;
          int velNodesX = velOrderX + 1;
          int velNodesY = velOrderY + 1;
          int velNodesPerElement = velNodesX * velNodesY;

          H5Gclose(firstPhaseVar);

          auto gllVelX = PSGetGLLNodes(velOrderX);
          auto gllVelY = PSGetGLLNodes(velOrderY);
          auto velWarpmToVTK = PSBuildWarpmToVTKLagrangeMapping(velOrderX, velOrderY);

          vtkNew<vtkPoints> velPoints;
          int velTotalNodes = velNx * velNy * velNodesPerElement;
          velPoints->SetNumberOfPoints(velTotalNodes);

          for (int ey = 0; ey < velNy; ++ey)
          {
            for (int ex = 0; ex < velNx; ++ex)
            {
              double elemVxMin = velXOffset + velXScale * ex;
              double elemVxMax = velXOffset + velXScale * (ex + 1);
              double elemVyMin = velYOffset + velYScale * ey;
              double elemVyMax = velYOffset + velYScale * (ey + 1);

              int elemStartIdx = (ey * velNx + ex) * velNodesPerElement;

              for (int ix = 0; ix < velNodesX; ++ix)
              {
                for (int iy = 0; iy < velNodesY; ++iy)
                {
                  int warpmLocalIdx = ix * velNodesY + iy;
                  int vtkLocalIdx = velWarpmToVTK[warpmLocalIdx];
                  int vtkPointIdx = elemStartIdx + vtkLocalIdx;

                  double xi = gllVelX[ix];
                  double eta = gllVelY[iy];
                  double vx = elemVxMin + (xi + 1.0) * 0.5 * (elemVxMax - elemVxMin);
                  double vy = elemVyMin + (eta + 1.0) * 0.5 * (elemVyMax - elemVyMin);

                  velPoints->SetPoint(vtkPointIdx, vx, vy, 0.0);
                }
              }
            }
          }

          vtkNew<vtkCellArray> velCells;
          for (int ey = 0; ey < velNy; ++ey)
          {
            for (int ex = 0; ex < velNx; ++ex)
            {
              int elemStartIdx = (ey * velNx + ex) * velNodesPerElement;
              std::vector<vtkIdType> cellPts(velNodesPerElement);
              for (int i = 0; i < velNodesPerElement; ++i)
              {
                cellPts[i] = elemStartIdx + i;
              }
              velCells->InsertNextCell(velNodesPerElement, cellPts.data());
            }
          }

          phaseOutput->SetPoints(velPoints);
          phaseOutput->SetCells(VTK_LAGRANGE_QUADRILATERAL, velCells);

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

            for (size_t i = 0; i < physDimIndices.size() && i < 2; ++i)
            {
              int idx = physDimIndices[i];
              int sliceVal = (i == 0) ? this->PhysicalSliceIndices[0] : this->PhysicalSliceIndices[1];
              sliceVal = std::min(sliceVal, static_cast<int>(varDims[idx]) - 1);
              sliceVal = std::max(sliceVal, 0);
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

            int physNodeIdx = std::min(this->PhysicalNodeIndex, physNodesPerElem - 1);
            physNodeIdx = std::max(physNodeIdx, 0);

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

            for (int vey = 0; vey < velNy; ++vey)
            {
              for (int vex = 0; vex < velNx; ++vex)
              {
                int velElemStartIdx = (vey * velNx + vex) * velNodesPerElement;

                for (int vix = 0; vix < velNodesX; ++vix)
                {
                  for (int viy = 0; viy < velNodesY; ++viy)
                  {
                    int velWarpmLocalIdx = vix * velNodesY + viy;
                    int velVtkLocalIdx = velWarpmToVTK[velWarpmLocalIdx];
                    int vtkPointIdx = velElemStartIdx + velVtkLocalIdx;

                    int fullNodeIdx = physNodeIdx * velNodesActual + velWarpmLocalIdx;
                    int dataIdx = ((vex * velNy + vey) * fullNodesPerElem + fullNodeIdx) * numComponents;

                    for (int c = 0; c < numComponents; ++c)
                    {
                      if (dataIdx + c < static_cast<int>(slicedData.size()))
                      {
                        velFieldArrays[c]->SetValue(vtkPointIdx, slicedData[dataIdx + c]);
                      }
                    }
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

  return 1;
}

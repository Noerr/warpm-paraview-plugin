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

#include "vtk_hdf5.h"

#include <algorithm>

vtkStandardNewMacro(vtkWARPMPhaseSpaceReader);

//----------------------------------------------------------------------------
vtkWARPMPhaseSpaceReader::vtkWARPMPhaseSpaceReader()
  : PhysicalNodeIndex(0)
  , FileHasPhaseSpace(false)
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
  std::string h5path = ResolveWARPMFile(fname);
  if (h5path.empty())
  {
    return 0;
  }
  // Phase space reader handles files WITH phase space variables
  return HasPhaseSpaceVariables(h5path.c_str()) ? 1 : 0;
}

//----------------------------------------------------------------------------
void vtkWARPMPhaseSpaceReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "PhysicalSliceIndices: " << this->PhysicalSliceIndices[0]
     << ", " << this->PhysicalSliceIndices[1]
     << ", " << this->PhysicalSliceIndices[2] << "\n";
  os << indent << "PhysicalNodeIndex: " << this->PhysicalNodeIndex << "\n";
  os << indent << "HasPhaseSpaceVariables: " << this->FileHasPhaseSpace << "\n";
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
  this->FileHasPhaseSpace = false;
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
            if (ReadStringAttribute(varGroup, "OnDomain", domainName))
            {
              hid_t domGroup = H5Gopen(domainsGroup, domainName.c_str(), H5P_DEFAULT);
              if (domGroup >= 0)
              {
                int numPhys = 0, numVel = 0;
                std::vector<std::string> coordNames;
                if (DetectPhaseSpaceDomain(domGroup, numPhys, numVel, coordNames))
                {
                  // Found a phase space variable - record the domain info
                  if (!this->FileHasPhaseSpace)
                  {
                    this->FileHasPhaseSpace = true;
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
      if (ReadStringAttribute(varGroup, "OnDomain", varDomainName) && domainsGroup >= 0)
      {
        hid_t domGroup = H5Gopen(domainsGroup, varDomainName.c_str(), H5P_DEFAULT);
        if (domGroup >= 0)
        {
          int numPhys = 0, numVel = 0;
          std::vector<std::string> coordNames;
          isPhaseSpace = DetectPhaseSpaceDomain(domGroup, numPhys, numVel, coordNames);
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

    ReadIntAttribute(domain, "ndims", ndims);
    ReadIntArrayAttribute(domain, "numCells", numCells);
    ReadIntArrayAttribute(domain, "startIndices", startIndices);
    ReadStringArrayAttribute(domain, "VertexCoordinateExpressions", coordExprs);

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

    ReadIntArrayAttribute(firstVarGroup, "ElementOrder", elementOrder);
    ReadIntAttribute(firstVarGroup, "EntriesPerElement", entriesPerElement);

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

    for (int d = 0; d < numPhysDims; ++d)
    {
      physNumCellsVec[d] = numCells[d];
      physOrders[d] = (d < static_cast<int>(elementOrder.size())) ? elementOrder[d] : 2;
      physVertices[d] = ComputeVertexPositions(coordExprs[d], startIndices[d], numCells[d]);
    }

    int nodesPerElement = ComputeTotalNodes(physOrders);
    int totalCells = 1;
    for (int nc : physNumCellsVec) totalCells *= nc;
    int totalNodes = totalCells * nodesPerElement;

    // Check if we can use cached geometry
    bool canUseCache = this->GeometryCached &&
                       this->CachedNodesPerElement == nodesPerElement &&
                       this->CachedNdims == numPhysDims &&
                       this->CachedDataDims == physNumCellsVec;

    if (!canUseCache)
    {
      // Build geometry using shared helper
      this->CachedPoints = vtkSmartPointer<vtkPoints>::New();
      this->CachedCells = vtkSmartPointer<vtkCellArray>::New();
      BuildLagrangeMesh(physVertices, physOrders, physNumCellsVec,
                        this->CachedPoints, this->CachedCells, this->CachedWarpmToVTK);

      this->CachedNdims = numPhysDims;
      this->CachedDataDims = physNumCellsVec;
      this->CachedNodesPerElement = nodesPerElement;
      this->GeometryCached = true;
    }

    output->SetPoints(this->CachedPoints);
    output->SetCells(GetLagrangeCellType(numPhysDims), this->CachedCells);

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
      ReadStringArrayAttribute(varGroup, "ComponentNames", componentNames);

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
        std::vector<int> cellIndices = DecomposeIndex(cellIdx, physNumCellsVec);

        int elemStartIdx = cellIdx * nodesPerElement;

        for (int nodeIdx = 0; nodeIdx < nodesPerElement; ++nodeIdx)
        {
          int vtkLocalIdx = this->CachedWarpmToVTK[nodeIdx];
          int vtkPointIdx = elemStartIdx + vtkLocalIdx;

          // Compute flat data index (row-major order in HDF5 data)
          int dataCellIdx = ComputeFlatIndex(cellIndices, physNumCellsVec);
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

    if (phaseOutput && domainsGroup >= 0 && !phaseSpaceDomainName.empty())
    {
      hid_t phaseDomGroup = H5Gopen(domainsGroup, phaseSpaceDomainName.c_str(), H5P_DEFAULT);
      if (phaseDomGroup >= 0)
      {
        std::vector<int> phaseNumCells;
        std::vector<int> phaseStartIndices;
        std::vector<std::string> phaseCoordExprs;
        ReadIntArrayAttribute(phaseDomGroup, "numCells", phaseNumCells);
        ReadIntArrayAttribute(phaseDomGroup, "startIndices", phaseStartIndices);
        ReadStringArrayAttribute(phaseDomGroup, "VertexCoordinateExpressions", phaseCoordExprs);

        std::vector<std::string> coordNames;
        ReadStringArrayAttribute(phaseDomGroup, "CoordinateNames", coordNames);

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

          // Get element orders from first phase space variable
          hid_t firstPhaseVar = H5Gopen(varsGroup, phaseSpaceVars[0].c_str(), H5P_DEFAULT);
          std::vector<int> phaseElemOrder;
          ReadIntArrayAttribute(firstPhaseVar, "ElementOrder", phaseElemOrder);
          H5Gclose(firstPhaseVar);

          // Compute vertex positions and setup for each velocity dimension
          for (int d = 0; d < numVelDims; ++d)
          {
            int idx = velDimIndices[d];
            velNumCellsVec[d] = phaseNumCells[idx];
            velOrders[d] = (idx < static_cast<int>(phaseElemOrder.size()))
                           ? phaseElemOrder[idx] : 2;
            velVertices[d] = ComputeVertexPositions(
              phaseCoordExprs[idx], phaseStartIndices[idx], phaseNumCells[idx]);
          }

          // Compute total cells and nodes
          int velNodesPerElement = ComputeTotalNodes(velOrders);
          int velTotalCells = 1;
          for (int nc : velNumCellsVec) velTotalCells *= nc;
          int velTotalNodes = velTotalCells * velNodesPerElement;

          // Build nodes-per-dimension array (needed for data loading loop)
          std::vector<int> velNodesPerDim(numVelDims);
          for (int d = 0; d < numVelDims; ++d)
          {
            velNodesPerDim[d] = velOrders[d] + 1;
          }

          // Build velocity mesh using shared helper
          vtkNew<vtkPoints> velPoints;
          vtkNew<vtkCellArray> velCells;
          std::vector<int> velWarpmToVTK;
          BuildLagrangeMesh(velVertices, velOrders, velNumCellsVec, velPoints, velCells, velWarpmToVTK);

          // Set mesh with appropriate cell type
          phaseOutput->SetPoints(velPoints);
          phaseOutput->SetCells(GetLagrangeCellType(numVelDims), velCells);

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
            ReadIntArrayAttribute(probePhaseVar, "ElementOrder", probeElemOrder);
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
              physVertices[d] = ComputeVertexPositions(
                phaseCoordExprs[idx], phaseStartIndices[idx], phaseNumCells[idx]);
              physGLL[d] = GetGLLNodes(physOrders[d]);

              // Get slice index with bounds checking
              int maxSlice = phaseNumCells[idx] - 1;
              sliceIndices[d] = std::max(0, std::min(this->PhysicalSliceIndices[d], maxSlice));
            }

            // Compute total nodes per physical cell and decompose node index
            int physNodesPerCell = ComputeTotalNodes(physOrders);
            std::vector<int> physNodesPerDim(numPhysDims);
            for (int d = 0; d < numPhysDims; ++d)
            {
              physNodesPerDim[d] = physOrders[d] + 1;
            }

            int physNodeIdx = this->PhysicalNodeIndex;
            int maxNodeIdx = physNodesPerCell - 1;
            physNodeIdx = std::max(0, std::min(physNodeIdx, maxNodeIdx));

            // Decompose node index into per-dimension indices
            std::vector<int> nodeIndices = DecomposeIndex(physNodeIdx, physNodesPerDim);

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
            ReadStringArrayAttribute(varGroup, "ComponentNames", componentNames);

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
            ReadIntArrayAttribute(varGroupReopen, "ElementOrder", phaseElemOrderRead);
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
              std::vector<int> cellIndices = DecomposeIndex(cellIdx, velNumCellsVec);

              int velElemStartIdx = cellIdx * velNodesPerElement;

              for (int nodeIdx = 0; nodeIdx < velNodesPerElement; ++nodeIdx)
              {
                // Decompose node index into per-dimension indices
                std::vector<int> nodeIndices = DecomposeIndex(nodeIdx, velNodesPerDim);

                int velVtkLocalIdx = velWarpmToVTK[nodeIdx];
                int vtkPointIdx = velElemStartIdx + velVtkLocalIdx;

                // Compute flat data index (row-major order in HDF5 data)
                // Data layout: velocity cells are in row-major order (last dim varies fastest)
                int velCellFlatIdx = ComputeFlatIndex(cellIndices, velNumCellsVec);
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

  return 1;
}

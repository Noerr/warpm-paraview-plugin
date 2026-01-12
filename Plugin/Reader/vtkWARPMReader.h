/**
 * @file vtkWARPMReader.h
 * @brief VTK reader for WARPM HDF5 simulation output
 *
 * This reader loads WARPM simulation data from HDF5 files. It supports:
 * - Direct .h5 file loading (single frame or file series via ParaView)
 * - .warpm metadata files (time series with file listing)
 *
 * High-order discontinuous Galerkin (DG) elements are mapped to VTK Lagrange
 * cell types for accurate visualization.
 *
 * For phase space (Vlasov-Maxwell) data, use vtkWARPMPhaseSpaceReader instead.
 */

#ifndef vtkWARPMReader_h
#define vtkWARPMReader_h

#include "WARPMReaderCoreModule.h" // For export macro
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkSmartPointer.h>

// Forward declaration of hid_t for protected helper method signatures.
// The actual HDF5 header is included in the .cxx files.
#ifndef H5public_H
#include <cstdint>
using hid_t = int64_t;
#endif
#include <string>
#include <vector>

class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;

class WARPMREADERCORE_EXPORT vtkWARPMReader : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkWARPMReader* New();
  vtkTypeMacro(vtkWARPMReader, vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  ///@{
  /**
   * Specify the WARPM data file (.h5 or .warpm).
   */
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  ///@}

  /**
   * Test whether the file can be read by this reader.
   * Validates WARPM-specific HDF5 structure (/domains, /variables, /timeData).
   * For .warpm files, parses to get the first .h5 file and validates that.
   */
  static int CanReadFile(const char* fname);

  ///@{
  /**
   * Get the data array selection object for point arrays (variables).
   */
  vtkDataArraySelection* GetPointDataArraySelection();
  ///@}

  /**
   * Methods for ParaView array selection UI.
   */
  int GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int GetPointArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, int status);

protected:
  vtkWARPMReader();
  ~vtkWARPMReader() override;

  /**
   * Called to report time steps and other metadata.
   */
  int RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector) override;

  /**
   * Called to generate the output data.
   */
  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector) override;

  int FillOutputPortInformation(int port, vtkInformation* info) override;

  char* FileName;

  /**
   * Parse the .warpm file to get list of HDF5 frame files.
   */
  bool ParseWarpmFile();

  /**
   * List of HDF5 frame file paths (absolute).
   */
  std::vector<std::string> FrameFiles;

  /**
   * Physical time values for each frame (read from /timeData/time).
   */
  std::vector<double> TimeValues;

  /**
   * Directory containing the .warpm file (used for relative paths).
   */
  std::string BaseDirectory;

  /**
   * Read time values from all HDF5 frame files.
   */
  bool ReadTimeValues();

  /**
   * Find the frame index closest to the requested time.
   */
  int FindClosestFrame(double requestedTime) const;

  // Geometry cache - mesh doesn't change between frames
  vtkSmartPointer<vtkPoints> CachedPoints;
  vtkSmartPointer<vtkCellArray> CachedCells;
  bool GeometryCached = false;

  // Cached mesh parameters for validation
  int CachedNdims = 0;
  std::vector<int> CachedDataDims;      // Number of cells in each dimension
  int CachedNodesPerElement = 0;
  std::vector<int> CachedWarpmToVTK;

  // Variable (point array) selection
  vtkSmartPointer<vtkDataArraySelection> PointDataArraySelection;

  ///@{
  /**
   * HDF5 attribute reading helpers.
   * These read attributes from HDF5 objects and return true on success.
   */
  static bool ReadIntAttribute(hid_t loc, const char* name, int& value);
  static bool ReadIntArrayAttribute(hid_t loc, const char* name, std::vector<int>& values);
  static bool ReadDoubleAttribute(hid_t loc, const char* name, double& value);
  static bool ReadStringAttribute(hid_t loc, const char* name, std::string& value);
  static bool ReadStringArrayAttribute(hid_t loc, const char* name, std::vector<std::string>& values);
  ///@}

  ///@{
  /**
   * Domain analysis helpers.
   */
  static bool IsPhysicalOnlyDomain(hid_t domainGroup);
  static bool IsVariableOnPhysicalDomain(hid_t file, hid_t varGroup);
  static bool DetectPhaseSpaceDomain(hid_t domainGroup, int& numPhysical, int& numVelocity,
                                     std::vector<std::string>& coordNames);
  ///@}

  ///@{
  /**
   * File validation helpers.
   */
  static bool IsWarpmHDF5File(const char* fname);
  static bool HasPhaseSpaceVariables(const char* fname);
  /**
   * Resolve a .h5 or .warpm file to its underlying HDF5 path.
   * For .h5 files, validates WARPM structure and returns the path.
   * For .warpm files, parses to find the first .h5 and validates it.
   * Returns empty string if file is invalid or not a WARPM file.
   */
  static std::string ResolveWARPMFile(const char* fname);
  ///@}

  ///@{
  /**
   * Coordinate and mesh computation helpers.
   */
  static std::vector<double> ComputeVertexPositions(const std::string& expression,
                                                    int startIndex, int numCells);
  static std::vector<double> GetGLLNodes(int order);
  ///@}

  ///@{
  /**
   * N-dimensional index helpers.
   */
  static std::vector<int> DecomposeIndex(int flatIdx, const std::vector<int>& sizes);
  static int ComputeFlatIndex(const std::vector<int>& indices, const std::vector<int>& sizes);
  static int ComputeTotalNodes(const std::vector<int>& orders);
  static int GetLagrangeCellType(int ndims);
  ///@}

  ///@{
  /**
   * WARPM-to-VTK node ordering mapping builders.
   * These convert from WARPM's row-major tensor ordering to VTK Lagrange ordering.
   */
  static std::vector<int> BuildWarpmToVTKMapping(const std::vector<int>& orders);
  static std::vector<int> BuildWarpmToVTKMapping1D(int order);
  static std::vector<int> BuildWarpmToVTKMapping2D(int orderX, int orderY);
  static std::vector<int> BuildWarpmToVTKMapping3D(int orderX, int orderY, int orderZ);
  ///@}

private:
  vtkWARPMReader(const vtkWARPMReader&) = delete;
  void operator=(const vtkWARPMReader&) = delete;
};

#endif

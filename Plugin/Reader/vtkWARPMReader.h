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
 */

#ifndef vtkWARPMReader_h
#define vtkWARPMReader_h

#include "WARPMReaderCoreModule.h" // For export macro
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkSmartPointer.h>

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

  char* FileName;

private:
  vtkWARPMReader(const vtkWARPMReader&) = delete;
  void operator=(const vtkWARPMReader&) = delete;

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
  int CachedNx = 0;
  int CachedNy = 0;
  int CachedNodesPerElement = 0;
  std::vector<int> CachedWarpmToVTK;

  // Variable (point array) selection
  vtkSmartPointer<vtkDataArraySelection> PointDataArraySelection;
};

#endif

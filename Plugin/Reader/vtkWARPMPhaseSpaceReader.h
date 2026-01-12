/**
 * @file vtkWARPMPhaseSpaceReader.h
 * @brief VTK reader for WARPM HDF5 phase space (Vlasov-Maxwell) data
 *
 * This reader extends vtkWARPMReader to support phase space data where
 * variables exist on domains with both physical (x, y, z) and velocity
 * (vx, vy, vz or wx, wy, wz) dimensions.
 *
 * Output ports:
 * - Port 0: Physical space mesh with physical-space variables
 * - Port 1: Velocity space mesh with phase-space variables (sliced at specified physical location)
 * - Port 2: Probe location (single point showing where the slice is taken in physical space)
 */

#ifndef vtkWARPMPhaseSpaceReader_h
#define vtkWARPMPhaseSpaceReader_h

#include "WARPMReaderCoreModule.h" // For export macro
#include "vtkWARPMReader.h"

#include <string>
#include <vector>

class WARPMREADERCORE_EXPORT vtkWARPMPhaseSpaceReader : public vtkWARPMReader
{
public:
  static vtkWARPMPhaseSpaceReader* New();
  vtkTypeMacro(vtkWARPMPhaseSpaceReader, vtkWARPMReader);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Test whether the file can be read by this reader.
   * Returns true ONLY if the file contains phase space variables
   * (variables on domains with velocity dimensions).
   */
  static int CanReadFile(const char* fname);

  ///@{
  /**
   * Physical slice indices for phase space variables.
   * Specifies which cell in the physical dimensions to extract.
   */
  vtkGetVector2Macro(PhysicalSliceIndices, int);
  vtkSetVector2Macro(PhysicalSliceIndices, int);
  ///@}

  ///@{
  /**
   * Node index within the physical cell for phase space variables.
   * Selects which GLL node in the physical slice to use.
   */
  vtkGetMacro(PhysicalNodeIndex, int);
  vtkSetMacro(PhysicalNodeIndex, int);
  ///@}

protected:
  vtkWARPMPhaseSpaceReader();
  ~vtkWARPMPhaseSpaceReader() override;

  int RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector) override;

  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector) override;

  int FillOutputPortInformation(int port, vtkInformation* info) override;

  // Phase space slice properties
  int PhysicalSliceIndices[2];
  int PhysicalNodeIndex;

  // Phase space detection (set in RequestInformation)
  bool HasPhaseSpaceVariables;
  int NumPhysicalDims;
  int NumVelocityDims;
  std::string PhaseSpaceDomainName;

private:
  vtkWARPMPhaseSpaceReader(const vtkWARPMPhaseSpaceReader&) = delete;
  void operator=(const vtkWARPMPhaseSpaceReader&) = delete;
};

#endif

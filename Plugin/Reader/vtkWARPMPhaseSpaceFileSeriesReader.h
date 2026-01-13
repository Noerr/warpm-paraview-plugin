/**
 * @file vtkWARPMPhaseSpaceFileSeriesReader.h
 * @brief File series reader subclass that auto-detects output port count
 *
 * vtkFileSeriesReader hardcodes SetNumberOfOutputPorts(1) in its constructor,
 * but the rest of the class is generic with respect to output count.
 *
 * This subclass overrides SetReader() to automatically match the output port
 * count of the internal reader, enabling file series support for readers
 * with any number of output ports.
 *
 * This approach is a prototype for a proposed upstream change to vtkFileSeriesReader.
 * See: https://discourse.paraview.org/t/vtkfileseriesreader-multi-output-port-support/17394
 */

#ifndef vtkWARPMPhaseSpaceFileSeriesReader_h
#define vtkWARPMPhaseSpaceFileSeriesReader_h

#include "WARPMReaderCoreModule.h" // For export macro
#include <vtkFileSeriesReader.h>

class vtkAlgorithm;

class WARPMREADERCORE_EXPORT vtkWARPMPhaseSpaceFileSeriesReader : public vtkFileSeriesReader
{
public:
  static vtkWARPMPhaseSpaceFileSeriesReader* New();
  vtkTypeMacro(vtkWARPMPhaseSpaceFileSeriesReader, vtkFileSeriesReader);

  /**
   * Redefine SetReader to automatically match output port count.
   * This queries the reader's GetNumberOfOutputPorts() and configures
   * this file series reader to have the same number.
   * Note: Can use 'override' once MR !7637 is merged into ParaView.
   */
  void SetReader(vtkAlgorithm* reader);

protected:
  vtkWARPMPhaseSpaceFileSeriesReader();
  ~vtkWARPMPhaseSpaceFileSeriesReader() override = default;

private:
  vtkWARPMPhaseSpaceFileSeriesReader(const vtkWARPMPhaseSpaceFileSeriesReader&) = delete;
  void operator=(const vtkWARPMPhaseSpaceFileSeriesReader&) = delete;
};

#endif

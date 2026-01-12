/**
 * @file vtkWARPMPhaseSpaceFileSeriesReader.h
 * @brief File series reader subclass for multi-output-port phase space reader
 *
 * vtkFileSeriesReader hardcodes SetNumberOfOutputPorts(1) in its constructor,
 * but the rest of the class is generic with respect to output count.
 * This subclass simply overrides the constructor to set 3 output ports,
 * enabling file series support for the phase space reader.
 *
 * See: https://discourse.paraview.org/t/file-series-support-for-custom-time-unaware-plugin/15566
 */

#ifndef vtkWARPMPhaseSpaceFileSeriesReader_h
#define vtkWARPMPhaseSpaceFileSeriesReader_h

#include "WARPMReaderCoreModule.h" // For export macro
#include <vtkFileSeriesReader.h>

class WARPMREADERCORE_EXPORT vtkWARPMPhaseSpaceFileSeriesReader : public vtkFileSeriesReader
{
public:
  static vtkWARPMPhaseSpaceFileSeriesReader* New();
  vtkTypeMacro(vtkWARPMPhaseSpaceFileSeriesReader, vtkFileSeriesReader);

protected:
  vtkWARPMPhaseSpaceFileSeriesReader();
  ~vtkWARPMPhaseSpaceFileSeriesReader() override = default;

private:
  vtkWARPMPhaseSpaceFileSeriesReader(const vtkWARPMPhaseSpaceFileSeriesReader&) = delete;
  void operator=(const vtkWARPMPhaseSpaceFileSeriesReader&) = delete;
};

#endif

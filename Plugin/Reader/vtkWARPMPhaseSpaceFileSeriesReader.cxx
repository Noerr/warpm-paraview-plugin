/**
 * @file vtkWARPMPhaseSpaceFileSeriesReader.cxx
 * @brief Implementation of file series reader for multi-output-port phase space reader
 */

#include "vtkWARPMPhaseSpaceFileSeriesReader.h"
#include <vtkObjectFactory.h>

vtkStandardNewMacro(vtkWARPMPhaseSpaceFileSeriesReader);

//----------------------------------------------------------------------------
vtkWARPMPhaseSpaceFileSeriesReader::vtkWARPMPhaseSpaceFileSeriesReader()
{
  // Override the base class's SetNumberOfOutputPorts(1) to support 3 outputs:
  // Port 0: Physical space mesh
  // Port 1: Velocity space mesh
  // Port 2: Probe location marker
  this->SetNumberOfOutputPorts(3);
}

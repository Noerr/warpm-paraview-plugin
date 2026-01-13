/**
 * @file vtkWARPMPhaseSpaceFileSeriesReader.cxx
 * @brief Implementation of file series reader that auto-detects output port count
 *
 * This is a prototype for a proposed upstream change to vtkFileSeriesReader.
 */

#include "vtkWARPMPhaseSpaceFileSeriesReader.h"
#include <vtkAlgorithm.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(vtkWARPMPhaseSpaceFileSeriesReader);

//----------------------------------------------------------------------------
vtkWARPMPhaseSpaceFileSeriesReader::vtkWARPMPhaseSpaceFileSeriesReader()
{
  // Don't hardcode output ports here - let SetReader() auto-detect from the
  // internal reader. The base class sets 1 port, which will be updated when
  // SetReader() is called.
}

//----------------------------------------------------------------------------
void vtkWARPMPhaseSpaceFileSeriesReader::SetReader(vtkAlgorithm* reader)
{
  // Replicate vtkSetObjectMacro logic since we can't call Superclass::SetReader
  // until ParaView is rebuilt with our MR (which adds SetReader as an override)
  if (this->Reader != reader)
  {
    vtkAlgorithm* old = this->Reader;
    this->Reader = reader;
    if (reader)
    {
      reader->Register(this);
    }
    if (old)
    {
      old->UnRegister(this);
    }
    this->Modified();
  }

  // Auto-detect output port count from the internal reader
  if (reader)
  {
    int numPorts = reader->GetNumberOfOutputPorts();
    if (numPorts != this->GetNumberOfOutputPorts())
    {
      this->SetNumberOfOutputPorts(numPorts);
    }
  }
}

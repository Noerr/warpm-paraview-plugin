set(classes
  vtkWARPMReader
  vtkWARPMPhaseSpaceReader)

vtk_module_add_module(WARPMReaderCore
  FORCE_STATIC
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS WARPMReader.xml)

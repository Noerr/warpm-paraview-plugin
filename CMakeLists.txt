cmake_minimum_required(VERSION 3.16)
project(WARPMReaderPlugin)

# Fix dynamic_cast errors when loading plugin in shipping ParaView
# See: https://gitlab.kitware.com/paraview/paraview/-/issues/14869
# The type_info for plugin interfaces must have public visibility
set(CMAKE_CXX_VISIBILITY_PRESET default)
set(CMAKE_VISIBILITY_INLINES_HIDDEN OFF)
# Force visibility flags directly in case VTK/ParaView overrides the above
add_compile_options(-fvisibility=default)

find_package(ParaView REQUIRED)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

set("_paraview_plugin_default_WARPMReader" ON)
paraview_plugin_scan(
  PLUGIN_FILES      "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/paraview.plugin"
  PROVIDES_PLUGINS  plugins
  REQUIRES_MODULES  required_modules)

foreach (module IN LISTS required_modules)
  if (NOT TARGET "${module}")
    message(FATAL_ERROR
      "Missing required module: ${module}")
  endif ()
endforeach ()

paraview_plugin_build(
  RUNTIME_DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY_SUBDIRECTORY "${PARAVIEW_PLUGIN_SUBDIR}"
  PLUGINS ${plugins})

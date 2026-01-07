# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**WARPM ParaView Reader Plugin**

A C++ ParaView reader plugin that loads native WARPM HDF5 simulation output with high-order Lagrange element support. Users access time series data via a `.warpm` metadata file.

### Goals
1. Eliminate the need for Python-based HDF5-to-VTK conversion
2. Support high-order DG elements using VTK Lagrange cell types
3. Provide native ParaView time series support

## Build System

### Prerequisites

**ParaView Development Environment**: The binary ParaView app doesn't include C++ development headers. We build ParaView from source:

```
~/Downloads/paraview_localbuild/paraview/       # ParaView source
~/Downloads/paraview_localbuild/paraview/build/ # Build directory
```

**Verified versions** (January 2025):
- ParaView 6.0
- VTK 9.5.2
- Python 3.14

**Homebrew Dependencies**:
```bash
brew install cmake ninja qt@6 hdf5 tbb
```

Note: `hdf5-mpi` (used by WARPM) is compatible with the build.

### Building the Plugin

```bash
cd /Users/noah/Documents/warpm/warpm_paraview_plugin
mkdir build
cmake -S . -B build -DParaView_DIR=$HOME/Downloads/paraview_localbuild/paraview/build -Wno-dev
make -C build -j4
```

Output: `build/lib/paraview-6.0/plugins/WARPMReader/WARPMReader.so`

### Loading in ParaView

**GUI:**
1. Tools > Manage Plugins
2. Load New > Browse to `build/lib/paraview-6.0/plugins/WARPMReader/WARPMReader.so`
3. Optionally check "Auto Load"

**pvpython:**
```python
import paraview.simple as pv
pv.LoadPlugin('/path/to/WARPMReader.so', remote=False)
reader = pv.WARPMReader(FileName='/path/to/file.warpm')
```

### macOS Security Note

First launch of ParaView from the developer build may require authorization via a macOS security dialog. Until authorized, pvpython may hang on `import vtk`.

## Plugin Architecture

### File Structure

ParaView 6.0 plugins require a specific structure with a VTK module subdirectory:

```
warpm_paraview_plugin/
├── CMakeLists.txt              # Top-level: paraview_plugin_scan/build
├── CLAUDE.md
├── README.md
└── Plugin/
    ├── paraview.plugin         # Plugin metadata (NAME, DESCRIPTION, REQUIRES_MODULES)
    ├── CMakeLists.txt          # Plugin definition: paraview_add_plugin()
    └── Reader/                 # VTK module subdirectory
        ├── vtk.module          # Module definition (NAME, DEPENDS, PRIVATE_DEPENDS)
        ├── CMakeLists.txt      # Module build: vtk_module_add_module()
        ├── vtkWARPMReader.h    # Reader class header
        ├── vtkWARPMReader.cxx  # Reader implementation
        └── WARPMReader.xml     # ServerManager XML (registers .warpm extension)
```

### Key Naming Constraint

**The plugin name and module name MUST be different.** Both create CMake targets, so identical names cause conflicts.

- Plugin name: `WARPMReader` (in `Plugin/CMakeLists.txt`)
- Module name: `WARPMReaderCore` (in `Plugin/Reader/vtk.module`)

The module name determines the export macro: `WARPMREADERCORE_EXPORT` and header `WARPMReaderCoreModule.h`.

### FileSeriesReader Wrapper Architecture

To support file series (opening multiple .h5 files as a time series), the XML defines **two proxies**:

1. **`WARPMReaderCore`** (`internal_sources` group) - Our actual `vtkWARPMReader` class
2. **`WARPMReader`** (`sources` group) - A `vtkFileSeriesReader` wrapper exposed to the GUI

The wrapper uses these key attributes:
```xml
<SourceProxy name="WARPMReader"
             class="vtkFileSeriesReader"
             si_class="vtkSIMetaReaderProxy"
             file_name_method="SetFileName">
  <SubProxy>
    <Proxy name="Reader" proxygroup="internal_sources" proxyname="WARPMReaderCore"/>
    <ExposedProperties>...</ExposedProperties>
  </SubProxy>
  <StringVectorProperty name="FileNames" command="AddFileName" .../>
</SourceProxy>
```

This enables:
- Opening file groups (e.g., `maxwell_2d_..h5`) as time series
- FileSeriesReader calls our reader once per file, aggregating time values
- Physical time values are read from each file's `/timeData/time` attribute

### CMake Structure

**Top-level CMakeLists.txt:**
```cmake
set("_paraview_plugin_default_WARPMReader" ON)  # Enable by default
paraview_plugin_scan(
  PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/paraview.plugin"
  PROVIDES_PLUGINS plugins
  REQUIRES_MODULES required_modules)
paraview_plugin_build(PLUGINS ${plugins} ...)
```

**Plugin/CMakeLists.txt:**
```cmake
paraview_add_plugin(WARPMReader
  REQUIRED_ON_CLIENT
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES WARPMReaderCore
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")
```

**Plugin/Reader/CMakeLists.txt:**
```cmake
vtk_module_add_module(WARPMReaderCore
  FORCE_STATIC
  CLASSES vtkWARPMReader)
paraview_add_server_manager_xmls(XMLS WARPMReader.xml)
```

### Key Classes

**vtkWARPMReader** (inherits `vtkUnstructuredGridAlgorithm`)
- `SetFileName()` / `GetFileName()` - Path to .h5 or .warpm file
- `CanReadFile()` - Validates WARPM HDF5 structure (/domains, /variables, /timeData)
- `RequestInformation()` - Reports available time steps to ParaView
- `RequestData()` - Reads HDF5 and populates VTK output
- `ParseWarpmFile()` - Parses .warpm metadata, supports glob patterns

## Opening WARPM Data

### Direct .h5 File Access

Open .h5 files directly in ParaView:
1. File > Open, navigate to your .h5 files
2. Select files using:
   - Single click for one file
   - Cmd/Ctrl-click to select multiple individual files
   - Shift-click to select a range
   - ParaView auto-detects numbered series (e.g., `maxwell_2d_..h5`)
3. Choose "WARPM HDF5 Data" as the reader
4. Optionally check "Set reader as default" for .h5 files

For multi-file selection, ParaView wraps the reader with FileSeriesReader to handle time series automatically.

The reader validates WARPM HDF5 structure via `CanReadFile()` (checks for `/domains`, `/variables`, `/timeData` groups).

### .warpm Metadata Files

Alternative to direct selection: a `.warpm` text file listing HDF5 frame files.

**Glob pattern** (matches all files):
```
# Comment lines start with #
maxwell_2d_*.h5
```

**Explicit listing** (one filename per line):
```
maxwell_2d_0000.h5
maxwell_2d_0001.h5
maxwell_2d_0002.h5
```

Paths are relative to the .warpm file location.

## WARPM HDF5 Structure

Each frame file contains:

```
/domains/physical_space_domain/
  - Kind: "rectilinearMesh_d" or "structuredQuadrilateral_d"
  - ndims: 2
  - numCells: [Nx, Ny]  (elements per dimension)
  - VertexCoordinateExpressions: parametric mesh expressions
  - PeriodicDirs: which dimensions are periodic

/variables/<variable_name>/
  - ElementOrder: [px, py]  (polynomial order per dimension)
  - ElementType: "WmGLTensorHypercubeRM" (row-major tensor)
  - EntriesPerElement: (px+1)*(py+1) total nodes per element
  - ComponentNames: ['Ex', 'Ey', 'Ez', ...]
  - /data: shape (Nx, Ny, nodes_per_element, num_components)

/timeData/
  - frame: integer frame number
  - time: physical time value

/run/
  - Software: "WarpM"
  - Version: "0.85"
```

### Data Layout Example (2D Maxwell, order-3)
- 5×5 grid of elements
- Each element: 4×4 = 16 Gauss-Lobatto nodes
- 6 components per node (Ex, Ey, Ez, Bx, By, Bz)
- Data shape: `(5, 5, 16, 6)`

## High-Order Element Support

### VTK Lagrange Cell Types
- 2D: `VTK_LAGRANGE_QUADRILATERAL` (type 70)
- 3D: `VTK_LAGRANGE_HEXAHEDRON` (type 72)

### Node Ordering

**WARPM uses row-major ordering** (y varies fastest, then x):
- Linear index = `x * nodesY + y`
- Reference: `../test_files/dg/nodal_dg_utils.py` (see `NodesQuadNDIM()`)

**VTK Lagrange cells require a specific ordering** (see `vtkHigherOrderQuadrilateral::PointIndexFromIJK`):
1. Corner vertices first (counter-clockwise from bottom-left: 0,0 → orderX,0 → orderX,orderY → 0,orderY)
2. Edge nodes - **ALL edges traverse in INCREASING parameter direction** (not counter-clockwise!):
   - Edge 0 (bottom, y=0): x from 1 to orderX-1
   - Edge 1 (right, x=orderX): y from 1 to orderY-1
   - Edge 2 (top, y=orderY): x from 1 to orderX-1 (NOT reversed)
   - Edge 3 (left, x=0): y from 1 to orderY-1 (NOT reversed)
3. Face nodes (3D only)
4. Interior nodes (x varies fastest, then y)

The plugin includes `BuildWarpmToVTKLagrangeMapping()` to convert between orderings.

**Key insight:** VTK's edge ordering does NOT follow a counter-clockwise traversal around the quad. All four edges enumerate nodes in the increasing parameter direction. This was verified against VTK source code.

### Reference for DG Node Utilities
- `../test_files/dg/nodal_dg_utils.py` - Current row-major utilities (GLL positions, interpolation)
- `../test_files/dg/nodal_dg_utils_ColumnMajor.py` - Legacy column-major version (historical)
- `../test_files/domains/structured_mesh_utils.py` - Mesh parsing utilities

## Implementation Phases

### Phase 1: Minimal Plugin Skeleton [COMPLETE ✓]
- Plugin compiles and loads in ParaView 6.0
- `.warpm` file type registered via ServerManager XML
- Placeholder mesh output (single quad with TestField)
- Verified working in both GUI and pvpython

### Phase 2: .warpm File Parsing [COMPLETE ✓]
- Parse metadata file for frame list
- Support glob patterns (e.g., `maxwell_2d_*.h5`) and explicit listings
- Report time steps to ParaView pipeline
- Time slider works in GUI

### Phase 3: HDF5 Reading [COMPLETE ✓]
- Added `VTK::hdf5` dependency to build system
- Read domain geometry from `/domains/physical_space_domain/`
- Read variable data arrays from `/variables/<name>/data`
- Compute physical coordinates from GLL reference nodes
- All field components (Ex, Ey, Ez, Bx, By, Bz) loaded as point data

### Phase 4: High-Order Lagrange Elements [COMPLETE ✓]
- Implemented `BuildWarpmToVTKLagrangeMapping()` for node reordering
- Using `VTK_LAGRANGE_QUADRILATERAL` (type 70) cells
- Order-3 elements render with 16 nodes per cell (4×4 GLL grid)
- GLL node positions for orders 1-8 included
- **Key fix:** VTK edge nodes all traverse in increasing parameter direction (verified against `vtkHigherOrderQuadrilateral::PointIndexFromIJK` source)

### Phase 5: Multi-Variable Support [COMPLETE ✓]
- Enumerate all variables from `/variables/` in HDF5
- Variable selection UI (checkboxes) in Properties panel - all enabled by default
- Load all enabled variables, each with its own component arrays (e.g., Ex, Ey, Ez...)
- Uses `vtkDataArraySelection` for standard ParaView array selection pattern

### Phase 6: Time Series Optimization [COMPLETE ✓]
- Read actual time values from `/timeData/time` (physical simulation time)
- Cache geometry across frames (points, cells, node mapping)
- Only reload field data when switching frames
- `FindClosestFrame()` maps requested time to correct frame index

### Validation & Comparison [COMPLETE ✓]
- Plugin output visually indistinguishable from Python conversion script (`convert_dg_to_structured_vtk.py`)
- Confirms correct node ordering and data mapping

## Test Data

```
/Users/noah/Downloads/test_macos_warpm/
├── maxwell_2d.warpm          # Metadata file
├── maxwell_2d_0000.h5        # Frame files
├── maxwell_2d_0001.h5
├── ...
└── maxwell_2d_0010.h5
```

## Reference Files

| Purpose | Path |
|---------|------|
| Python conversion script | `../test_files/dg/convert_dg_to_structured_vtk.py` |
| DG node utilities (current) | `../test_files/dg/nodal_dg_utils.py` |
| DG node utilities (legacy) | `../test_files/dg/nodal_dg_utils_ColumnMajor.py` |
| Mesh utilities | `../test_files/domains/structured_mesh_utils.py` |
| ParaView build | `~/Downloads/paraview_localbuild/paraview/build` |

## Legacy Workflow (for reference)

The old Python-based conversion:
```bash
source ~/Documents/warpm/warpm_venv/bin/activate
python ../test_files/dg/convert_dg_to_structured_vtk.py maxwell_2d_0000.h5
```

This produces `.vts` files with linear elements (high-order broken into sub-elements).

## ParaView Installation

**Developer Build** (required for plugin development):
```
~/Downloads/paraview_localbuild/paraview/build/bin/
├── paraview.app/     # GUI application
├── pvpython          # Python interpreter with VTK/ParaView
├── pvbatch           # Batch processing
└── pvserver          # Server for client-server mode
```

**Binary Distribution** (for reference, lacks dev headers):
```
/Applications/ParaView-6.0.1.app/Contents/bin/
```

## Known Limitations

### Plugin Only Works with Developer-Built ParaView

**The plugin does NOT work with the pre-built ParaView.app from Kitware** (tested with both signed and unsigned macOS builds). The plugin appears to load successfully, but fails at runtime with:

```
ERROR: Algorithm vtkWARPMReader did not create output for port 0 when asked
by REQUEST_DATA_OBJECT and does not specify any DATA_TYPE_NAME.
```

**Solution:** Use the ParaView built from source (`~/Downloads/paraview_localbuild/paraview/build/bin/paraview.app`).

### Root Cause: C++ Symbol Visibility

The actual error (visible in macOS Console.app, not the ParaView console) is:

```
dynamic_cast error 2: One or more of the following type_info's has hidden
visibility or is defined in more than one translation unit. They should all
have public visibility. 11vtkPVPlugin, 17WARPMReaderPlugin,
38vtkPVDynamicInitializerPluginInterface.
```

This is [GitLab issue #14869](https://gitlab.kitware.com/paraview/paraview/-/issues/14869). The shipping ParaView binaries are built with `-fvisibility=hidden`, which causes `dynamic_cast` across shared library boundaries to fail. When ParaView tries to cast our plugin to `vtkPVServerManagerPluginInterface`, the cast fails because the type_info symbols have mismatched visibility.

**Why we can't fix this on our side:**
- We tried building with `-fvisibility=default` and `CMAKE_CXX_VISIBILITY_PRESET=default`
- The `dynamic_cast` requires **both sides** to have matching visibility
- The shipping ParaView's interfaces have hidden visibility baked in
- Only Kitware can fix this by building releases with public visibility for plugin interfaces

**Diagnostic approach:**
- macOS Console.app shows the `dynamic_cast` errors (filter by process name "paraview")
- Local ParaView build: no Console errors, debug output appears in ParaView
- Shipping ParaView: `dynamic_cast error 2` in Console, no debug output (code never executes)

### Linux Alternative

For Linux distribution, Kitware provides the `kitware/paraview_org-plugin-devel` Docker image for building binary-compatible plugins. No equivalent exists for macOS.

## Debugging Tips

1. **Plugin won't load**: Check that it was built against the same ParaView version
2. **File type not recognized**: Verify `WARPMReader.xml` has correct `ReaderFactory` hints
3. **HDF5 issues**: Use `h5dump` or Python h5py to inspect file structure
4. **Node ordering wrong**: Compare against `nodal_dg_utils.py` (row-major, y varies fastest)
5. **CMake target name conflict**: Plugin and module names must differ
6. **pvpython hangs on import**: Authorize ParaView via GUI first (macOS security)
7. **Plugin fails in shipping ParaView**: Must use developer-built ParaView (see Known Limitations)

## Development Reference

### ParaView Plugin Examples

Located in ParaView source at `Examples/Plugins/`:
- `SMMyProxy/` - Plugin with custom VTK module (good reference)
- `MyTiffWriter/` - Simple writer plugin
- `AnalyzeNIfTIReaderWriter/` - Full reader/writer with module structure

### Testing with pvpython

```bash
~/Downloads/paraview_localbuild/paraview/build/bin/pvpython -c "
import paraview.simple as pv
pv.LoadPlugin('/path/to/WARPMReader.so', remote=False)
reader = pv.WARPMReader(FileName='/path/to/file.warpm')
reader.UpdatePipeline()
output = pv.servermanager.Fetch(reader)
print('Points:', output.GetNumberOfPoints())
print('Cells:', output.GetNumberOfCells())
"
```

## Future Work

### Error Checking and Robustness
- Add comprehensive error handling for malformed HDF5 files
- Validate HDF5 structure before reading (check for required groups/attributes)
- Graceful degradation when optional metadata is missing
- Better error messages for user troubleshooting

### HDF5 I/O Optimization
- Currently opens/closes every H5 file to read time values and other metadata
- Consider caching file handles or reading metadata in bulk
- Investigate HDF5 file handle pooling for repeated access patterns

### Parallel ParaView / Parallel HDF5
- Support for distributed datasets larger than single pvserver memory
- Integration with ParaView's parallel pipeline (multiple pvserver ranks)
- Leverage parallel HDF5 for concurrent reads
- Domain decomposition awareness (map WARPM patches to pvserver ranks)

### Vector Variable Support
- Currently all components are scalar arrays (Ex, Ey, Ez, Bx, By, Bz)
- ParaView can display vectors with magnitude/direction glyphs
- Need to group components: {Ex, Ey, Ez} → vector E, {Bx, By, Bz} → vector B
- May require encoding vector groupings in WARPM .h5 output
- Investigate VTK's `vtkAssignAttribute` or native 3-component arrays

### Linux Distribution
- Use Kitware's `kitware/paraview_org-plugin-devel` Docker image
- Build binary-compatible plugins for official Linux ParaView releases
- CI/CD pipeline for automated plugin builds

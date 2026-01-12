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
        ├── vtkWARPMReader.h    # Base reader class header
        ├── vtkWARPMReader.cxx  # Base reader implementation
        ├── vtkWARPMPhaseSpaceReader.h    # Phase space reader header
        ├── vtkWARPMPhaseSpaceReader.cxx  # Phase space reader implementation
        └── WARPMReader.xml     # ServerManager XML (registers both readers)
```

### Key Naming Constraint

**The plugin name and module name MUST be different.** Both create CMake targets, so identical names cause conflicts.

- Plugin name: `WARPMReader` (in `Plugin/CMakeLists.txt`)
- Module name: `WARPMReaderCore` (in `Plugin/Reader/vtk.module`)

The module name determines the export macro: `WARPMREADERCORE_EXPORT` and header `WARPMReaderCoreModule.h`.

### Two Reader Classes

The plugin provides two reader classes for different data types:

1. **`vtkWARPMReader`** - Base reader for physical-space data (single output port)
   - Handles standard electromagnetic field data (E, B fields)
   - Single output: physical mesh with field variables

2. **`vtkWARPMPhaseSpaceReader`** - Derived reader for Vlasov-Maxwell phase space data (three output ports)
   - Port 0 ("Physical Space"): Physical mesh with physical-space variables
   - Port 1 ("Velocity Space"): Velocity mesh at user-specified physical location
   - Port 2 ("Slice Position"): Single-point marker showing where velocity slice is taken
   - Slice controls select which physical cell/node to extract velocity data from

**Smart file detection:** The phase space reader's `CanReadFile()` returns true ONLY if the file contains phase space variables (domains with velocity coordinates). This prevents the reader selection dialog from appearing for typical physical-only files.

**Domain type detection:** Domains are identified by their `CoordinateNames` attribute, not by name:
- Physical dimensions: coordinates like "x", "y", "z"
- Velocity dimensions: coordinates starting with "v" or "w"

This makes the reader robust to arbitrary domain naming conventions.

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
  CLASSES vtkWARPMReader vtkWARPMPhaseSpaceReader)
paraview_add_server_manager_xmls(XMLS WARPMReader.xml)
```

### Key Classes

**vtkWARPMReader** (inherits `vtkUnstructuredGridAlgorithm`)
- `SetFileName()` / `GetFileName()` - Path to .h5 or .warpm file
- `CanReadFile()` - Validates WARPM HDF5 structure (/domains, /variables, /timeData)
- `RequestInformation()` - Reports available time steps to ParaView
- `RequestData()` - Reads HDF5 and populates VTK output
- `ParseWarpmFile()` - Parses .warpm metadata, supports glob patterns
- `IsPhysicalOnlyDomain()` - Checks if domain has no velocity coordinates
- `IsVariableOnPhysicalDomain()` - Checks if variable's domain is physical-only

**vtkWARPMPhaseSpaceReader** (inherits `vtkWARPMReader`)
- Three output ports: Physical Space (port 0), Velocity Space (port 1), Slice Position (port 2)
- `SetPhysicalSliceIndices()` - Cell indices for slicing phase space (currently 2 indices for 2D physical)
- `SetPhysicalNodeIndex()` - Node index within physical cell for slice
- `CanReadFile()` - Returns true ONLY if phase space variables exist
- Builds velocity-space mesh from phase space domain geometry
- Uses HDF5 hyperslab selection to efficiently read sliced data
- Port 2 outputs a single vtkPolyData point showing the slice location in physical space

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
/domains/<domain_name>/
  - Kind: "rectilinearMesh_d" or "structuredQuadrilateral_d"
  - ndims: 2 or 4 (physical + velocity dimensions)
  - numCells: [Nx, Ny] or [Nx, Ny, Nvx, Nvy]
  - CoordinateNames: ["x", "y"] or ["x", "y", "wx", "wy"]
  - VertexCoordinateExpressions: parametric mesh expressions
  - PeriodicDirs: which dimensions are periodic

/variables/<variable_name>/
  - OnDomain: "<domain_name>"
  - ElementOrder: [px, py, ...]  (polynomial order per dimension)
  - ElementType: "WmGLTensorHypercubeRM" (row-major tensor)
  - EntriesPerElement: product of (order+1) for all dimensions
  - ComponentNames: ['Ex', 'Ey', 'Ez', ...] or ['proton_pdf']
  - /data: shape (cells..., nodes_per_element, num_components)

/timeData/
  - frame: integer frame number
  - time: physical time value

/run/
  - Software: "WarpM"
  - Version: "0.85"
```

### Domain Type Detection

Domain types are identified by the `CoordinateNames` attribute (not by domain name):
- **Physical domain**: All coordinates are spatial (x, y, z)
- **Phase space domain**: Mix of physical (x, y) and velocity (vx, vy or wx, wy) coordinates

Velocity coordinates are identified by names starting with 'v' or 'w'.

**Multiple domain handling:** If a file contains variables on multiple physical domains (or multiple phase space domains), only variables from the first domain encountered are loaded. Variables on other domains are skipped with a warning.

### Data Layout Example (2D Maxwell, order-3)
- 5×5 grid of elements
- Each element: 4×4 = 16 Gauss-Lobatto nodes
- 6 components per node (Ex, Ey, Ez, Bx, By, Bz)
- Data shape: `(5, 5, 16, 6)`

### Phase Space Data Layout Example (2D-2V, order-2)
- Physical: 3×3 cells, Velocity: 6×6 cells
- Each element: 3×3×3×3 = 81 nodes (order-2 in 4D)
- Data shape: `(3, 3, 6, 6, 81, num_components)`

### Coordinate Expressions (VertexCoordinateExpressions)

WARPM stores mesh geometry as parametric expressions rather than explicit vertex coordinates. The `VertexCoordinateExpressions` attribute is an array of strings, one per dimension, defining how vertex positions are computed from cell indices.

**Format:** Each expression is a function `x(k)` where `k` is the vertex index. The expressions use WARPM's HIP-based expression language, which supports:
- Standard C arithmetic operators: `+`, `-`, `*`, `/`
- Math functions: `sin`, `cos`, `exp`, `log`, `sqrt`, `pow`, etc.
- The variable `k` representing the vertex index

**Example expressions:**
```
# Physical space (polynomial mesh stretching):
x=-5.42127659574468e-04 + 1.11523404255319e-03*k + -9.29361702127659e-05*k*k + 6.19574468085106e-05*k*k*k;

# Velocity space (symmetric about zero):
x= k*(1.48148148148148e-02 *k*k + 2.00000000000000e-01 );
```

**startIndices attribute:** The `k` values don't always start at 0. The `startIndices` attribute specifies the starting index for each dimension:
- Physical dimensions typically start at 0
- Velocity dimensions often use symmetric ranges (e.g., `startIndices = [-3, -3]` for a 6-cell velocity grid spanning k = -3 to +3)

**Evaluation:** The plugin uses VTK's embedded Python interpreter with numpy for vectorized expression evaluation:
```python
import numpy as np
from math import *
k = np.arange(startIndex, startIndex + numCells + 1, dtype=np.float64)
result = <expression_rhs>  # e.g., k*(0.0148*k*k + 0.2)
```

This approach:
- Handles arbitrary mathematical expressions (not just polynomials)
- Supports all standard math functions via `from math import *`
- Uses numpy vectorization for efficiency (single Python call per dimension)
- Correctly handles non-zero startIndices for symmetric velocity grids

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

### Phase 7: Phase Space Reader [COMPLETE ✓]
- Derived `vtkWARPMPhaseSpaceReader` class with three output ports
- Port 0: Physical space mesh with physical-only variables
- Port 1: Velocity space mesh at user-specified physical location
- Port 2: Slice position probe marker (single point in physical space)
- Smart `CanReadFile()` returns true only for files with velocity dimensions
- Coordinate-based domain detection (not hardcoded domain names)
- Slice controls: Physical Slice Indices and Physical Node Index
- HDF5 hyperslab selection for efficient sliced data reading
- Warning when variables span multiple domains (only first domain loaded)
- Python-based coordinate expression evaluation (handles polynomial and math expressions)

**Current limitation:** The phase space reader currently assumes 2D physical space + 2D velocity space (4D total). The code has hardcoded assumptions about 2 physical dimensions and 2 velocity dimensions in several places:
- `PhysicalSliceIndices` is a 2-element array
- Mesh building loops assume 2D (nested x/y loops)
- VTK cell types are 2D (VTK_LAGRANGE_QUADRILATERAL)

### Phase 8: Arbitrary Dimension Support [PLANNED]
Generalize from 2D+2D to support all valid combinations:
- Physical dimensions: 1D, 2D, or 3D
- Velocity dimensions: 1D, 2D, or 3D
- Phase space total: 2D to 6D

This requires:
- Dynamic `PhysicalSliceIndices` size based on detected physical dimensions
- Generalized mesh building for 1D/2D/3D (lines, quads, hexes)
- Appropriate VTK cell types: `VTK_LAGRANGE_CURVE`, `VTK_LAGRANGE_QUADRILATERAL`, `VTK_LAGRANGE_HEXAHEDRON`
- Generalized node ordering mappings for each dimension count
- Dynamic hyperslab selection based on actual dimension counts

## Test Data

```
/Users/noah/Downloads/test_macos_warpm/
├── maxwell_2d.warpm          # Metadata file for physical-space data
├── maxwell_2d_0000.h5        # Physical-space frame files
├── maxwell_2d_0001.h5
├── ...
├── maxwell_2d_0010.h5
└── phase_space_vars_0000.h5  # Phase space test file (2D-2V Vlasov-Maxwell)
```

The phase space test file contains:
- `EM_field_n`: Physical-space electromagnetic field data
- `pdf_n`: Phase space distribution function (proton PDF)

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

### Per-Output-Port Property Visibility Not Supported

The Phase Space Reader has properties (Physical Slice Indices, Physical Node Index) that only affect the Velocity Space output (Port 1), not the Physical Space output (Port 0). Ideally, these controls would only appear when the user selects the Velocity Space output in the Pipeline Browser.

**Current behavior:** The slice properties are always visible in the Properties panel, grouped under "Velocity Space Location" with documentation noting they only affect Port 1.

**Desired behavior:** Properties would be conditionally shown/hidden based on which output port is selected.

**Why this isn't possible:** ParaView's Properties panel is source-centric, not port-centric. The ServerManager XML schema doesn't support nesting properties inside `<OutputPort>` elements or using PropertyWidgetDecorator to check which output port is selected. All properties belong to the source proxy as a whole. This is a fundamental limitation of ParaView's architecture.

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

### Arbitrary Phase Space Dimensions
The current implementation assumes 2D physical + 2D velocity (4D phase space). Generalization needed:
- Support 1D, 2D, or 3D physical space
- Support 1D, 2D, or 3D velocity space
- Total phase space dimensions: 2 to 6
- See "Phase 8" in Implementation Phases for details

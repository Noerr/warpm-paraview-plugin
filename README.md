# WARPM ParaView Reader Plugin

A C++ ParaView reader plugin that loads native [WARPM](https://github.com/Noerr/warpm) HDF5 simulation output with high-order Lagrange element support.

## Features

- **Direct HDF5 loading** - Open `.h5` files directly, no post-processing format conversion needed
- **High-order elements** - DG elements rendered as VTK Lagrange cells (not broken into linear sub-elements)
- **Time series support** - File groups, multi-select, or `.warpm` metadata files
- **Physical time values** - Reads actual simulation time from each frame
- **Variable selection** - Choose which fields to load via Properties panel
- **Geometry caching** - Optimized frame switching (only reloads field data)
- **Arbitrary dimensions** - Supports 1D, 2D, and 3D physical/velocity spaces
- **Phase space visualization** - Dedicated reader for phase space data (probability distribution functions) in up to six dimensions (3D+3V)

## Two Reader Classes

The plugin provides two readers that automatically detect which file type they can handle:

### WARPM HDF5 Reader (`vtkWARPMReader`)

For standard physical-space simulation data (electromagnetic fields, etc.):
- Single output: unstructured grid with field variables as point data
- Supports 1D, 2D, or 3D meshes
- Automatically selected for files without phase space variables

### WARPM Phase Space Reader (`vtkWARPMPhaseSpaceReader`)

For kinetic simulations with probability distribution functions:
- **Three output ports:**
  - **Port 0 - Physical Space:** Physical mesh with EM field variables
  - **Port 1 - Velocity Space:** Velocity mesh showing distribution function at a selected physical node location
  - **Port 2 - Slice Position:** Single point marker showing where the velocity slice is taken
- Supports arbitrary dimension combinations (1-3D physical × 1-3D velocity)
- Slice control indices in Properties panel select which physical cell and node to visualize

**Automatic reader selection:** When opening a WARPM file, ParaView automatically chooses the correct reader. The phase space reader only claims files that contain phase-space variables (identified by included coordinate names starting with 'v' or 'w'). This prevents the "choose reader" dialog from appearing.

### Working with Multiple Output Ports

The Phase Space Reader produces three outputs accessible in the Pipeline Browser:

1. After loading a phase space file, expand the reader in the Pipeline Browser
2. Select which output to view:
   - "Physical Space" - Standard field visualization
   - "Velocity Space" - Distribution function f(vx, vy) or f(vx, vy, vz)
   - "Slice Position" - Marker point (useful for overlaying on Physical Space)
3. Use the **Physical Slice Indices** and **Physical Node Index** properties to change where the velocity slice is taken

## Requirements

- **ParaView 6.0** built from source (developer build required)
- CMake 3.16+
- HDF5 library
- Qt6

> **Note:** The plugin does not work with pre-built MacOS ParaView binaries due to C++ symbol visibility issues. See [Known Limitations](#known-limitations).

## Building

```bash
# Clone the repository
git clone https://github.com/Noerr/warpm-paraview-plugin.git
cd warpm-paraview-plugin

# Build (adjust ParaView_DIR to your build location)
mkdir build && cd build
cmake .. -DParaView_DIR=$HOME/path/to/paraview/build
make -j4
```

Output: `build/lib/paraview-6.0/plugins/WARPMReader/WARPMReader.so`

## Usage

### Loading the Plugin

1. **GUI:** Tools > Manage Plugins > Load New > Browse to `WARPMReader.so`
2. **pvpython:**
   ```python
   import paraview.simple as pv
   pv.LoadPlugin('/path/to/WARPMReader.so', remote=False)
   ```

### Opening WARPM Data

**Direct .h5 files:**
- File > Open, select one or more `.h5` files
- The appropriate reader is selected automatically based on file contents
- For file groups (e.g., `simulation_*.h5`), ParaView handles time series automatically

**Using .warpm metadata files:**
```
# my_simulation.warpm
simulation_*.h5
```

### pvpython Examples

**Physical-space data:**
```python
import paraview.simple as pv
pv.LoadPlugin('/path/to/WARPMReader.so', remote=False)

reader = pv.WARPMReader(FileName='maxwell_2d_0000.h5')
reader.UpdatePipeline()
print(f"Points: {reader.GetDataInformation().GetNumberOfPoints()}")
```

**Phase space data:**
```python
reader = pv.WARPMPhaseSpaceReader(FileName='vlasov_2d2v_0000.h5')
reader.PhysicalSliceIndices = [2, 3, 0]  # Slice at cell (2,3) in physical space
reader.PhysicalNodeIndex = 4             # Use node 4 within that cell
reader.UpdatePipeline()

# Access different output ports
physical_output = pv.servermanager.Fetch(reader, idx=0)  # Physical mesh
velocity_output = pv.servermanager.Fetch(reader, idx=1)  # Velocity mesh
probe_output = pv.servermanager.Fetch(reader, idx=2)     # Slice position marker
```

## WARPM HDF5 Format

The reader expects HDF5 files with this structure:

```
/domains/<domain_name>/
  - ndims: 1, 2, or 3 (physical) or 2-6 (phase space)
  - numCells: [Nx, ...] or [Nx, Ny, Nvx, Nvy, ...]
  - CoordinateNames: ["x", "y"] or ["x", "y", "wx", "wy"]
  - VertexCoordinateExpressions: parametric mesh definition

/variables/<name>/
  - OnDomain: "<domain_name>"
  - ElementOrder: [px, py, ...]
  - ComponentNames: ['Ex', 'Ey', ...] or ['pdf']
  - /data: field values

/timeData/
  - time: physical simulation time
```

**Domain type detection:** Domains are classified by their `CoordinateNames` attribute:
- Physical coordinates: `x`, `y`, `z`
- Velocity coordinates: names starting with `v` or `w` (e.g., `vx`, `wy`)

## Known Limitations

### Requires Developer-Built ParaView

The plugin does not work with pre-built ParaView.app from Kitware due to C++ symbol visibility issues (`-fvisibility=hidden`). This causes `dynamic_cast` failures when loading the plugin.

**Workaround:** Build ParaView from source with default visibility settings.

See [ParaView GitLab Issue #14869](https://gitlab.kitware.com/paraview/paraview/-/issues/14869) for details.

### Single Domain Per Type

If a file contains multiple physical domains or multiple phase space domains, only variables from the first domain of each type are loaded. A warning is displayed for skipped variables.

## License

BSD 3-Clause License. See [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! Please open an issue or pull request.

## Acknowledgments

- [ParaView](https://www.paraview.org/) and [VTK](https://vtk.org/) by Kitware
- Development assisted by [Claude Code](https://claude.ai/code)

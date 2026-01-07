# WARPM ParaView Reader Plugin

A C++ ParaView reader plugin that loads native [WARPM](https://github.com/Noerr/warpm) HDF5 simulation output with high-order Lagrange element support.

## Features

- **Direct HDF5 loading** - Open `.h5` files directly, no format conversion needed
- **High-order elements** - DG elements rendered as VTK Lagrange cells (not broken into linear sub-elements)
- **Time series support** - File groups, multi-select, or `.warpm` metadata files
- **Physical time values** - Reads actual simulation time from each frame
- **Variable selection** - Choose which fields to load via Properties panel
- **Geometry caching** - Optimized frame switching (only reloads field data)

## Requirements

- **ParaView 6.0** built from source (developer build required)
- CMake 3.16+
- HDF5 library
- Qt6

> **Note:** The plugin does not work with pre-built ParaView binaries due to C++ symbol visibility issues. See [Known Limitations](#known-limitations).

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
- Choose "WARPM HDF5 Data" as the reader
- For file groups (e.g., `simulation_..h5`), ParaView handles time series automatically

**Using .warpm metadata files:**
```
# my_simulation.warpm
simulation_*.h5
```

## WARPM HDF5 Format

The reader expects HDF5 files with this structure:

```
/domains/physical_space_domain/
  - numCells: [Nx, Ny]
  - VertexCoordinateExpressions: mesh definition

/variables/<name>/
  - ElementOrder: [px, py]
  - ComponentNames: ['Ex', 'Ey', ...]
  - /data: field values

/timeData/
  - time: physical simulation time
```

## Known Limitations

### Requires Developer-Built ParaView

The plugin does not work with pre-built ParaView.app from Kitware due to C++ symbol visibility issues (`-fvisibility=hidden`). This causes `dynamic_cast` failures when loading the plugin.

**Workaround:** Build ParaView from source with default visibility settings.

See [ParaView GitLab Issue #14869](https://gitlab.kitware.com/paraview/paraview/-/issues/14869) for details.

## License

BSD 3-Clause License. See [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! Please open an issue or pull request.

## Acknowledgments

- [ParaView](https://www.paraview.org/) and [VTK](https://vtk.org/) by Kitware
- Development assisted by [Claude Code](https://claude.ai/code)

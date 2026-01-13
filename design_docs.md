# WARPM ParaView Plugin - Architecture Documentation

This document provides detailed architecture diagrams and design documentation for the WARPM ParaView Reader Plugin.

For general usage and build instructions, see [README.md](README.md).
For development guidelines and implementation details, see [CLAUDE.md](CLAUDE.md).

## Architecture Diagrams

### Class Hierarchy

```
VTK/ParaView Base Classes                    WARPM Plugin Classes
========================                     ====================

vtkObjectBase
    │
    └── vtkObject
            │
            └── vtkAlgorithm
                    │
                    ├── vtkDataObjectAlgorithm
                    │       │
                    │       └── vtkFileSeriesReader ◄─────────┐
                    │               │                         │ inherits
                    │               └── vtkWARPMPhaseSpaceFileSeriesReader
                    │                   (just overrides constructor
                    │                    to set 3 output ports)
                    │
                    └── vtkUnstructuredGridAlgorithm
                            │
                            └── vtkWARPMReader (base reader)
                                    │  - 1 output port
                                    │  - Physical-space data only
                                    │  - All static helpers defined here
                                    │
                                    └── vtkWARPMPhaseSpaceReader
                                        - 3 output ports
                                        - Overrides CanReadFile, RequestInfo, RequestData
                                        - Phase space + physical data
```

### ServerManager Proxy Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        WARPMReader.xml                                  │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ProxyGroup: "internal_sources" (not user-visible)                      │
│  ┌─────────────────────────────────────────────────────────────────┐    │
│  │  WARPMReaderCore (class=vtkWARPMReader)                         │    │
│  │  WARPMPhaseSpaceReaderCore (class=vtkWARPMPhaseSpaceReader)     │    │
│  └─────────────────────────────────────────────────────────────────┘    │
│                              ▲                                          │
│                              │ SubProxy reference                       │
│                              │                                          │
│  ProxyGroup: "sources" (user-visible in ParaView)                       │
│  ┌─────────────────────────────────────────────────────────────────┐    │
│  │  WARPMReader                                                    │    │
│  │    class=vtkFileSeriesReader                                    │    │
│  │    si_class=vtkSIMetaReaderProxy                                │    │
│  │    → wraps WARPMReaderCore for multi-file time series           │    │
│  ├─────────────────────────────────────────────────────────────────┤    │
│  │  WARPMPhaseSpaceReader                                          │    │
│  │    class=vtkWARPMPhaseSpaceFileSeriesReader                     │    │
│  │    si_class=vtkSIMetaReaderProxy                                │    │
│  │    → wraps WARPMPhaseSpaceReaderCore (3 output ports)           │    │
│  └─────────────────────────────────────────────────────────────────┘    │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### VTK Pipeline Lifecycle

When ParaView interacts with a reader, methods are called in this order:

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     READER INSTANTIATION                                │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  1. User: File → Open → selects .h5 or .warpm file(s)                   │
│                                                                         │
│  2. ParaView queries all registered readers:                            │
│     ┌──────────────────────────────────────────────────────────────┐    │
│     │  CanReadFile(filename)  [STATIC - no instance yet]           │    │
│     │  • vtkWARPMReader::CanReadFile()                             │    │
│     │    - Calls ResolveWARPMFile() to validate HDF5 structure     │    │
│     │    - Returns 1 if valid WARPM file WITHOUT phase space vars  │    │
│     │  • vtkWARPMPhaseSpaceReader::CanReadFile()                   │    │
│     │    - Returns 1 ONLY if file HAS phase space variables        │    │
│     │  → Exactly one reader claims the file (no dialog needed)     │    │
│     └──────────────────────────────────────────────────────────────┘    │
│                                                                         │
│  3. ParaView instantiates the winning reader:                           │
│     ┌──────────────────────────────────────────────────────────────┐    │
│     │  vtkWARPMReader::New() or vtkWARPMPhaseSpaceReader::New()    │    │
│     │  • Constructor initializes member variables                  │    │
│     │  • Sets number of output ports (1 or 3)                      │    │
│     │  • Creates vtkDataArraySelection for variable checkboxes     │    │
│     └──────────────────────────────────────────────────────────────┘    │
│                                                                         │
│  4. If file series, ParaView wraps with FileSeriesReader:               │
│     ┌──────────────────────────────────────────────────────────────┐    │
│     │  vtkFileSeriesReader or vtkWARPMPhaseSpaceFileSeriesReader   │    │
│     │  • SetReader(internal_reader) - stores our reader instance   │    │
│     │  • AddFileName() called for each selected file               │    │
│     └──────────────────────────────────────────────────────────────┘    │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────┐
│                     PIPELINE INFORMATION PHASE                          │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  5. ParaView requests metadata (before any data is loaded):             │
│     ┌──────────────────────────────────────────────────────────────┐    │
│     │  RequestInformation(request, inputVec, outputVec)            │    │
│     │                                                              │    │
│     │  Base reader (vtkWARPMReader):                               │    │
│     │  • ParseWarpmFile() - expands glob patterns, builds file list│    │
│     │  • ReadTimeValues() - opens each .h5, reads /timeData/time   │    │
│     │  • Populates PointDataArraySelection with variable names     │    │
│     │  • Sets TIME_STEPS and TIME_RANGE on output info             │    │
│     │                                                              │    │
│     │  Phase space reader (vtkWARPMPhaseSpaceReader):              │    │
│     │  • Calls Superclass::RequestInformation() first              │    │
│     │  • Detects phase space domains via DetectPhaseSpaceDomain()  │    │
│     │  • Adds phase space variables to PointDataArraySelection     │    │
│     │  • Sets TIME_STEPS on ALL 3 output ports                     │    │
│     └──────────────────────────────────────────────────────────────┘    │
│                                                                         │
│  6. ParaView shows Properties panel with:                               │
│     • Variable checkboxes (from PointDataArraySelection)                │
│     • Time slider (from TIME_STEPS)                                     │
│     • Slice controls (phase space reader only)                          │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────┐
│                     PIPELINE DATA PHASE                                 │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  7. User clicks "Apply" or changes time step:                           │
│     ┌──────────────────────────────────────────────────────────────┐    │
│     │  RequestData(request, inputVec, outputVec)                   │    │
│     │                                                              │    │
│     │  • FindClosestFrame(requestedTime) - maps time to file index │    │
│     │  • Opens HDF5 file for that frame                            │    │
│     │  • Checks geometry cache (reuses if mesh unchanged)          │    │
│     │                                                              │    │
│     │  If cache miss (first load or mesh changed):                 │    │
│     │  • ComputeVertexPositions() - evaluates coordinate exprs     │    │
│     │  • GetGLLNodes() - gets Gauss-Lobatto-Legendre positions     │    │
│     │  • BuildLagrangeMesh() - creates VTK points and cells        │    │
│     │  • BuildWarpmToVTKMapping() - computes node reordering       │    │
│     │                                                              │    │
│     │  For each enabled variable:                                  │    │
│     │  • H5Dread() - reads field data from HDF5                    │    │
│     │  • Reorders nodes from WARPM to VTK Lagrange ordering        │    │
│     │  • Creates vtkDoubleArray, adds to output PointData          │    │
│     │                                                              │    │
│     │  Phase space reader additionally:                            │    │
│     │  • Builds separate velocity mesh for Port 1                  │    │
│     │  • Uses H5 hyperslab to read only sliced physical location   │    │
│     │  • Creates probe point for Port 2                            │    │
│     └──────────────────────────────────────────────────────────────┘    │
│                                                                         │
│  8. Output data flows downstream to filters/renderers                   │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────┐
│                     SUBSEQUENT UPDATES                                  │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Time step change:                                                      │
│  • RequestData() called again with new UPDATE_TIME_STEP                 │
│  • Geometry cache HIT (mesh doesn't change between frames)              │
│  • Only field data is re-read from new HDF5 file                        │
│                                                                         │
│  Property change (e.g., slice indices):                                 │
│  • RequestData() called again                                           │
│  • May rebuild velocity mesh if slice location changed                  │
│                                                                         │
│  Variable selection change:                                             │
│  • RequestData() called again                                           │
│  • Geometry cache HIT                                                   │
│  • Only reads newly-enabled variables                                   │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### Static Method Callers

**Public static methods** (called by ParaView, no instance required):

```
┌────────────────────────────────────────────────────────────────────────┐
│  PUBLIC STATIC METHOD               │  CALLED BY                       │
├─────────────────────────────────────┼──────────────────────────────────┤
│  vtkWARPMReader::CanReadFile(fname) │  ParaView ReaderFactory          │
│    • Called during File→Open to determine which reader handles file    │
│    • Returns 1 for WARPM files WITHOUT phase space variables           │
│    • Must not instantiate reader (static context)                      │
├─────────────────────────────────────┼──────────────────────────────────┤
│  vtkWARPMPhaseSpaceReader::CanReadFile(fname)                          │
│    • Called during File→Open                                           │
│    • Returns 1 ONLY for files WITH phase space variables               │
│    • Ensures exactly one reader claims each file type                  │
└─────────────────────────────────────┴──────────────────────────────────┘
```

**Protected static methods** (internal helpers, called from within reader classes):

```
┌────────────────────────────────────────────────────────────────────────┐
│  PROTECTED STATIC METHOD            │  CALLED BY                       │
├─────────────────────────────────────┼──────────────────────────────────┤
│  File Validation (used by CanReadFile):                                │
│  ResolveWARPMFile(fname)            │  CanReadFile                     │
│  IsWarpmHDF5File(fname)             │  ResolveWARPMFile                │
│  HasPhaseSpaceVariables(fname)      │  CanReadFile (phase space)       │
├─────────────────────────────────────┼──────────────────────────────────┤
│  Domain Analysis (used by RequestInformation):                         │
│  IsPhysicalOnlyDomain(domGroup)     │  RequestInfo, CanReadFile        │
│  IsVariableOnPhysicalDomain(...)    │  RequestInfo                     │
│  DetectPhaseSpaceDomain(...)        │  RequestInfo                     │
├─────────────────────────────────────┼──────────────────────────────────┤
│  HDF5 Attribute Readers (used throughout):                             │
│  ReadIntAttribute(loc, name, val)   │  RequestInfo, RequestData        │
│  ReadDoubleAttribute(...)           │  RequestInfo, RequestData        │
│  ReadStringAttribute(...)           │  RequestInfo, RequestData        │
│  ReadIntArrayAttribute(...)         │  RequestInfo, RequestData        │
│  ReadStringArrayAttribute(...)      │  RequestInfo, RequestData        │
├─────────────────────────────────────┼──────────────────────────────────┤
│  Mesh Building (used by RequestData):                                  │
│  ComputeVertexPositions(expr,...)   │  RequestData                     │
│  GetGLLNodes(order)                 │  RequestData                     │
│  BuildLagrangeMesh(...)             │  RequestData                     │
├─────────────────────────────────────┼──────────────────────────────────┤
│  Node Ordering (used by RequestData):                                  │
│  BuildWarpmToVTKMapping(orders)     │  RequestData                     │
│    └── BuildWarpmToVTKMapping1D()   │    (dispatched by dimension)     │
│    └── BuildWarpmToVTKMapping2D()   │                                  │
│    └── BuildWarpmToVTKMapping3D()   │                                  │
├─────────────────────────────────────┼──────────────────────────────────┤
│  N-Dimensional Helpers (used by RequestData):                          │
│  DecomposeIndex(flat, sizes)        │  RequestData (loop iteration)    │
│  ComputeFlatIndex(indices, sizes)   │  RequestData (data indexing)     │
│  ComputeTotalNodes(orders)          │  RequestData (allocation sizing) │
│  GetLagrangeCellType(ndims)         │  RequestData (VTK cell type)     │
└─────────────────────────────────────┴──────────────────────────────────┘
```

### Object Lifetime

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     PARAVIEW SESSION LIFETIME                           │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Plugin Load (Tools → Manage Plugins → Load)                            │
│  ├── Plugin .so loaded into memory                                      │
│  ├── ServerManager XML parsed, proxies registered                       │
│  └── ReaderFactory updated with file extensions                         │
│                                                                         │
│  File Open                                                              │
│  ├── CanReadFile() called (static, no instance)                         │
│  ├── Reader instance created ─────────────────────┐                     │
│  │   └── Constructor called                       │                     │
│  │       └── Output ports configured              │                     │
│  │       └── PointDataArraySelection created      │                     │
│  │                                                │                     │
│  │   FileSeriesReader wrapper created ◄───────────┤ (if multi-file)     │
│  │   └── Holds reference to our reader            │                     │
│  │                                                │                     │
│  ├── RequestInformation() called                  │                     │
│  │   └── File list populated                      │  Instance           │
│  │   └── Time values read                         │  Lifetime           │
│  │   └── Variables enumerated                     │                     │
│  │                                                │                     │
│  ├── User clicks Apply                            │                     │
│  │   └── RequestData() called                     │                     │
│  │       └── Geometry built (cached)              │                     │
│  │       └── Field data loaded                    │                     │
│  │                                                │                     │
│  ├── Time animation / property changes            │                     │
│  │   └── RequestData() called repeatedly          │                     │
│  │       └── Geometry cache reused                │                     │
│  │       └── Only field data changes              │                     │
│  │                                                │                     │
│  └── Pipeline deleted (user deletes source)       │                     │
│      └── Destructor called ◄──────────────────────┘                     │
│          └── Cached geometry released                                   │
│          └── PointDataArraySelection released                           │
│                                                                         │
│  Plugin Unload (rare, usually stays loaded)                             │
│  └── Plugin .so unloaded from memory                                    │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### Instance Count Clarification

For a time series (multiple .h5 files or .warpm with glob pattern):
- **ONE reader instance** (vtkWARPMReader or vtkWARPMPhaseSpaceReader)
- **ONE FileSeriesReader wrapper** holds reference to that reader
- The wrapper calls `SetFileName()` on our reader for each time step
- Our reader's `FrameFiles` vector holds all file paths
- Geometry cache is effective: built once, reused across all frames

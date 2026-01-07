# README.md

** WARPM ParaView Reader Plugin ** 
The goal of this project is to create a new ParaView Reader plugin that allows ParaView to read and work with datasets in the natural output format of WARPM.  (one .h5 HDF5 output file per time frame)  This avoids the expensive and time-consuming conversion of WARPM simulation output into an alternative format that ParaView natively understands.

This is a new project, with the only existing material being the old (working) scripts in Python that convert from native WARPM .h5 output to VTK .vts files.
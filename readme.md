Software Requirements:
======================

CMake [cmake.org](https://cmake.org/):
* minimum version 3.13
    
Paraview [www.paraview.org](https://www.paraview.org/):
* any version, for example the most recent
    
C++ Compiler:
* with C++11 support, for example gcc6.3 or Visual C++ 14.0
    
How to get VirtualFluidsCPU:
==========================

Option 1: use git
1. checkout out https://git.irmb.bau.tu-bs.de/VirtualFluids/VirtualFluidsCPU.git with your credentials

Option 2: without git
1. go to git.irmb.tu-bs.de
2. Log in with your credentials
3. click on VirtualFluids/VirtualFluidsCPU
4. click on the download symbol on the top right and download zip/tar.gz file

How to build VirtualFluidsCPU:
============================

1. CMake the project
2. set the output path in Applications/LidDrivenCavity/LidDrivenCavity.cpp
3. build the project ("compile")
4. run the generated executable (usually in <build directory>/Applications/LidDrivenCavity)

VirtualFluidsCPU results files:
===============================

VirtualFluidsCPU generates a set of multiple output directories in the prescribed output path. The flow fields can be found in the _mq_ directory. To view the flow fields, it is most conveniant to open the _mq_collection.pvd_ file in Paraview. The _bc_ directory contains the boundary condition information, the _geo_ directory contains information on the geometry of the flow domain and the _blocks_ directory contains the block grid.

Documentation:
==============

The doxygen generated documentation can be found [here](https://git.irmb.bau.tu-bs.de/doku/CPU).

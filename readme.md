Software Requirements:
======================

CMake [cmake.org](https://cmake.org/):
* minimum version 3.13

CUDA [developer.nvidia.com/cuda-zone](https://developer.nvidia.com/cuda-zone):
* Minimum CUDA Version 9.0
* Minimum Compute Capability 3.0, because of maximal number of Blocks in x direction
* Recommended Compute Capability 6.0, because of atomics for double precision floating point data (GKS only)
    
Paraview [www.paraview.org](https://www.paraview.org/):
* any version, for example the most recent
    
C++ Compiler:
* with C++11 support, for example gcc6.3 or Visual C++ 14.0
    
How to get VirtualFluidsGPU:
==========================

Option 1: use git
1. checkout out https://git.irmb.bau.tu-bs.de/VirtualFluids/VirtualFluidsGPU.git with your credentials

Option 2: without git
1. go to git.irmb.tu-bs.de
2. Log in with your credentials
3. click on VirtualFluids/VirtualFluidsGPU
4. click on the download symbol on the top right and download zip/tar.gz file

How to build VirtualFluidsGPU:
============================

1. CMake the project
2. set the output path in targets/apps/LidDrivenCavity/LidDrivenCavity.cpp
3. build the project ("compile")
4. run the generated executable (usually in <build directory>/bin/)

VirtualFluidsGPU results files:
===============================

VirtualFluidsGPU generates a the time series of output files directly in the output path. In Paraview these time series can be read directly.
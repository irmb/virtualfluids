Software Requirements:
======================

CMake [cmake.org](https://cmake.org/):
* minimum version 3.13
    
Paraview [www.paraview.org](https://www.paraview.org/):
* any version, for example the most recent
    
C++ Compiler:
* with C++11 support, for example gcc6.3 or Visual C++ 14.0
    
How to get Virtual Fluids:
==========================

Option 1: use git
1. checkout out https://git.irmb.bau.tu-bs.de/VirtualFluids/VirtualFluidsCPU.git with your credentials

Option 2: without git
1. go to git.irmb.tu-bs.de
2. Log in with your credentials
3. click on VirtualFluids/VirtualFluidsCPU
4. click on the download symbol on the top right and download zip/tar.gz file

How to build Virtual Fluids:
============================

1. CMake the project
2. set the output path in source/Applications/LidDrivenCavity/LidDrivenCavity.cpp
3. build the project ("compile")
4. run the generated executable (usually in <build directory>/Applications/LidDrivenCavity)

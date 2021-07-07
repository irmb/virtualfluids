![VirtualFluids](docs/img/VF_logo.png)

VirtualFluids (VF) is a research code developed at the Institute for Computational Modeling in Civil Engineering (iRMB). The code is a Computational Fluid Dynamics (CFD) solver based on the Lattice Boltzmann Method (LBM) for turbulent, thermal, multiphase and multicomponent flow problems as well as for multi-field problems such as Fluid-Structure-interaction including distributed pre- and postprocessing capabilities for simulations with more than 100 billion degrees of freedom.

## Getting Start
### Suported Platforms
VirtualFluids has been used on a variety of platforms:
 - Linux
 - Mac OS X
 - Windows
 - Cygwin
### Software Requirements
 
 - [CMake](https://cmake.org/) (minimum version 3.15)
 - C++ compiler with C++14 support
 - [Paraview](https://www.paraview.org/) (most recent version)

with usage of the GPU:  
 - CUDA [developer.nvidia.com/cuda-zone](https://developer.nvidia.com/cuda-zone):
    * Minimum CUDA Version 9.0
    * Minimum Compute Capability 3.0, because of maximal number of Blocks in x direction
    * Recommended Compute Capability 6.0, because of atomics for double precision floating point data (GKS only)
    

### Contributing
To contribute to VirtualFluids please follow these [instructions](CONTRIBUTING.md).

### Build VirtualFluids
```shell
$ mkdir build
$ cd build
```
Pass the relevant [options](#options) to cmake.
E.g. for the CPU part:
```shell
$ cmake .. -DBUILD_VF_CPU=ON
$ make
```
Alternatively enable the options via the cmake-gui.

### <a id="options"></a> Options
- BUILD_VF_CPU
  - Build VirtualFluids CPU variant
- BUILD_VF_GPU 
  - Build VirtualFluids GPU variant
- BUILD_VF_UNIT_TESTS
  -  Build VirtualFluids unit tests
- VF_DOUBLE_ACCURACY 
    - GPU change between Double and Single Precision

### Result Files
The output files can be found in `<build directory>/bin/output`. As there is an usually high amount of data, you might want to change the output path in the main function.

The CPU part generates a set of multiple output directories in the prescribed output path. The flow fields can be found in the _mq_ directory. To view the flow fields, it is most conveniant to open the _mq_collection.pvd_ file in Paraview. The _bc_ directory contains the boundary condition information, the _geo_ directory contains information on the geometry of the flow domain and the _blocks_ directory contains the block grid.

A GPU computation generates a the time series of output files directly in the output path. In Paraview these time series can be read directly.



## Documentation
The doxygen generated documentation can be found [here](https://irmb.gitlab-pages.rz.tu-bs.de/VirtualFluids_dev).


## Known Issues
If you notice any problems on your platform, please report an issue. 


## Authors
A list of the developers of VirtualFluids is available [here](AUTHORS.md).

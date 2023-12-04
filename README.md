![VirtualFluids](docs/img/VF_logo.png)

VirtualFluids (VF) is a research code developed at the Institute for Computational Modeling in Civil Engineering (iRMB). The code is a Computational Fluid Dynamics (CFD) solver based on the Lattice Boltzmann Method (LBM) for turbulent, thermal, multiphase and multicomponent flow problems as well as for multi-field problems such as Fluid-Structure-interaction including distributed pre- and postprocessing capabilities for simulations with more than 100 billion degrees of freedom.

## Getting Started
VirtualFluids is mainly supported on these two platforms:
 - Linux
 - Windows

VirtualFluids can also be build and used in a Docker image. An ubuntu development environment is located in the [container registry](https://git.rz.tu-bs.de/irmb/virtualfluids/container_registry).
An extensive guide about the usage and development in VirtualFluids with docker can be found [here](https://git.rz.tu-bs.de/irmb/virtualfluids/-/wikis/Getting-Started-with-the-Development-of-VirtualFluids).


The following is a brief explanation of how to use it without Docker:
### Software Requirements

 - [CMake](https://cmake.org/) (minimum version 3.15)
 - C++ compiler with C++14 support
 - [Paraview](https://www.paraview.org/) for visualizations (most recent version)


with usage of the GPU:
 - CUDA [developer.nvidia.com/cuda-zone](https://developer.nvidia.com/cuda-zone):
    * Minimum CUDA Version 9.0
    * Minimum Compute Capability 3.0, because of maximal number of Blocks in x direction


### Build VirtualFluids
```shell
$ mkdir build && cd build
```
Pass the relevant [options](#options) to cmake.
E.g. for the CPU part:
```shell
$ cmake .. -DVF_ENABLE_CPU=ON
$ make
```
Alternatively enable the options via the cmake-gui.

### <a id="options"></a> Options
- VF_ENABLE_CPU
  - Build VirtualFluids CPU variant
- VF_ENABLE_GPU
  - Build VirtualFluids GPU variant
- VF_ENABLE_UNIT_TESTS
  -  Build VirtualFluids unit tests
- VF_ENABLE_DOUBLE_ACCURACY
    - GPU change between Double and Single Precision

### Result Files
The output files can be found in `<build directory>/bin/output`. As there is an usually high amount of data, you might want to change the output path in the main function.

The CPU part generates a set of multiple output directories in the prescribed output path. The flow fields can be found in the _mq_ directory. To view the flow fields, it is most conveniant to open the _mq_collection.pvd_ file in Paraview. The _bc_ directory contains the boundary condition information, the _geo_ directory contains information on the geometry of the flow domain and the _blocks_ directory contains the block grid.

A GPU computation generates a the time series of output files directly in the output path. In Paraview these time series can be read directly.

## Contributing
To contribute to VirtualFluids please follow these [instructions](CONTRIBUTING.md).

## Documentation
The doxygen generated documentation can be found [here](https://irmb.gitlab-pages.rz.tu-bs.de/VirtualFluids_dev).


## Known Issues
If you notice any problems on your platform, please report an [issue](https://git.rz.tu-bs.de/irmb/virtualfluids/-/issues/new).


## Authors
A list of the developers of VirtualFluids is available [here](AUTHORS.md).

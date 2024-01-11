<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# Build and Run

This guide describes how to start using and developing VirtualFluids in the terminal. Alternativly you can use Visual Studio to build an run it.

## Build

The necessary packages needs to be installed. Either we are working in a <!-- DOXYGEN_MAKE_REF -->[Docker container](Getting-Started-with-Docker.md) or have installed the packages <!-- DOXYGEN_MAKE_REF -->[manually](Getting-Started-Not-Using-Docker.md). 
We can check if the necessary packages are installed by running the following commands in the terminal:
```
   git --version
   cmake --version
   mpic++ --version
   g++ --version
```

Then we can clone the repository onto the machine:
```
   git clone https://git.rz.tu-bs.de/irmb/virtualfluids.git
```

Naviagting to the project folder:
```
   cd virtualfluids
```

and build the code
```
   mkdir build && cd build
   cmake --preset=all_make -DCMAKE_CUDA_ARCHITECTURE=70 ..
   make -j 8
```

The option "-DCMAKE_CUDA_ARCHITECTURE" is only necessary when the GPU version is used (CMAKE_CUDA_ARCHITECTURE should correspond to the [compute capability](https://en.wikipedia.org/wiki/CUDA#GPUs_supported) of your GPU).
VirtualFluids cmake provides the following presets:

- all_make: build all targets with Make
- cpu_make: build only the CPU targets with Make
- gpu_make: build only the GPU targets with Make
- all_msvc: build all targets with MSVC
- cpu_msvc: build only the CPU targets with MSVC
- gpu_msvc: build only the GPU targets with MSVC

Additionaly, the following options can be passed to cmake:
- -DVF_ENABLE_UNIT_TESTS=ON: enable unit tests (included in all_make and all_msvc)
- -DVF_ENABLE_DOUBLE_ACCURACY=ON: enable double precision (included in all_make, all_msvc, cpu_make and cpu_msvc)
- -DVF_ENABLE_PYTHON_BINDINGS=ON: enable python bindings (included in all_make and all_msvc)

## Run the examples

VirtualFluids project comes with a list of examples. The source code of the examples are located in the folder `./apps/`. Most of the apps requires a configuration file, which lays next to source code example.

For instance a simulation on the GPU containing a flow around a sphere can be started with the following command:
```
   ./build/bin/SphereInChannel ./apps/gpu/SphereInChannel/sphere_1level.cfg

```

The result files of this simulation are usually stored in: `./output/`
The result files of VirtualFluids are mostly in [VTK](https://kitware.github.io/vtk-examples/site/VTKFileFormats/) format. These files can be visualised with the free software [Paraview](https://www.paraview.org/).

The CPU part generates a set of multiple output directories in the prescribed output path. The flow fields can be found in the _mq_ directory. To view the flow fields, it is most conveniant to open the _mq_collection.pvd_ file in Paraview. The _bc_ directory contains the boundary condition information, the _geo_ directory contains information on the geometry of the flow domain and the _blocks_ directory contains the block grid.

A GPU computation generates a the time series of output files directly in the output path. In Paraview these time series can be read directly.
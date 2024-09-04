<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->
# Project Structure

The general code structure of the project is as follows:

## Folder

- `3rdparty/`: Contains the third party libraries: Metis, cuda_samples and MuParser
    - CUDA installation is required to build the project with gpu support.
    - MPI and OpenMP are required to build the project with cpu support.
    - Other third party libraries (spdlog, boost, yaml-cpp, googletest) are downloaded on demand during the cmake configuration process.
- `apps/`: Contains the source code of the examples.
- `CMake/`: Contains the cmake scripts.
- `Containers/`: Contains the Dockerfiles
- `docs/`: Contains the documentation.
- `Python/`: Contains the python applications
- `pythonbindings/`: Contains the python bindings.
- `src/`: Contains the source code of the library.
- `tests/`: Contains the unit, regression and performance tests


## Libraries
The source code in /src is divided into several libraries:
- `basics`
- `core/gpu`
- `core/cpu`
- `lbm`
- `logger`
- `parallel`

`core/cpu` and `core/gpu` contain the main data structures for the CPU and GPU parallelization, respectively. `basics` contains some functions used by the other libraries, like geometries or vtk writer `lbm` contains the Lattice Boltzmann Method specific functions. `logger` contains the logging functions. `parallel` contains a MPI abstraction layer.


## CMake Build-System

The cmake configuration is divided into several files:
- `CMakeLists.txt`: The main cmake file. it contains mainly project options. It includes the common subdirectories (basics, parallel, lbm, logger) and either `gpu.cmake` or `cpu.cmake` depending on the build type.
- `gpu.cmake` and `cpu.cmake` are including the `core/gpu` or `core/cpu` subdirectories, respectively and their specific applications.
- every library contains its own cmake file e.g. `src/lbm/CMakeLists.txt`. In this file usually a macro `vf_add_library()`is called.
- The cmake structure contains two interface projects: `project_warnings` and `project_options`: These projects are included in every library and application. They contain compiler options and warnings.
- The `CMake/` directory contains our custom cmake macros and functions. During the configuration process, for every library a `vf_add_library()` macro is called. This macro is defined in `CMake/VirtualFluidsMacros.cmake`. This macro is responsible for adding the library to the project, setting the include directories, linking the libraries and setting the compiler options. It also traverses the current folder and alle subfolders recursively to find all source files and add them to the library.
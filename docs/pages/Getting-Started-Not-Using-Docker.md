<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# Getting Started Not Using Docker

VirtualFluids is supported on Linux and Windows.
We recommend using VirtualFludis with Docker. If, for whatever reason, Docker cannot be used, VirtualFluids can also be installed manually. The following packages are necessary for this:

- [Git](https://git-scm.com/) is necessary to get access to VirtualFluids. Therefore at least some basic knowledge in Git is required. If you are new to Git, this [free book](https://git-scm.com/book/en/v2) is a good starting point
- [CMake](https://cmake.org/) (minimum version 3.20)
- C++ compiler with C++17 support
    - [GCC](https://gcc.gnu.org/)
    - or [Clang](https://clang.llvm.org/)
    - or [MSVC](https://visualstudio.microsoft.com/de/vs/features/cplusplus/) (only on Windows)
- MPI
    - [OpenMPI](https://www.open-mpi.org/)
    - or [MPICH](https://www.mpich.org/)
    - or [Microsoft MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi) (only on Windows)
- [Paraview](https://www.paraview.org/) for visualizations (most recent version)
 - with usage of the GPU: CUDA [developer.nvidia.com/cuda-zone](https://developer.nvidia.com/cuda-zone):
    * Minimum CUDA Version 9.0
    * Minimum Compute Capability 3.0, because of maximal number of Blocks in x direction
- a development environment:
    * [Visual Studio](https://visualstudio.microsoft.com/de/vs/features/cplusplus/) (only on Windows)
    * or [Clion](https://www.jetbrains.com/de-de/clion/)
    * or [VS Code](https://code.visualstudio.com/)

After the installation of the necessary software, the VirtualFluids project folder can be either either directly opend in vscode or clion or a Visual Studio solution can be generated with cmake.


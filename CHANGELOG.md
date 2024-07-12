<!-- SPDX-License-Identifier: CC-BY-4.0 -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->
# Changelog

## [0.1.2](https://git.rz.tu-bs.de/irmb/VirtualFluids/-/milestones/3) - 2024-07-09

### Added
- [ALL] MetaData Writer class was added for writing the simulation parameters to a YAML file.
- [ALL] Adds pipeline for performance tests with weekly schedule
- [DOCS] Adds Coding Guidelines to the documentation. The guidelines can be found [here](https://irmb.gitlab-pages.rz.tu-bs.de/VirtualFluids/coding-guidelines.html).
- [GPU] Adds laminar pipe flow app with regression test and performance test
- [GPU] Adds helper Class for creating subdomains on multiple GPUs
- [ALL] new docker image for CUDA 12.4.1 (Docker Image Version 1.3) 

### Fixes
- [GPU] Pressure BoundaryCondition was not working correctly for all directions. This is now fixed.
- [GPU] Density calculation in PressureNonEquilibriumCompressible was incorrect and not conditioned. This is now fixed. 
- [GPU] Fix weekly schedule for multi-GPU testcases on high-performance compute cluster
- [GPU] Fix NUPS calculation at the end of the simulation

### Refactoring
- [CPU] move MetisPartitioner to parallel
- [GPU] Probes and Interactors

## [0.1.1](https://git.rz.tu-bs.de/irmb/VirtualFluids/-/milestones/2) - 2024-01-12

### Added
- [DOCS] We reworked the doxygen documentation. It is now more structured and easier to read. You can find it [here](https://irmb.gitlab-pages.rz.tu-bs.de/VirtualFluids/). The documentation theme is based on [awesome doxygen](https://jothepro.github.io/doxygen-awesome-css/). Additionally, it contains the [wiki pages](https://irmb.gitlab-pages.rz.tu-bs.de/VirtualFluids/documentation.html) and a reworked [module](https://irmb.gitlab-pages.rz.tu-bs.de/VirtualFluids/topics.html) overview of the code. 
- [TESTS] All tests are now located in the `tests` folder. The `tests` folder is located in the root directory of the project. It contains unit-tests, regression-tests and numerical tests running on the gpu.
- [DOCS] the project is now [reuse](https://reuse.software/) compliant. All files are licensed and the licenses can be found in the `LICENSES` folder.
- [GPU] the PrecollisionIntercator class does not have to be freed manually anymore. The destructor of the classes will do this for you.

### Fixes
- [CPU] we removed a virtual function call in the cpu collision kernel.
- [GPU] Rework app SphereInChannel, as there was an error in calculating the viscosity.
- [BUILD] code coverage is now working again.

## [0.1.0](https://git.rz.tu-bs.de/irmb/VirtualFluids/-/milestones/1) - 2023-12-07

### Added
- Initial release

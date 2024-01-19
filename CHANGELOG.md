<!-- SPDX-License-Identifier: CC-BY-4.0 -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->
# Changelog

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

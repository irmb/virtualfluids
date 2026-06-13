<!-- SPDX-License-Identifier: CC-BY-4.0 -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->
# Changelog
## [0.3.0](https://git.rz.tu-bs.de/irmb/VirtualFluids/-/milestones/5) - 2026-06-08
### Added
- [ALL] Added fast winding algorithm for discretization of STLs for the numerical grids.  
- [GPU] Added Hub and Tower drag to ActuatorFarm class.
- [ALL] Update devcontainer to ubuntu 24.04.
- [ALL] Added grid refinement for advection diffusion.
- [GPU] Add option for spatially heterogeneous roughness wall model.
- [GPU] Add precursor and outflow BCs for advection diffusion.
- [GPU] Add vertical axis wind turbine app.
- [GPU] Add actuator disc app.
- [GPU] Add homogeneous forest in atmospheric boundary layer.

### Fixes
- [GPU] Bug Fix non reflected pressure correction outflow.
- [GPU] Fix periodic stopper neighbors for GPU.
- [GPU] Fix planar average probe, when refinement covers whole domain.
- [GPU] Fix device memory allocator (AdvectionDiffusionFluxBC).
- [GPU] Fix rounding errors in precursor.
- [GPU] Fix rescaling of buoyancy factor in refined areas.
- [GPU] Several fixes in python bindings.
- [GPU] Fix rounding bug in communication indices.
- [GPU] Fix coefficient in moeng turbulent diffusivity model.
- [ALL] Fix bug in combination with newer CMake Versions.

### Refactoring
- [ALL] Compile all VirtualFluids libraries into one python binding.
- [GPU] Move gpu library into namespace vf::gpu.  

## [0.2.1]() - 2025-12-11
### Added
- [CPU] Added comprehensive unit tests to verify the VTI reading capabilities.
- [CPU] Introduced a new NoSlipInterpolatedRelaxed boundary condition, providing more flexibility in simulations.
- [CPU] Added extensive unit tests for the setQuadricLimiter method to ensure correct behavior with valid and invalid inputs.
- [CPU] Added new CPU application: Flow Around Cylinder Benchmark.
- [GPU] Performance improvement by using atimic exchange in runcollision to read and write the forces in one step.
- [GPU] Add a regression test (Gaussian Hill) for a multi gpu domain with periodic BCs.

### Fixes
- [GPU] Fix the path of the file name for the timeseries output file in the probe.
- [GPU] Fix bug that only appeared when using more than 4 gpus, and some cleanup in the index rearrangement.
- [GPU] Fix bug that appeared when using 2x2 grid of 4 gpus and periodic BCs in both directions.
- [CPU] Fix config files for laminar plane flow app.

### Refactoring
- [CPU] VTI File support in GbVoxelMatrix3D
- [CPU] Improved the CalculateForcesSimulationObserver by refactoring it to use a dedicated ForceCalculator and removing redundant swap calls.
- [CPU] Applied conditioning to interpolation coefficients for improved numerical stability.

## [0.2.0]() - 2025-10-23
### Added
- [ALL] More refined CMake presets including build presets
- [GPU] Checks and recommendations for viscosity and velocity for K17 kernel
- [GPU] Ability to sample concentration with samplers
- [GPU] Adds multi-GPU communication hiding for advection-diffusion
- [GPU] Adds turbulent diffusivity models for advection-diffusion

- [GPU] Adds interactors for buoyancy and coriolis force computation
- [GPU] Adds example apps for advection-diffusion: GaussianHillOfConcentration and HeatedCube
- [GPU] Adds new "surface layer" boundary condition to apply Monin-Obukhov similarity theory for thermally stratified atmospheric boundary layers

### Fixes
- [GPU] fix in multiple Grid Builder Facade for Multi GPU and periodic boundary conditions
- [GPU] fix bug in probe for cases with refinement
- [GPU] fix computation of temporal averages in WallModel and PlanarAverage probes

### Refactoring
- [GPU] Make boundary conditions for advection diffusion solver more flexible
- [GPU] Make use of template programming and constant expressions to simplify many boundary conditions
- [GPU] Body force field contains level wise forces in lattice units of related level instead of coarse level

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

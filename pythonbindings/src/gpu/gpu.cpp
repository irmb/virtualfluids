//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file gpu.cpp
//! \ingroup gpu
//! \author Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include "submodules/pre_collision_interactor.cpp"
#include "submodules/simulation.cpp"
#include "submodules/parameter.cpp"
#include "submodules/boundary_conditions.cpp"
#include "submodules/communicator.cpp"
#include "submodules/cuda_memory_manager.cpp"
#include "submodules/probes.cpp"
#include "submodules/precursor_writer.cpp"
#include "submodules/grid_provider.cpp"
#include "submodules/grid_generator.cpp"
#include "submodules/turbulence_models.cpp"
#include "submodules/transient_bc_setter.cpp"
#include "submodules/actuator_farm.cpp"
#include "submodules/grid_scaling_factory.cpp"

namespace gpu
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module gpuModule = parentModule.def_submodule("gpu");
        simulation::makeModule(gpuModule);
        parameter::makeModule(gpuModule);
        pre_collision_interactor::makeModule(gpuModule);
        actuator_farm::makeModule(gpuModule);
        boundary_conditions::makeModule(gpuModule);
        transient_bc_setter::makeModule(gpuModule);
        communicator::makeModule(gpuModule); 
        cuda_memory_manager::makeModule(gpuModule);
        probes::makeModule(gpuModule);
        precursor_writer::makeModule(gpuModule);
        grid_generator::makeModule(gpuModule);
        grid_provider::makeModule(gpuModule);
        turbulence_model::makeModule(gpuModule);
        grid_scaling_factory::makeModule(gpuModule);
        return gpuModule;
    }
}
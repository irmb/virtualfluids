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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \author Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include "submodules/pre_collision_interactor.cpp"
#include "submodules/simulation.cpp"
#include "submodules/parameter.cpp"
#include "submodules/boundary_conditions.cpp"
#include "submodules/cuda_memory_manager.cpp"
#include "submodules/probes.cpp"
#include "submodules/precursor_writer.cpp"
#include "submodules/grid_provider.cpp"
#include "submodules/grid_generator.cpp"
#include "submodules/turbulence_models.cpp"
#include "submodules/transient_bc_setter.cpp"
#include "submodules/actuator_farm.cpp"
#include "submodules/grid_scaling_factory.cpp"
#include "submodules/kernel.cpp"

namespace gpu_bindings
{
PYBIND11_MODULE(gpu, m)
{
    simulation::makeModule(m);
    parameter::makeModule(m);
    pre_collision_interactor::makeModule(m);
    actuator_farm::makeModule(m);
    boundary_conditions::makeModule(m);
    transient_bc_setter::makeModule(m);
    cuda_memory_manager::makeModule(m);
    probes::makeModule(m);
    precursor_writer::makeModule(m);
    grid_generator::makeModule(m);
    grid_provider::makeModule(m);
    turbulence_model::makeModule(m);
    grid_scaling_factory::makeModule(m);
    kernel::makeModule(m);
}
} // namespace gpu_bindings

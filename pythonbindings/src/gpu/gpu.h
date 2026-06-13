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
#include "submodules/pre_collision_interactor.h"
#include "submodules/simulation.h"
#include "submodules/parameter.h"
#include "submodules/boundary_conditions.h"
#include "submodules/cuda_memory_manager.h"
#include "submodules/probes.h"
#include "submodules/precursor_writer.h"
#include "submodules/grid_provider.h"
#include "submodules/grid_generator.h"
#include "submodules/turbulence_models.h"
#include "submodules/transient_bc_setter.h"
#include "submodules/actuator_farm.h"
#include "submodules/grid_scaling_factory.h"
#include "submodules/kernel.h"
#include "submodules/sampler.h"
#include "submodules/forest.h"
#include "submodules/coriolis_force.h"
#include "submodules/damping_layer.h"
#include "submodules/buoyancy_provider.h"

namespace gpu_bindings
{
        
namespace py = pybind11;

inline void makeModule(py::module_& pyfluids)
{
    auto gpu = pyfluids.def_submodule("gpu");
    simulation::makeModule(gpu);
    parameter::makeModule(gpu);
    pre_collision_interactor::makeModule(gpu);
    actuator_farm::makeModule(gpu);
    forest::makeModule(gpu);
    coriolis_force::makeModule(gpu);
    boundary_conditions::makeModule(gpu);
    transient_bc_setter::makeModule(gpu);
    cuda_memory_manager::makeModule(gpu);
    sampler::makeModule(gpu);
    probes::makeModule(gpu);
    precursor_writer::makeModule(gpu);
    grid_generator::makeModule(gpu);
    grid_provider::makeModule(gpu);
    turbulence_model::makeModule(gpu);
    grid_scaling_factory::makeModule(gpu);
    kernel::makeModule(gpu);
    damping_layer::makeModule(gpu);
    buoyancy_provider::makeModule(gpu);
}
} // namespace gpu_bindings

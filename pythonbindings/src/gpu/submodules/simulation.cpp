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
#include <gpu/core/Calculation/Simulation.h>
#include <gpu/core/Kernel/KernelFactory/KernelFactory.h>
#include <gpu/core/PreProcessor/PreProcessorFactory/PreProcessorFactory.h>
#include <gpu/core/DataStructureInitializer/GridProvider.h>
#include <gpu/core/Parameter/Parameter.h>
#include <gpu/core/Cuda/CudaMemoryManager.h>
#include <gpu/core/Output/DataWriter.h>
#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "parallel/Communicator.h"

namespace simulation
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        // missing setFactories and setDataWriter, not possible to wrap these functions as long as they take unique ptr arguments
        py::class_<Simulation>(parentModule, "Simulation")
        .def(py::init<std::shared_ptr<Parameter>,
                      std::shared_ptr<GridBuilder>,
                      BoundaryConditionFactory*>(),
             py::arg("parameter"), py::arg("grid_builder"), py::arg("bc_factory"))
        .def(py::init<std::shared_ptr<Parameter>,
                      std::shared_ptr<GridBuilder>,
                      BoundaryConditionFactory*,
                      GridScalingFactory*>(),
             py::arg("parameter"), py::arg("grid_builder"), py::arg("bc_factory"), py::arg("grid_scaling_factory"))

        .def(py::init<std::shared_ptr<Parameter>,
                      std::shared_ptr<CudaMemoryManager>,
                      std::shared_ptr<GridBuilder>,
                      BoundaryConditionFactory*>(),
             py::arg("parameter"), py::arg("memory_manager"), py::arg("grid_builder"), py::arg("bc_factory"))
        .def(py::init<std::shared_ptr<Parameter>,
                      std::shared_ptr<CudaMemoryManager>,
                      std::shared_ptr<GridBuilder>,
                      BoundaryConditionFactory*,
                      GridScalingFactory*>(),
             py::arg("parameter"), py::arg("memory_manager"), py::arg("grid_builder"), py::arg("bc_factory"),
             py::arg("scaling_factory"))

        .def(py::init<std::shared_ptr<Parameter>,
                      std::shared_ptr<GridBuilder>,
                      BoundaryConditionFactory*,
                      std::shared_ptr<TurbulenceModelFactory>>(),
             py::arg("parameter"), py::arg("grid_builder"), py::arg("bc_factory"), py::arg("tm_factory"))
        .def(py::init<std::shared_ptr<Parameter>,
                      std::shared_ptr<GridBuilder>,
                      BoundaryConditionFactory*,
                      std::shared_ptr<TurbulenceModelFactory>,
                      GridScalingFactory*>(),
             py::arg("parameter"), py::arg("grid_builder"), py::arg("bc_factory"), py::arg("tm_factory"),
             py::arg("grid_scaling_factory"))

        
        .def(py::init<std::shared_ptr<Parameter>,
                      std::shared_ptr<CudaMemoryManager>,
                      std::shared_ptr<GridBuilder>,
                      BoundaryConditionFactory*,
                      std::shared_ptr<TurbulenceModelFactory>>(),
             py::arg("parameter"), py::arg("memory_manager"), py::arg("grid_builder"), py::arg("bc_factory"), py::arg("tm_factory"))
        .def(py::init<std::shared_ptr<Parameter>,
                      std::shared_ptr<CudaMemoryManager>,
                      std::shared_ptr<GridBuilder>,
                      BoundaryConditionFactory*,
                      std::shared_ptr<TurbulenceModelFactory>,
                      GridScalingFactory*>(),
             py::arg("parameter"), py::arg("memory_manager"), py::arg("grid_builder"), py::arg("bc_factory"), py::arg("tm_factory"),
             py::arg("grid_scaling_factory"))
        .def("run", &Simulation::run)
        .def("init_timers", &Simulation::initTimers)
        .def("calculate_timestep", &Simulation::calculateTimestep, py::arg("time_step"))
        .def("addKineticEnergyAnalyzer", &Simulation::addKineticEnergyAnalyzer, py::arg("t_analyse"))
        .def("addEnstrophyAnalyzer", &Simulation::addEnstrophyAnalyzer, py::arg("t_analyse"));
    }
}
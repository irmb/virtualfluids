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
//! \file simulation.cpp
//! \ingroup submodules
//! \author Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include <gpu/core/Calculation/Simulation.h>
#include <gpu/core/Kernel/KernelFactory/KernelFactory.h>
#include <gpu/core/PreProcessor/PreProcessorFactory/PreProcessorFactory.h>
#include <gpu/core/DataStructureInitializer/GridProvider.h>
#include <gpu/core/Parameter/Parameter.h>
#include <gpu/core/Cuda/CudaMemoryManager.h>
#include <gpu/core/DataStructureInitializer/GridProvider.h>
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
        .def(py::init<  std::shared_ptr<Parameter>,
                        std::shared_ptr<CudaMemoryManager>,
                        vf::parallel::Communicator &,
                        GridProvider &,
                        BoundaryConditionFactory*,
                        GridScalingFactory*>(), 
                        py::arg("parameter"),
                        py::arg("memoryManager"),
                        py::arg("communicator"),
                        py::arg("gridProvider"),
                        py::arg("bcFactory"),
                        py::arg("gridScalingFactory"))
        .def(py::init<  std::shared_ptr<Parameter>,
                        std::shared_ptr<CudaMemoryManager>,
                        vf::parallel::Communicator &,
                        GridProvider &,
                        BoundaryConditionFactory*>(), 
                        py::arg("parameter"),
                        py::arg("memoryManager"),
                        py::arg("communicator"),
                        py::arg("gridProvider"),
                        py::arg("bcFactory"))
        .def(py::init<  std::shared_ptr<Parameter>,
                        std::shared_ptr<CudaMemoryManager>,
                        vf::parallel::Communicator &,
                        GridProvider &,
                        BoundaryConditionFactory*,
                        std::shared_ptr<TurbulenceModelFactory>,
                        GridScalingFactory*>(), 
                        py::arg("parameter"),
                        py::arg("memoryManager"),
                        py::arg("communicator"),
                        py::arg("gridProvider"),
                        py::arg("bcFactory"),
                        py::arg("tmFactory"),
                        py::arg("gridScalingFactory"))
        .def("run", &Simulation::run)
        .def("addKineticEnergyAnalyzer", &Simulation::addKineticEnergyAnalyzer, py::arg("t_analyse"))
        .def("addEnstrophyAnalyzer", &Simulation::addEnstrophyAnalyzer, py::arg("t_analyse"));
    }
}
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
//! \file simulationconfig.cpp
//! \ingroup submodules
//! \author Sven Marcus, Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include <simulationconfig/Simulation.h>

namespace simulation
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<CPUSimulation, std::shared_ptr<CPUSimulation>>(parentModule, "Simulation")
                .def(py::init())
                .def("set_writer", &CPUSimulation::setWriterConfiguration)
                .def("set_grid_parameters", &CPUSimulation::setGridParameters)
                .def("set_physical_parameters", &CPUSimulation::setPhysicalParameters)
                .def("set_runtime_parameters", &CPUSimulation::setRuntimeParameters)
                .def("set_kernel_config", &CPUSimulation::setKernelConfiguration)
                .def("add_object", &CPUSimulation::addObject)
                .def("add_bc_adapter", &CPUSimulation::addBCAdapter)
                .def("run_simulation", &CPUSimulation::run);
    }

}
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
//! \author Mohammad Mehdi Mohammadi
//=======================================================================================

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <gpu/core/PreCollisionInteractor/Forest.h>
#include <gpu/core/PreCollisionInteractor/PreCollisionInteractor.h>

namespace gpu_bindings::forest
{
    namespace py = pybind11;
    using namespace vf::gpu;

    void makeModule(py::module_ &parentModule)
    {
        py::enum_<PlantAreaDensityFunction>(parentModule, "PlantAreaDensityFunction")
            .value("Flat", PlantAreaDensityFunction::Flat)
            .value("BottomHeavy", PlantAreaDensityFunction::BottomHeavy)
            .value("TopHeavy", PlantAreaDensityFunction::TopHeavy);

        py::class_<Forest, PreCollisionInteractor, std::shared_ptr<Forest>>(parentModule, "Forest")
            .def(py::init<
                SPtr<Parameter>,
                SPtr<CudaMemoryManager>,
                const real,
                const real,
                const real,
                const int,
                PlantAreaDensityFunction
            >(),
                py::arg("para"),
                py::arg("cuda_memory_manager"),
                py::arg("height"),
                py::arg("dragCoeff"),
                py::arg("plantAreaIndex"),
                py::arg("level"),
                py::arg("plantAreaDensityFunction")
            )
            .def_property_readonly("number_of_indices", &Forest::getNumberOfIndices);
    }
}
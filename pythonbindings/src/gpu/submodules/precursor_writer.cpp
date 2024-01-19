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
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <gpu/core/PreCollisionInteractor/PreCollisionInteractor.h>
#include <gpu/core/PreCollisionInteractor/PrecursorWriter.h>

namespace precursor_writer
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::enum_<OutputVariable>(parentModule, "OutputVariable")
        .value("Velocities", OutputVariable::Velocities)
        .value("Distributions", OutputVariable::Distributions);

        py::class_<PrecursorWriter, PreCollisionInteractor, std::shared_ptr<PrecursorWriter>>(parentModule, "PrecursorWriter")
        .def(py::init < std::string,
                        std::string,
                        real,
                        real, real,
                        real, real,
                        uint, uint, 
                        OutputVariable, 
                        uint>(),
                        py::arg("filename"),
                        py::arg("output_path"), 
                        py::arg("x_pos"),
                        py::arg("y_min"), py::arg("y_max"),
                        py::arg("z_min"), py::arg("z_max"),
                        py::arg("t_start_out"), py::arg("t_save"), 
                        py::arg("output_variable"), 
                        py::arg("max_timesteps_per_file"));
    }
}
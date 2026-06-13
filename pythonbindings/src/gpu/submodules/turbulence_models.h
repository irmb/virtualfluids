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
#include "pybind11/pybind11.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"
#include <lbm/collision/TurbulentViscosity.h>
#include <lbm/advectionDiffusion/TurbulentDiffusivity.h>
#include <string>

namespace gpu_bindings::turbulence_model
{
    namespace py = pybind11;
    using namespace vf::lbm;
    using namespace vf::gpu;

    inline void makeModule(py::module_ &parentModule)
    {
        py::class_<TurbulenceModelFactory, std::shared_ptr<TurbulenceModelFactory>>(parentModule, "TurbulenceModelFactory")
        .def(py::init< std::shared_ptr<Parameter>>(), py::arg("para"))
        .def("set_turbulence_model", py::overload_cast<std::string>(&TurbulenceModelFactory::setTurbulenceModel), py::arg("turbulence_model"))
        .def("set_turbulence_model", py::overload_cast<TurbulenceModel>(&TurbulenceModelFactory::setTurbulenceModel), py::arg("turbulence_model"))
        .def("set_turbulence_model_advection_diffusion", py::overload_cast<std::string>(&TurbulenceModelFactory::setAdvectionDiffusionTurbulenceModel), py::arg("turbulence_model"))
        .def("set_turbulence_model_advection_diffusion", py::overload_cast<advection_diffusion::TurbulenceModel>(&TurbulenceModelFactory::setAdvectionDiffusionTurbulenceModel), py::arg("turbulence_model"))
        .def("set_model_constant", &TurbulenceModelFactory::setModelConstant, py::arg("model_constant"))
        .def("read_config_file", &TurbulenceModelFactory::readConfigFile, py::arg("config_data"));
    }
}
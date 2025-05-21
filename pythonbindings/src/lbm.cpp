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
#include <lbm/collision/TurbulentViscosity.h>
#include <lbm/advectionDiffusion/TurbulentDiffusivity.h>

namespace lbm_bindings
{
    namespace py = pybind11;

    PYBIND11_MODULE(lbm, m)
    {

        py::enum_<vf::lbm::TurbulenceModel>(m, "TurbulenceModel")
        .value("None", vf::lbm::TurbulenceModel::None)
        .value("Smagorinsky", vf::lbm::TurbulenceModel::Smagorinsky)
        .value("QR", vf::lbm::TurbulenceModel::QR)
        .value("AMD", vf::lbm::TurbulenceModel::AMD)
        .value("AMDStratified", vf::lbm::TurbulenceModel::AMDStratified);

        py::enum_<vf::lbm::advection_diffusion::TurbulenceModel>(m, "TurbulenceModelAdvectionDiffusion")
        .value("None", vf::lbm::advection_diffusion::TurbulenceModel::None)
        .value("Default", vf::lbm::advection_diffusion::TurbulenceModel::Default)
        .value("Moeng", vf::lbm::advection_diffusion::TurbulenceModel::Moeng)
        .value("AMDStratified", vf::lbm::advection_diffusion::TurbulenceModel::AMDStratified);
    }
}
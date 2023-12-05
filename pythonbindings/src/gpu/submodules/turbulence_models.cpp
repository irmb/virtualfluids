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
//! \file turbulence_models.cpp
//! \ingroup submodules
//! \author Henry Korb
//=======================================================================================
#include "pybind11/pybind11.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"
#include "gpu/core/Calculation/Calculation.h"

namespace turbulence_model
{
    namespace py = pybind11;
    using namespace vf::lbm;

    void makeModule(py::module_ &parentModule)
    {
        py::enum_<TurbulenceModel>(parentModule, "TurbulenceModel")
        .value("Smagorinsky", TurbulenceModel::Smagorinsky)
        .value("AMD", TurbulenceModel::AMD)
        .value("QR", TurbulenceModel::QR)
        .value("None", TurbulenceModel::None);

        py::class_<TurbulenceModelFactory, std::shared_ptr<TurbulenceModelFactory>>(parentModule, "TurbulenceModelFactory")
        .def(py::init< std::shared_ptr<Parameter>>(), py::arg("para"))
        .def("set_turbulence_model", &TurbulenceModelFactory::setTurbulenceModel, py::arg("turbulence_model"))
        .def("set_model_constant", &TurbulenceModelFactory::setModelConstant, py::arg("model_constant"))
        .def("read_config_file", &TurbulenceModelFactory::readConfigFile, py::arg("config_data"));

    }
}
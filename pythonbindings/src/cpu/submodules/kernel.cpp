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
//! \author Sven Marcus, Henry Korb
//=======================================================================================
#include <memory>
#include <pybind11/pybind11.h>
#include <simulationconfig/KernelFactory.h>
#include <simulationconfig/KernelConfigStructs.h>

namespace kernel
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::module kernelModule = parentModule.def_submodule("kernel");

        py::enum_<KernelFactory::KernelType>(kernelModule, "KernelType")
                .value("BGK", KernelFactory::BGK)
                .value("CompressibleCumulantFourthOrderViscosity",
                       KernelFactory::COMPRESSIBLE_CUMULANT_4TH_ORDER_VISCOSITY);

        py::class_<LBMKernelConfiguration, std::shared_ptr<LBMKernelConfiguration>>(kernelModule, "LBMKernel")
                .def(py::init<KernelFactory::KernelType>())
                .def_readwrite("type", &LBMKernelConfiguration::kernelType)
                .def_readwrite("use_forcing", &LBMKernelConfiguration::useForcing)
                .def_readwrite("forcing_in_x1", &LBMKernelConfiguration::forcingX1)
                .def_readwrite("forcing_in_x2", &LBMKernelConfiguration::forcingX2)
                .def_readwrite("forcing_in_x3", &LBMKernelConfiguration::forcingX3)
                .def("set_forcing", [](LBMKernelConfiguration &kernelConfig, real x1, real x2, real x3)
                {
                    kernelConfig.forcingX1 = x1;
                    kernelConfig.forcingX2 = x2;
                    kernelConfig.forcingX3 = x3;
                })
                .def("__repr__", [](LBMKernelConfiguration &kernelConfig)
                {
                    std::ostringstream stream;
                    stream << "<" << kernelConfig.kernelType << std::endl
                           << "Use forcing: " << kernelConfig.useForcing << std::endl
                           << "Forcing in x1: " << kernelConfig.forcingX1 << std::endl
                           << "Forcing in x2: " << kernelConfig.forcingX2 << std::endl
                           << "Forcing in x3: " << kernelConfig.forcingX3 << ">" << std::endl;

                    return stream.str();
                });
    }

}
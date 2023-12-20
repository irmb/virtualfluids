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
#include <gpu/core/Kernel/KernelTypes.h>

namespace kernel
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        auto kernel_module = parentModule.def_submodule("Kernel", "Kernel types");
        auto compressible = kernel_module.def_submodule("compressible", "Compressible Kernel types");
        auto incompressible = kernel_module.def_submodule("incompressible", "Incompressible Kernel types");

        compressible.attr("BGK") = vf::collisionKernel::compressible::BGK;
        compressible.attr("BGKPlus") = vf::collisionKernel::compressible::BGKPlus;
        compressible.attr("K17CompressibleNavierStokes") = vf::collisionKernel::compressible::K17CompressibleNavierStokes;
        compressible.attr("K15CompressibleNavierStokes") = vf::collisionKernel::compressible::K15CompressibleNavierStokes;

        incompressible.attr("BGK") = vf::collisionKernel::incompressible::BGK;
        incompressible.attr("BGKPlus") = vf::collisionKernel::incompressible::BGKPlus;
        incompressible.attr("CumulantK15") = vf::collisionKernel::incompressible::CumulantK15;
    }
}

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
#include <LBM/LBMKernel.h>
#include <LBM/K17CompressibleNavierStokes.h>
#include <LBM/B92IncompressibleNavierStokes.h>
#include <simulationconfig/D3Q27LBMSystem.h>
#include "simulationconfig/KernelFactory.h"

std::shared_ptr<LBMKernel> KernelFactory::makeKernel(KernelType kernelType)
{
    switch (kernelType) {
        case BGK:
            return std::shared_ptr<LBMKernel>(new B92IncompressibleNavierStokes());
        case COMPRESSIBLE_CUMULANT_4TH_ORDER_VISCOSITY:
            return std::shared_ptr<LBMKernel>(new K17CompressibleNavierStokes());
        default:
            throw std::logic_error("No such kernel type");
    }
}

std::shared_ptr<AbstractLBMSystem> KernelFactory::makeLBMSystem(KernelType type)
{
    return std::shared_ptr<AbstractLBMSystem>(new D3Q27LBMSystem());
}

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
//! \file TurbulenceModelManager.cpp
//! \ingroup KernelManager
//! \author Henrik Asmuth
//=======================================================================================
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <iostream>
#include <stdexcept>
#include <string>

#include "TurbulenceModelManager.h"
#include "TurbulenceModels/TurbulenceModelFactory.h"
#include "Parameter/Parameter.h"

TurbulenceModelManager::TurbulenceModelManager(SPtr<Parameter> parameter, TurbulenceModelFactory *turbulenceModelFactory) : para(parameter)
{
    this->turbulenceModelKernel = turbulenceModelFactory->getTurbulenceModelKernel();

    // checkBoundaryCondition(this->velocityBoundaryConditionPost, this->para->getParD(0)->velocityBC,
    //                        "velocityBoundaryConditionPost");
}

void TurbulenceModelManager::runTurbulenceModelKernel(const int level) const
{
    this->turbulenceModelKernel(para.get(), level);
}

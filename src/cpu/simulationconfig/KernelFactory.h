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
#ifndef VIRTUALFLUIDSPYTHONBINDINGS_KERNELFACTORY_H
#define VIRTUALFLUIDSPYTHONBINDINGS_KERNELFACTORY_H

#include <LBM/LBMKernel.h>
#include "AbstractLBMSystem.h"


class KernelFactory {
public:
    enum KernelType {
        BGK,
        COMPRESSIBLE_CUMULANT_4TH_ORDER_VISCOSITY /*,
        COMPRESSIBLE_CUMULANT,
        CUMULANT_K17,
        INCOMPRESSIBLE_CUMULANT,
        INCOMPRESSIBLE_CUMULANT_WITH_SPONGE_LAYER,
        INIT_DENSITIY,
        ET_D3Q27_BGK
        */
    };

    KernelFactory() = default;
    virtual ~KernelFactory() = default;

    std::shared_ptr<LBMKernel> makeKernel(KernelType kernelType);

    std::shared_ptr<AbstractLBMSystem> makeLBMSystem(KernelType type);
};


#endif //VIRTUALFLUIDSPYTHONBINDINGS_KERNELFACTORY_H

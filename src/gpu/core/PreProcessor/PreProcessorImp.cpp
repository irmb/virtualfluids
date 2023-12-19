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
//! \addtogroup gpu_PreProcessor PreProcessor
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "PreProcessorImp.h"

#include "PreProcessorStrategy/PreProcessorStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<PreProcessorImp> PreProcessorImp::getNewInstance()
{
    return std::shared_ptr<PreProcessorImp>(new PreProcessorImp());
}

void PreProcessorImp::addStrategy(std::shared_ptr<PreProcessorStrategy> strategy)
{
    strategies.push_back(strategy);
}

void PreProcessorImp::init(std::shared_ptr<Parameter> para, int level)
{
    para->getParD(level)->isEvenTimestep = false;
    for (std::size_t i = 0; i < strategies.size(); i++)
        strategies.at(i)->init(level);

    para->getParD(level)->isEvenTimestep = true;
    for (std::size_t i = 0; i < strategies.size(); i++)
        strategies.at(i)->init(level);
}

bool PreProcessorImp::checkParameter()
{
    for (std::size_t i = 0; i < strategies.size(); i++) {
        if (!strategies.at(i)->checkParameter())
            return false;
    }
    return true;
}

PreProcessorImp::PreProcessorImp()
{
    strategies.resize(0);
}

//! \}

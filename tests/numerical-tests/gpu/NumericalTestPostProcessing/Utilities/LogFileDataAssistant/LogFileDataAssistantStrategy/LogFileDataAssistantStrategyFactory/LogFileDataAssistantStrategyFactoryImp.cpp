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
//! \addtogroup NumericalTestPostProcessing
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "LogFileDataAssistantStrategyFactoryImp.h"

#include "Simulation/ShearWave/LogFileDataAssistantStrategy/ShearWaveLogFileDataAssistantStrategy.h"
#include "Simulation/TaylorGreenVortexUx/LogFileDataAssistantStrategy/TaylorGreenVortexUxLogFileDataAssistantStrategy.h"
#include "Simulation/TaylorGreenVortexUz/LogFileDataAssistantStrategy/TaylorGreenVortexUzLogFileDataAssistantStrategy.h"

std::shared_ptr<LogFileDataAssistantStrategyFactory> LogFileDataAssistantStrategyFactoryImp::getNewInstance()
{
    return std::shared_ptr<LogFileDataAssistantStrategyFactory>(new LogFileDataAssistantStrategyFactoryImp());
}

std::shared_ptr<LogFileDataAssistantStrategy> LogFileDataAssistantStrategyFactoryImp::makeLogFileDataAssistantStrategy(BasicSimulation sim)
{
    std::shared_ptr<LogFileDataAssistantStrategy> assistentStrategy;
    switch (sim)
    {
    case ShearWave:
        assistentStrategy = ShearWaveLogFileDataAssistantStrategy::getNewInstance();
        break;
    case TaylorGreenVortexUx:
        assistentStrategy = TaylorGreenVortexUxLogFileDataAssistantStrategy::getNewInstance();
        break;
    case TaylorGreenVortexUz:
        assistentStrategy = TaylorGreenVortexUzLogFileDataAssistantStrategy::getNewInstance();
        break;
    default:
        break;
    }
    return assistentStrategy;
}

LogFileDataAssistantStrategyFactoryImp::LogFileDataAssistantStrategyFactoryImp()
{
}

//! \}

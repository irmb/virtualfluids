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
#include "TaylorGreenVortexUxLogFileDataAssistantStrategy.h"

#include "Utilities/LogFileData/LogFileData.h"
#include "Simulation/TaylorGreenVortexUx/LogFileData/TaylorGreenVortexUxLogFileData.h"

std::shared_ptr<LogFileDataAssistantStrategy> TaylorGreenVortexUxLogFileDataAssistantStrategy::getNewInstance()
{
    return std::shared_ptr<LogFileDataAssistantStrategy>(new TaylorGreenVortexUxLogFileDataAssistantStrategy());
}

std::string TaylorGreenVortexUxLogFileDataAssistantStrategy::getSimulationName()
{
    return simName;
}

bool TaylorGreenVortexUxLogFileDataAssistantStrategy::checkSimulationParameter(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
    if (!equalDouble(logFileData1->getTaylorGreenVortexUxLogFileData()->getUx().at(0), logFileData2->getTaylorGreenVortexUxLogFileData()->getUx().at(0)))
        return false;
    if (!equalDouble(logFileData1->getTaylorGreenVortexUxLogFileData()->getAmplitude().at(0), logFileData2->getTaylorGreenVortexUxLogFileData()->getAmplitude().at(0)))
        return false;
    if (logFileData1->getTaylorGreenVortexUxLogFileData()->getL0().at(0) != logFileData2->getTaylorGreenVortexUxLogFileData()->getL0().at(0))
        return false;

    return true;
}

TaylorGreenVortexUxLogFileDataAssistantStrategy::TaylorGreenVortexUxLogFileDataAssistantStrategy()
{
    this->simName = "TaylorGreenVortexUx";
}

//! \}

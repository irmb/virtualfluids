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
//! \addtogroup NumericalTests
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "TimeImp.h"

#include <sstream>

std::shared_ptr<TimeImp> TimeImp::getNewInstance()
{
    return std::shared_ptr<TimeImp>(new TimeImp());
}

TimeImp::TimeImp()
{

}

void TimeImp::setSimulationStartTime()
{
    simulationStartTime = time(NULL);
}

void TimeImp::setSimulationEndTime()
{
    simulationEndTime = time(NULL);
}

void TimeImp::setTestStartTime()
{
    testStartTime = clock();
}

void TimeImp::setTestEndTime()
{
    testEndTime = clock();
}

void TimeImp::setAnalyticalResultWriteStartTime()
{
    analyticalResultWriteStartTime = time(NULL);
}

void TimeImp::setAnalyticalResultWriteEndTime()
{
    analyticalResultWriteEndTime = time(NULL);
}

void TimeImp::setResultCheckStartTime()
{
    resultCheckStartTime = clock();
}

void TimeImp::setResultCheckEndTime()
{
    resultCheckEndTime = clock();
}

std::string TimeImp::getSimulationTime()
{
    std::ostringstream oss;
    oss << calcSimulationTime() << "sec";
    return oss.str();
}

std::string TimeImp::getResultCheckTime()
{
    std::ostringstream oss;
    oss << calcResultCheckTime() << "sec";
    return oss.str();
}

std::string TimeImp::getTestTime()
{
    std::ostringstream oss;
    oss << calcTestTime() << "sec";
    return oss.str();
}

std::string TimeImp::getAnalyticalResultWriteTime()
{
    std::ostringstream oss;
    oss << calcAnalyticalResultWriteTime() << "sec";
    return oss.str();
}

double TimeImp::calcSimulationTime()
{
    return difftime(simulationEndTime, simulationStartTime);
}

float TimeImp::calcResultCheckTime()
{
    float timeInMiliSec = ((float)(resultCheckEndTime - resultCheckStartTime) / CLOCKS_PER_SEC);
    return timeInMiliSec;
}

float TimeImp::calcTestTime()
{
    float timeInMiliSec = ((float)(testEndTime - testStartTime) / CLOCKS_PER_SEC);
    return timeInMiliSec;
}

double TimeImp::calcAnalyticalResultWriteTime()
{
    return difftime(analyticalResultWriteEndTime, analyticalResultWriteStartTime);
}

//! \}

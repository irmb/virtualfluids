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
#include "SimulationInfoImp.h"

#include "Utilities/Time/TimeInfo.h"

#include <iomanip>
#include <sstream>

void SimulationInfoImp::setTimeInfo(std::shared_ptr<TimeInfo> timeInfo)
{
    this->timeInfo = timeInfo;
}

std::string SimulationInfoImp::getKernelName()
{
    return kernelName;
}

double SimulationInfoImp::getViscosity()
{
    return viscosity;
}

std::string SimulationInfoImp::getSimulationName()
{
    return simulationName;
}

std::string SimulationInfoImp::getSimulationParameterString()
{
    return simulationParameterString;
}

int SimulationInfoImp::getLx()
{
    return lx;
}

int SimulationInfoImp::getNumberOfSimulations()
{
    return numberOfSimulations;
}

int SimulationInfoImp::getSimulationID()
{
    return simID;
}

std::string SimulationInfoImp::getRunTimeOutput()
{
    std::ostringstream oss;
    oss << "SimulationTime_" << lx << "=" << timeInfo->getSimulationTime() << std::endl;
    oss << "ResultsCheckTime_" << lx << "=" << timeInfo->getResultCheckTime() << std::endl;
    oss << "TestTime_" << lx << "=" << timeInfo->getTestTime() << std::endl;
    oss << "AnalyticalVTKFileWritingTime_" << lx << "=" << timeInfo->getAnalyticalResultWriteTime() << std::endl;
    return oss.str();
}

std::vector<std::string> SimulationInfoImp::getDataToCalcTests()
{
    return dataToCalcTests;
}

SimulationInfoImp::SimulationInfoImp(int simID, std::string kernel, double viscosity, int lx, int numberOfSimulations, std::string simulationName, std::vector<std::string> dataToCalcTests)
    : simID(simID), lx(lx), viscosity(viscosity), numberOfSimulations(numberOfSimulations), simulationName(simulationName), dataToCalcTests(dataToCalcTests)
{
    this->kernelName = kernel;
}

//! \}

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
#ifndef SIMULATION_INFO_IMP_H
#define SIMULATION_INFO_IMP_H

#include "SimulationInfo.h"

#include <memory>

class TimeInfo;

class SimulationInfoImp : public SimulationInfo
{
public:
    void setTimeInfo(std::shared_ptr<TimeInfo> timeInfo);

    std::string getKernelName();
    double getViscosity();
    std::string getSimulationName();
    std::string getSimulationParameterString();
    int getLx();
    int getNumberOfSimulations();
    int getSimulationID();
    std::string getRunTimeOutput();
    std::vector<std::string> getDataToCalcTests();

protected:
    SimulationInfoImp() {};
    SimulationInfoImp(int simID, std::string kernel, double viscosity, int lx, int numberOfSimulations, std::string simulationName, std::vector<std::string> dataToCalcTests);

    double viscosity;
    std::string kernelName;
    std::string simulationName;
    std::string simulationParameterString;
    int lx;
    int numberOfSimulations, simID;
    std::shared_ptr<TimeInfo> timeInfo;
    std::vector<std::string> dataToCalcTests;

private:

};
#endif 
//! \}

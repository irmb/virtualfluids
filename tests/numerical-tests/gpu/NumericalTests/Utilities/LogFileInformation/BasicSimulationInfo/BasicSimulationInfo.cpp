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
#include "BasicSimulationInfo.h"

std::shared_ptr<BasicSimulationInfo> BasicSimulationInfo::getNewInstance(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, std::string kernel)
{
    return std::shared_ptr<BasicSimulationInfo>(new BasicSimulationInfo(numberOfTimeSteps, viscosity, basicTimeStepLength, kernel));
}

std::string BasicSimulationInfo::getOutput()
{
    makeCenterHead("Basic Simulation Information");
    oss << "Kernel=" << kernelName << std::endl;
    oss << "NumberOfTimeSteps=" << numberOfTimeSteps << std::endl;
    oss << "Viscosity=" << viscosity << std::endl;
    oss << "BasisTimeStepLength=" << basicTimeStepLength << std::endl;
    oss << std::endl;
    return oss.str();
}

BasicSimulationInfo::BasicSimulationInfo(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, std::string kernel)
  : numberOfTimeSteps(numberOfTimeSteps), viscosity(viscosity), basicTimeStepLength(basicTimeStepLength)
{
    kernelName = kernel;
}
//! \}

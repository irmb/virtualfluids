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
#include "SimulationResults.h"

#include "Utilities/SimulationParameter/SimulationParameter.h"

#define _USE_MATH_DEFINES
#include <math.h>

SimulationResults::SimulationResults(std::shared_ptr<SimulationParameter> simPara) : ResultsImp(simPara->getL0())
{
    this->xNodes = simPara->getLx();
    this->yNodes = 1;
    this->zNodes = simPara->getLz();
    this->numberOfNodes = xNodes*yNodes*zNodes;
    this->timeStepLength = simPara->getTimeStepLength();
    this->numberOfTimeSteps = 0;
}

std::shared_ptr<SimulationResults> SimulationResults::getNewInstance(std::shared_ptr<SimulationParameter> simPara)
{
    return std::shared_ptr<SimulationResults>(new SimulationResults(simPara));
}

void SimulationResults::addTimeStep(unsigned int timeStep, unsigned int time, std::vector<unsigned int> level, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> vx, std::vector<double> vy, std::vector<double> vz, std::vector<double> press, std::vector<double> rho)
{
    this->timeStep.push_back(timeStep);
    this->time.push_back(time);
    this->x.push_back(x);
    this->y.push_back(y);
    this->z.push_back(z);
    this->vx.push_back(vx);
    this->vy.push_back(vy);
    this->vz.push_back(vz);
    this->press.push_back(press);
    this->rho.push_back(rho);
    this->level.push_back(level);
    numberOfTimeSteps++;
}
//! \}

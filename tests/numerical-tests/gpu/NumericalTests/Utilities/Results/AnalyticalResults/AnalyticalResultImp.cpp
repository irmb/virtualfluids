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
#include "AnalyticalResultImp.h"

#include "../SimulationResults/SimulationResults.h"

bool AnalyticalResultsImp::isCalculated()
{
    return calculated;
}

int AnalyticalResultsImp::getL0()
{
    return l0;
}

AnalyticalResultsImp::AnalyticalResultsImp()
{
    calculated = false;
}

AnalyticalResultsImp::AnalyticalResultsImp(int l0)
{
    this->l0 = l0;
    calculated = false;
}

void AnalyticalResultsImp::init(std::shared_ptr<SimulationResults> simResults)
{
    this->xNodes = simResults->getNumberOfXNodes();
    this->yNodes = simResults->getNumberOfYNodes();
    this->zNodes = simResults->getNumberOfZNodes();
    this->numberOfNodes = xNodes*yNodes*zNodes;
    this->timeStepLength = simResults->getTimeStepLength();
    this->numberOfTimeSteps = simResults->getNumberOfTimeSteps();
    this->timeStep = simResults->getTimeSteps();
    this->time = simResults->getTime();
    this->x = simResults->getXNodes();
    this->y = simResults->getYNodes();
    this->z = simResults->getZNodes();
    this->level = simResults->getLevels();
    this->l0 = simResults->getL0();

    this->vx.resize(numberOfTimeSteps);
    this->vy.resize(numberOfTimeSteps);
    this->vz.resize(numberOfTimeSteps);
    this->press.resize(numberOfTimeSteps);
    this->rho.resize(numberOfTimeSteps);

    for (int i = 0; i < numberOfTimeSteps; i++) {
        this->vx.at(i).resize(numberOfNodes);
        this->vy.at(i).resize(numberOfNodes);
        this->vz.at(i).resize(numberOfNodes);
        this->press.at(i).resize(numberOfNodes);
        this->rho.at(i).resize(numberOfNodes);
    }
}
//! \}

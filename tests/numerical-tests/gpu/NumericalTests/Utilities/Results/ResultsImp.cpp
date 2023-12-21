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
#include "ResultsImp.h"

#include <cmath>
#include <iostream>

void logInvalidSimulationData(const std::string &quantity)
{
    std::cout << "done." << std::endl;
    std::cout << "Invalid quantity: " << quantity << std::endl;
    std::cout << "Simulation Result Data contains failure data." << std::endl;
    std::cout << "Testing not possible." << std::endl;
}

bool isValid(const double quantity, const std::string &quantityName)
{
    if (std::isnan(quantity)) {
        logInvalidSimulationData(quantityName);
        return false;
    }
    return true;
}

int ResultsImp::getNumberOfTimeSteps()
{
    return numberOfTimeSteps;
}

std::vector<std::vector<double>> ResultsImp::getVx()
{
    return vx;
}

std::vector<std::vector<double>> ResultsImp::getVy()
{
    return vy;
}

std::vector<std::vector<double>> ResultsImp::getVz()
{
    return vz;
}

int ResultsImp::getNumberOfXNodes()
{
    return xNodes;
}

int ResultsImp::getNumberOfYNodes()
{
    return yNodes;
}

int ResultsImp::getNumberOfZNodes()
{
    return zNodes;
}

std::vector<std::vector<double>> ResultsImp::getXNodes()
{
    return x;
}

std::vector<std::vector<double>> ResultsImp::getYNodes()
{
    return y;
}

std::vector<std::vector<double>> ResultsImp::getZNodes()
{
    return z;
}

int ResultsImp::getTimeStepLength()
{
    return timeStepLength;
}

std::vector<unsigned int> ResultsImp::getTimeSteps()
{
    return timeStep;
}

std::vector<int> ResultsImp::getTime()
{
    return time;
}

std::vector<std::vector<unsigned int>> ResultsImp::getLevels()
{
    return level;
}

std::vector<std::vector<double>> ResultsImp::getPress()
{
    return press;
}

std::vector<std::vector<double>> ResultsImp::getRho()
{
    return rho;
}

int ResultsImp::getL0()
{
    return l0;
}

bool ResultsImp::checkYourData()
{
    std::cout << "checking Simulation Results Data...";
    for (int i = 0; i < vx.size(); i++) {
        for (int j = 0; j < vx.at(i).size(); j++) {
            bool valid = isValid(vx.at(i).at(j), "Vx") && isValid(vy.at(i).at(j), "Vy") &&
                         isValid(vz.at(i).at(j), "Vz") && isValid(rho.at(i).at(j), "Rho") &&
                         isValid(press.at(i).at(j), "Pressure");

            if (!valid)
                return false;
        }
    }
    std::cout << "done." << std::endl;
    std::cout << "Simulation Result Data contains no failure data." << std::endl;
    return true;
}

ResultsImp::ResultsImp(int l0)
{
    this->l0 = l0;
}

//! \}

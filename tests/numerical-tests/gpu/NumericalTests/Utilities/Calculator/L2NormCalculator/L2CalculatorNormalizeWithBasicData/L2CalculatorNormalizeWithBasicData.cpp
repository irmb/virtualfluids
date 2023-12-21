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
#include "L2CalculatorNormalizeWithBasicData.h"

#include <cmath>

std::shared_ptr<L2NormCalculator> L2CalculatorNormalizeWithBasicData::getInstance()
{
    static std::shared_ptr<L2NormCalculator> uniqueInstance;
    if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<L2NormCalculator>(new L2CalculatorNormalizeWithBasicData());
    return uniqueInstance;
}

double L2CalculatorNormalizeWithBasicData::calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz, double l0)
{
    double counter = calcCounter(basicData, divergentData, level, lx, lz);
    double denominator = 0.0;
    for (int i = 0; i < basicData.size(); i++) {
        double area = (1 / std::pow(2.0, level.at(i))) * (1 / std::pow(2.0, level.at(i)));
        denominator += (basicData.at(i)*basicData.at(i)) * area;
    }
    if (equalDouble(denominator, 0.0))
        return -1.0;

    return std::sqrt(counter / denominator);
}

L2CalculatorNormalizeWithBasicData::L2CalculatorNormalizeWithBasicData() : L2NormCalculatorImp("Test could not run. BasicData is zero. Normalization of the data is not possible.")
{

}

//! \}

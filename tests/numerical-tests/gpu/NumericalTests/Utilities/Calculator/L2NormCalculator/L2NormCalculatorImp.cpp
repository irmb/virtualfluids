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
#include "L2NormCalculatorImp.h"

#include "Utilities/Calculator/AlmostEquals.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

L2NormCalculatorImp::L2NormCalculatorImp(std::string errorMessage) : errorMessage(errorMessage)
{
}

bool L2NormCalculatorImp::equalDouble(double num1, double num2)
{
    const FloatingPoint<double> lhs(num1), rhs(num2);

    if (lhs.AlmostEquals(rhs))
        return true;
    return false;
}

double L2NormCalculatorImp::calcCounter(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz)
{
    double counter = 0.0;
    for (int i = 0; i < basicData.size(); i++) {
        double area = (1 / pow(2.0, level.at(i))) * (1 / pow(2.0, level.at(i)));
        counter += ((divergentData.at(i) - basicData.at(i))*(divergentData.at(i) - basicData.at(i))) * area;
    }
    return counter;
}

std::string L2NormCalculatorImp::getErrorMessage()
{
    return errorMessage;
}

L2NormCalculatorImp::L2NormCalculatorImp()
{
    
}
//! \}

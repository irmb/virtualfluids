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
#include "PhiLogFileDataImp.h"

std::shared_ptr<PhiLogFileDataImp> PhiLogFileDataImp::getNewInstance()
{
    return std::shared_ptr<PhiLogFileDataImp>(new PhiLogFileDataImp());
}

std::vector<double> PhiLogFileDataImp::getBasicGridLengths()
{
    return basicGridLengths;
}

int PhiLogFileDataImp::getStartTimeStepCalculation()
{
    return startTimeStepCalculation;
}

int PhiLogFileDataImp::getEndTimeStepCalculation()
{
    return endTimeStepCalculation;
}

std::string PhiLogFileDataImp::getDataToCalc()
{
    return dataToCalc;
}

std::vector<double> PhiLogFileDataImp::getPhiDiff()
{
    return phiDiff;
}

std::vector<std::vector<double>> PhiLogFileDataImp::getOrderOfAccuracy()
{
    return orderOfAccuracy;
}

void PhiLogFileDataImp::setBasicGridLengths(std::vector<double> basicGridLengths)
{
    this->basicGridLengths = basicGridLengths;
}

void PhiLogFileDataImp::setStartTimeStepCalculation(int startTimeStepCalculation)
{
    this->startTimeStepCalculation = startTimeStepCalculation;
}

void PhiLogFileDataImp::setEndTimeStepCalculation(int endTimeStepCalculation)
{
    this->endTimeStepCalculation = endTimeStepCalculation;
}

void PhiLogFileDataImp::setDataToCalc(std::string dataToCalc)
{
    this->dataToCalc = dataToCalc;
}

void PhiLogFileDataImp::setPhiDiff(std::vector<double> phiDiff)
{
    this->phiDiff = phiDiff;
}

void PhiLogFileDataImp::setOrderOfAccuracy(std::vector<std::vector<double>> orderOfAccuracy)
{
    this->orderOfAccuracy = orderOfAccuracy;
}

PhiLogFileDataImp::PhiLogFileDataImp()
{
}

//! \}

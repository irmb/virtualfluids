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
#include "NyLogFileDataImp.h"

std::shared_ptr<NyLogFileDataImp> NyLogFileDataImp::getNewInstance()
{
    return std::shared_ptr<NyLogFileDataImp>(new NyLogFileDataImp());
}

std::vector<double> NyLogFileDataImp::getBasicGridLengths()
{
    return basicGridLengths;
}

int NyLogFileDataImp::getStartTimeStepCalculation()
{
    return startTimeStepCalculation;
}

int NyLogFileDataImp::getEndTimeStepCalculation()
{
    return endTimeStepCalculation;
}

std::string NyLogFileDataImp::getDataToCalc()
{
    return dataToCalc;
}

std::vector<double> NyLogFileDataImp::getNy()
{
    return ny;
}

std::vector<double> NyLogFileDataImp::getNyDiff()
{
    return nyDiff;
}

std::vector<std::vector<double>> NyLogFileDataImp::getOrderOfAccuracy()
{
    return orderOfAccuracy;
}

void NyLogFileDataImp::setBasicGridLengths(std::vector<double> basicGridLengths)
{
    this->basicGridLengths = basicGridLengths;
}

void NyLogFileDataImp::setStartTimeStepCalculation(int startTimeStepCalculation)
{
    this->startTimeStepCalculation = startTimeStepCalculation;
}

void NyLogFileDataImp::setEndTimeStepCalculation(int endTimeStepCalculation)
{
    this->endTimeStepCalculation = endTimeStepCalculation;
}

void NyLogFileDataImp::setDataToCalc(std::string dataToCalc)
{
    this->dataToCalc = dataToCalc;
}

void NyLogFileDataImp::setNy(std::vector<double> ny)
{
    this->ny = ny;
}

void NyLogFileDataImp::setNyDiff(std::vector<double> nyDiff)
{
    this->nyDiff = nyDiff;
}

void NyLogFileDataImp::setOrderOfAccuracy(std::vector<std::vector<double>> orderOfAccuracy)
{
    this->orderOfAccuracy = orderOfAccuracy;
}

NyLogFileDataImp::NyLogFileDataImp()
{
}

//! \}

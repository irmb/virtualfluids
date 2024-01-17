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
#include "L2NormLogFileDataImp.h"

std::shared_ptr<L2NormLogFileDataImp> L2NormLogFileDataImp::getNewInstance()
{
    return std::shared_ptr<L2NormLogFileDataImp>(new L2NormLogFileDataImp());
}

std::vector<double> L2NormLogFileDataImp::getBasicGridLengths()
{
    return basicGridLengths;
}

std::string L2NormLogFileDataImp::getDataToCalc()
{
    return dataToCalc;
}

std::string L2NormLogFileDataImp::getNormalizeData()
{
    return normalizeData;
}

int L2NormLogFileDataImp::getBasicTimeStep()
{
    return basicTimeStep;
}

int L2NormLogFileDataImp::getDivergentTimeStep()
{
    return divergentTimeStep;
}

std::vector<double> L2NormLogFileDataImp::getL2NormForBasicTimeStep()
{
    return l2NormForBasicTimeStep;
}

std::vector<double> L2NormLogFileDataImp::getL2NormForDivergentTimeStep()
{
    return l2NormForDivergentTimeStep;
}

std::vector<double> L2NormLogFileDataImp::getL2NormDiff()
{
    return l2NormDiff;
}

void L2NormLogFileDataImp::setBasicGridLengths(std::vector<double> basicGridLengths)
{
    this->basicGridLengths = basicGridLengths;
}

void L2NormLogFileDataImp::setDataToCalc(std::string dataToCalc)
{
    this->dataToCalc = dataToCalc;
}

void L2NormLogFileDataImp::setNormalizeData(std::string normalizeData)
{
    this->normalizeData = normalizeData;
}

void L2NormLogFileDataImp::setBasicTimeStep(int basicTimeStep)
{
    this->basicTimeStep = basicTimeStep;
}

void L2NormLogFileDataImp::setDivergentTimeStep(int divergentTimeStep)
{
    this->divergentTimeStep = divergentTimeStep;
}

void L2NormLogFileDataImp::setL2NormForBasicTimeStep(std::vector<double> l2NormForBasicTimeStep)
{
    this->l2NormForBasicTimeStep = l2NormForBasicTimeStep;
}

void L2NormLogFileDataImp::setL2NormForDivergentTimeStep(std::vector<double> l2NormForDivergentTimeStep)
{
    this->l2NormForDivergentTimeStep = l2NormForDivergentTimeStep;
}

void L2NormLogFileDataImp::setL2NormDiff(std::vector<double> l2NormDiff)
{
    this->l2NormDiff = l2NormDiff;
}

L2NormLogFileDataImp::L2NormLogFileDataImp()
{
}

L2NormLogFileDataImp::~L2NormLogFileDataImp()
{
}

//! \}

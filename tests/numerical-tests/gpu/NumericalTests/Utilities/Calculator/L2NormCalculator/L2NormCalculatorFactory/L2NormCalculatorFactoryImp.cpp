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
#include "L2NormCalculatorFactoryImp.h"

#include "Utilities/Calculator/L2NormCalculator/L2CalculatorNormalizeWithAmplitude/L2CalculatorNormalizeWithAmplitude.h"
#include "Utilities/Calculator/L2NormCalculator/L2CalculatorNormalizeWithBasicData/L2CalculatorNormalizeWithBasicData.h"

std::shared_ptr<L2NormCalculatorFactory> L2NormCalculatorFactoryImp::getInstance()
{
    static std::shared_ptr<L2NormCalculatorFactory> uniqueInstance;
    if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<L2NormCalculatorFactory>(new L2NormCalculatorFactoryImp());
    return uniqueInstance;
}

std::shared_ptr<L2NormCalculator> L2NormCalculatorFactoryImp::makeL2NormCalculator(std::string type)
{
    if(type == "Amplitude")
        return L2CalculatorNormalizeWithAmplitude::getInstance();
    if (type == "BasicData")
        return L2CalculatorNormalizeWithBasicData::getInstance();
}

L2NormCalculatorFactoryImp::L2NormCalculatorFactoryImp()
{
}

//! \}

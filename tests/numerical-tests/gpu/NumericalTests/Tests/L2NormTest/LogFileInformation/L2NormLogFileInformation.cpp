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
#include "L2NormLogFileInformation.h"

#include "Tests/L2NormTest/L2NormTest.h"
#include "Tests/L2NormTest/L2NormTestParameterStruct.h"

#include <iomanip>
#include <sstream>

std::shared_ptr<L2NormInformation> L2NormInformation::getNewInstance(std::vector<std::shared_ptr<L2NormTest> > tests, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::string> dataToCalcTests)
{
    return std::shared_ptr<L2NormInformation>(new L2NormInformation(tests, testParameter, dataToCalcTests));
}

std::string L2NormInformation::getOutput()
{
    std::ostringstream headName;
    headName << " L2Norm Test";
    makeCenterHead(headName.str());

    oss << "BasicTimeStep_L2Norm=" << basicTimeStep << std::endl;
    oss << "DivergentTimeStep_L2Norm=" << divergentTimeStep << std::endl;
    oss << "DataToCalc_L2Norm=\"";
    for (int i = 0; i < dataToCalc.size(); i++) {
        oss << dataToCalc.at(i);
        if (i < dataToCalc.size() - 1)
            oss << " ";
        else
            oss << "\"" << std::endl;
    }
    oss << "NormalizeData_L2Norm=\"";
    for (int i = 0; i < normalizeData.size(); i++) {
        oss << normalizeData.at(i);
        if (i < normalizeData.size() - 1)
            oss << " ";
        else
            oss << "\"" << std::endl;
    }
    oss << std::endl;

    std::ostringstream failMessage;
    failMessage << "FailTests_L2Norm=\"";
    for (int i = 0; i < tests.size(); i++) {
        if (tests.at(i)->getTestStatus() == passed || tests.at(i)->getTestStatus() == failed)
            oss << tests.at(i)->getLogFileOutput();
        if (tests.at(i)->getTestStatus() == test_error || tests.at(i)->getTestStatus() == simulationCrashed)
            failMessage << tests.at(i)->getErrorLogFileOutput() << " ";
    }
    std::string fail = failMessage.str();
    if(fail.back() == ' ')
        fail = fail.substr(0, fail.size() - 1);
    failMessage.str(std::string());
    failMessage << fail << "\"";
    oss << failMessage.str() << std::endl << std::endl;        
        
    return oss.str();
}

L2NormInformation::L2NormInformation(std::vector<std::shared_ptr<L2NormTest> > tests, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::string> dataToCalcTests) : tests(tests), dataToCalc(dataToCalcTests)
{
    basicTimeStep = testParameter->basicTimeStep;
    divergentTimeStep = testParameter->divergentTimeStep;
    normalizeData = testParameter->normalizeData;
}
//! \}

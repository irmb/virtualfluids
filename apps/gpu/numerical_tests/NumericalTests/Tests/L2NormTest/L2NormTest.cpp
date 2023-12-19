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
#include "L2NormTest.h"

#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities/SimulationInfo/SimulationInfo.h"

#include "Tests/L2NormTest/PostProcessingStrategy/PostProcessingStrategyL2NormTest.h"
#include "Tests/L2NormTest/L2NormTestParameterStruct.h"

#include <iomanip>
#include <sstream>

std::shared_ptr<L2NormTest> L2NormTest::getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate, double maxL2NormDiff, std::string normalizeData)
{
    return std::shared_ptr<L2NormTest>(new L2NormTest(colorOutput, testParameter, dataToCalculate, maxL2NormDiff, normalizeData));
}

void L2NormTest::update()
{
    TestImp::update();
}

void L2NormTest::addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormPostProcessingStrategy> postProStrategy)
{
    TestImp::addSimulation(sim, simInfo, postProStrategy);
    l2NormPostProStrategies.push_back(postProStrategy);
}

void L2NormTest::evaluate()
{
    std::vector<double> results = l2NormPostProStrategies.at(0)->getL2Norm(dataToCalculate, normalizeData);
        
    resultBasicTimestep = results.at(0);
    resultDivergentTimeStep = results.at(1);
    diffL2Norm = resultDivergentTimeStep - resultBasicTimestep;

    if (resultBasicTimestep < 0 || resultDivergentTimeStep < 0) {
        testStatus = test_error;
    }
    else
    {
        testPassed = maxL2NormDiff > diffL2Norm;
        if (testPassed)
            testStatus = passed;
        else
            testStatus = failed;
    }
    
    makeConsoleOutput();
}

std::string L2NormTest::getLogFileOutput()
{
    std::ostringstream oss;
    oss << "L2Norm_BasicTimeStep_L" << l2NormPostProStrategies.at(0)->getNumberOfXNodes() << "_" << dataToCalculate << "_" << normalizeData << "=" << resultBasicTimestep << std::endl;
    oss << "L2Norm_DivergentTimeStep_L" << l2NormPostProStrategies.at(0)->getNumberOfXNodes() << "_" << dataToCalculate << "_" << normalizeData << "=" << resultDivergentTimeStep << std::endl;
    oss << "L2Norm_Diff_L" << l2NormPostProStrategies.at(0)->getNumberOfXNodes() << "_" << dataToCalculate << "_" << normalizeData << "=" << diffL2Norm << std::endl << std::endl;
    return oss.str();
}

std::string L2NormTest::getErrorLogFileOutput()
{
    std::ostringstream oss;
    oss << "L" << l2NormPostProStrategies.at(0)->getNumberOfXNodes() << "_" << dataToCalculate << "_" << normalizeData;
    return oss.str();
}

L2NormTest::L2NormTest(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate, double maxL2NormDiff, std::string normalizeData)
    : TestImp(colorOutput), dataToCalculate(dataToCalculate), normalizeData(normalizeData)
{
    basicTimeStep = testParameter->basicTimeStep;
    divergentTimeStep = testParameter->divergentTimeStep;
    this->maxL2NormDiff = maxL2NormDiff;
}

std::vector<std::string> L2NormTest::buildTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    oss << "L2Norm BasicTimeStep: " << resultBasicTimestep;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "L2Norm DivergentTimeStep: " << resultDivergentTimeStep;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "L2NormDiff: " << diffL2Norm;
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}

std::vector<std::string> L2NormTest::buildBasicTestOutput()
{
    std::vector<std::string> output;
    std::ostringstream oss;

    output.push_back("L2 Norm Test");

    oss << "Kernel: " << simInfos.at(0)->getKernelName();
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "Viscosity: " << simInfos.at(0)->getViscosity();
    output.push_back(oss.str());
    oss.str(std::string());

    output.push_back(oss.str());

    oss << simInfos.at(0)->getSimulationName();
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "L: " << simInfos.at(0)->getLx() << simInfos.at(0)->getSimulationParameterString();
    output.push_back(oss.str());
    oss.str(std::string());

    output.push_back(oss.str());

    oss << "DataToCalculate: " << dataToCalculate;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "NormalizeData: " << normalizeData;
    output.push_back(oss.str());
    oss.str(std::string());

    output.push_back(oss.str());

    oss << "BasicTimeStep: " << basicTimeStep;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "DivergentTimeStep: " << divergentTimeStep;
    output.push_back(oss.str());
    oss.str(std::string());

    output.push_back(oss.str());

    return output;
}

std::vector<std::string> L2NormTest::buildErrorTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    oss << "Error Message: " << l2NormPostProStrategies.at(0)->getErrorMessage(normalizeData);
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}

//! \}

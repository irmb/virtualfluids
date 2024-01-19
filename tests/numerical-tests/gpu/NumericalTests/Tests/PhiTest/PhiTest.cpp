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
#include "PhiTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities/TestSimulation/TestSimulation.h"
#include "Utilities/SimulationInfo/SimulationInfo.h"

#include "Tests/PhiTest/PostProcessingStrategy/PhiTestPostProcessingStrategy.h"
#include "Tests/PhiTest/PhiTestParameterStruct.h"

#include "basics/DataTypes.h"

#include <iomanip>
#include <cmath>

std::shared_ptr<PhiTest> PhiTest::getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiTestParameterStruct> testPara, std::string dataToCalculate)
{
    return std::shared_ptr<PhiTest>(new PhiTest(colorOutput, viscosity, testPara, dataToCalculate));
}

void PhiTest::evaluate()
{
    for (uint i = 0; i < postProStrategies.size(); i++)
        phiDiff.push_back(postProStrategies.at(i)->getPhiDiff(dataToCalculate));
    
    orderOfAccuracy = calcOrderOfAccuracy(phiDiff);
    testStatus = checkTestPassed(orderOfAccuracy);
    
    makeConsoleOutput();
}

void PhiTest::update()
{
    TestImp::update();
}

void PhiTest::addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<PhiTestPostProcessingStrategy> postProStrategy)
{
    TestImp::addSimulation(sim, simInfo, postProStrategy);
    postProStrategies.push_back(postProStrategy);
    lx.push_back(postProStrategy->getNumberOfXNodes());
}

std::string PhiTest::getDataToCalculate()
{
    return dataToCalculate;
}

std::vector<int> PhiTest::getLx()
{
    std::vector<int> lxINT;
    for (uint i = 0; i < lx.size(); i++)
        lxINT.push_back((int)lx.at(i));
    return lxINT;
}

std::vector<double> PhiTest::getPhiDiff()
{
    return phiDiff;
}

double PhiTest::getOrderOfAccuracy()
{
    return orderOfAccuracy;
}

PhiTest::PhiTest(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiTestParameterStruct> testPara, std::string dataToCalculate)
    : TestImp(colorOutput), dataToCalculate(dataToCalculate)
{
    minOrderOfAccuracy = testPara->minOrderOfAccuracy;
    startStepCalculation = testPara->startTimeStepCalculation;
    endStepCalculation = testPara->endTimeStepCalculation;

    lx.resize(0);
    phiDiff.resize(0);
}

double PhiTest::calcOrderOfAccuracy(std::vector<double> data)
{
    double ooa = std::log(data.at(0) / data.at(1)) / std::log(lx.at(1) / lx.at(0));
    
    return ooa;
}

TestStatus PhiTest::checkTestPassed(double orderOfAccuracy)
{
    if (orderOfAccuracy > minOrderOfAccuracy)
        return passed;
    else
        return failed;
}

std::vector<std::string> PhiTest::buildTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    for (uint i = 0; i < phiDiff.size(); i++) {
        oss << "PhiDiff" << simInfos.at(i)->getLx() << ": " << phiDiff.at(i);
        output.push_back(oss.str());
        oss.str(std::string());
    }
    oss << "OrderOfAccuracy: " << orderOfAccuracy;
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}

std::vector<std::string> PhiTest::buildBasicTestOutput()
{
    std::vector<std::string> output;
    std::ostringstream oss;

    output.push_back("Phi Test");

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

    for (uint i = 0; i < simInfos.size(); i++) {
        oss << "L: " << std::setfill(' ') << std::right << std::setw(4) << simInfos.at(i)->getLx() << simInfos.at(i)->getSimulationParameterString();
        output.push_back(oss.str());
        oss.str(std::string());
    }

    output.push_back(oss.str());

    oss << "DataToCalculate: " << dataToCalculate;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "StartTimeStep: " << startStepCalculation;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "EndTimeStep: " << endStepCalculation;
    output.push_back(oss.str());
    oss.str(std::string());

    output.push_back(oss.str());

    return output;
}

std::vector<std::string> PhiTest::buildErrorTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    oss << "Error Message: ";
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}

//! \}

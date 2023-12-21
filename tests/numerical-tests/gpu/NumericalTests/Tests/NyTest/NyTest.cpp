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
#include "NyTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities/TestSimulation/TestSimulation.h"
#include "Utilities/SimulationInfo/SimulationInfo.h"

#include "Tests/NyTest/PostProcessingStrategy/NyTestPostProcessingStrategy.h"
#include "Tests/NyTest/NyTestParameterStruct.h"

#include <iomanip>
#include <cmath>

std::shared_ptr<NyTest> NyTest::getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<NyTestParameterStruct> testPara, std::string dataToCalculate)
{
    return std::shared_ptr<NyTest>(new NyTest(colorOutput, viscosity, testPara, dataToCalculate));
}

void NyTest::evaluate()
{
    for (int i = 0; i < postProStrategies.size(); i++)
        ny.push_back(postProStrategies.at(i)->getNy(dataToCalculate));
    
    if (checkNy(ny)) {
        nyDiff = calcNyDiff(ny);
        orderOfAccuracy = calcOrderOfAccuracy(nyDiff);
        testStatus = checkTestPassed(orderOfAccuracy);
    }
    else
        testStatus = test_error;
    

    makeConsoleOutput();
}

void NyTest::update()
{
    TestImp::update();
}

void NyTest::addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<NyTestPostProcessingStrategy> postProStrategy)
{
    TestImp::addSimulation(sim, simInfo, postProStrategy);
    postProStrategies.push_back(postProStrategy);
    lx.push_back(postProStrategy->getNumberOfXNodes());
}

std::string NyTest::getDataToCalculate()
{
    return dataToCalculate;
}

std::vector<int> NyTest::getLx()
{
    std::vector<int> lxINT;
    for (int i = 0; i < lx.size(); i++)
        lxINT.push_back((int)lx.at(i));
    return lxINT;
}

std::vector<double> NyTest::getNy()
{
    return ny;
}

std::vector<double> NyTest::getNyDiff()
{
    return nyDiff;
}

double NyTest::getOrderOfAccuracyNyDiff()
{
    return orderOfAccuracy;
}

NyTest::NyTest(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<NyTestParameterStruct> testPara, std::string dataToCalculate)
    : TestImp(colorOutput), viscosity(viscosity), dataToCalculate(dataToCalculate)
{
    minOrderOfAccuracy = testPara->minOrderOfAccuracy;
    startStepCalculation = testPara->startTimeStepCalculation;
    endStepCalculation = testPara->endTimeStepCalculation;

    lx.resize(0);
    nyDiff.resize(0);
}

double NyTest::calcOrderOfAccuracy(std::vector<double> data)
{
    double ooa = std::log(data.at(0) / data.at(1)) / std::log(lx.at(1) / lx.at(0));
    
    return ooa;
}

TestStatus NyTest::checkTestPassed(double orderOfAccuracy)
{
    if (orderOfAccuracy > minOrderOfAccuracy)
        return passed;
    else
        return failed;
}

bool NyTest::checkNy(std::vector<double> ny)
{
    for(int i = 0; i < ny.size(); i++)
        if(ny.at(i) < 0.0)
            return false;
    return true;
}

std::vector<double> NyTest::calcNyDiff(std::vector<double> ny)
{
    std::vector<double> results;
    for (int i = 0; i < ny.size(); i++)
        results.push_back(std::fabs((ny.at(i) - viscosity) / viscosity));
    return results;
}

std::vector<std::string> NyTest::buildTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    for (int i = 0; i < ny.size(); i++) {
        oss << "Ny" << simInfos.at(i)->getLx() << ": " << ny.at(i);
        output.push_back(oss.str());
        oss.str(std::string());

        oss << "NyDiff" << simInfos.at(i)->getLx() << ": " << nyDiff.at(i);
        output.push_back(oss.str());
        oss.str(std::string());
    }
    oss << "OrderOfAccuracy: " << orderOfAccuracy;
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}

std::vector<std::string> NyTest::buildBasicTestOutput()
{
    std::vector<std::string> output;
    std::ostringstream oss;

    output.push_back("Ny Test");

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

    for (int i = 0; i < simInfos.size(); i++) {
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

std::vector<std::string> NyTest::buildErrorTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    oss << "Error Message: Ny < 0";
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}

std::vector<std::string> NyTest::buildSimulationFailedTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    oss << "Simulation crashed!";
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}

//! \}

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
#include "L2NormTestBetweenKernels.h"

#include "Utilities/Calculator/L2NormCalculator/L2NormCalculator.h"
#include "Utilities/Calculator/L2NormCalculator/L2NormCalculatorFactory/L2NormCalculatorFactory.h"
#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities/TestSimulation/TestSimulation.h"
#include "Utilities/SimulationInfo/SimulationInfo.h"
#include "Tests/L2NormTestBetweenKernels/PostProcessingStrategy/L2NormBetweenKernelPostProcessingStrategy.h"

#include <iomanip>

std::shared_ptr<L2NormTestBetweenKernels> L2NormTestBetweenKernels::getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep, std::string normalizeWith, std::shared_ptr<L2NormCalculatorFactory> factory)
{
    return std::shared_ptr<L2NormTestBetweenKernels>(new L2NormTestBetweenKernels(colorOutput, dataToCalculate, timeStep, normalizeWith, factory));
}

void L2NormTestBetweenKernels::update()
{
    TestImp::update();
}

void L2NormTestBetweenKernels::evaluate()
{
    basicPostProcessingStrategy->evaluate();
    divergentPostProcessingStrategy->evaluate();

    int tS = calcTimeStepInResults(timeStep);

    basicL2Result = basicPostProcessingStrategy->getL2Norm(dataToCalculate, normalizeData, timeStep);
    divergentL2Result = divergentPostProcessingStrategy->getL2Norm(dataToCalculate, normalizeData, timeStep);

    if (dataToCalculate == "Vx") 
        resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getVx().at(tS), divergentSimResults->getVx().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getL0());
    if (dataToCalculate == "Vy") 
        resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getVy().at(tS), divergentSimResults->getVy().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getL0());
    if (dataToCalculate == "Vz")
        resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getVz().at(tS), divergentSimResults->getVz().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getL0());
    if (dataToCalculate == "Press")
        resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getPress().at(tS), divergentSimResults->getPress().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getL0());
    if (dataToCalculate == "Rho")
        resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getRho().at(tS), divergentSimResults->getRho().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getL0());
    
    
    
    if (basicL2Result < 0 || divergentL2Result < 0 || resultL2ToBasicKernel < 0)
        testStatus = test_error;
    else if (basicL2Result <= divergentL2Result)
        testStatus = passed;
    else
        testStatus = failed;

    makeConsoleOutput();
}

std::string L2NormTestBetweenKernels::getLogFileOutput()
{
    std::ostringstream oss;
    oss << "L2Norm_BasicKernel_"     << "L" << basicPostProcessingStrategy->getNumberOfXNodes() << "_"<< dataToCalculate << "_TimeStep_" << timeStep << "_" << normalizeData << "=" << basicL2Result << std::endl;
    oss << "L2Norm_DivergentKernel_" << "L" << basicPostProcessingStrategy->getNumberOfXNodes() << "_"<< dataToCalculate << "_TimeStep_" << timeStep << "_" << normalizeData << "=" << divergentL2Result << std::endl;
    oss << "L2Norm_Between_Kernels_" << "L" << basicPostProcessingStrategy->getNumberOfXNodes() << "_"<< dataToCalculate << "_TimeStep_" << timeStep << "_" << normalizeData << "=" << resultL2ToBasicKernel << std::endl << std::endl;
    return oss.str();
}

std::string L2NormTestBetweenKernels::getErrorLogFileOutput()
{
    std::ostringstream oss;
    oss << "L" << basicPostProcessingStrategy->getNumberOfXNodes() << "_"<< dataToCalculate << "_TimeStep_" << timeStep << "_" << normalizeData;
    return oss.str();
}

double L2NormTestBetweenKernels::getBasicL2Result()
{
    return basicL2Result;
}

void L2NormTestBetweenKernels::setBasicSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy)
{
    TestImp::addSimulation(sim, simInfo, postProcessingStrategy);
    this->basicSim = sim;
    this->basicSimInfo = simInfo;
    this->basicPostProcessingStrategy = postProcessingStrategy;
    this->basicSimResults = basicPostProcessingStrategy->getSimulationResult();
}

void L2NormTestBetweenKernels::setDivergentKernelSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy)
{
    TestImp::addSimulation(sim, simInfo, postProcessingStrategy);
    this->divergentSim = sim;
    this->divergentSimInfo = simInfo;
    this->divergentPostProcessingStrategy = postProcessingStrategy;
    this->divergentSimResults = divergentPostProcessingStrategy->getSimulationResult();
}

L2NormTestBetweenKernels::L2NormTestBetweenKernels(std::shared_ptr<ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep, std::string normalizeWith, std::shared_ptr<L2NormCalculatorFactory> factory)
    : TestImp(colorOutput), timeStep(timeStep), dataToCalculate(dataToCalculate), normalizeData(normalizeWith)
{
    l2Normcalculator = factory->makeL2NormCalculator(normalizeWith);
}

int L2NormTestBetweenKernels::calcTimeStepInResults(unsigned int timeStep)
{
    for (int i = 0; i < basicSimResults->getTimeSteps().size(); i++) {
        if (timeStep == basicSimResults->getTimeSteps().at(i))
            return basicSimResults->getTimeSteps().at(i);
    }
}

std::vector<std::string> L2NormTestBetweenKernels::buildTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    oss << "L2Norm BasicKernel: " << basicL2Result;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "L2Norm DivergentKernel: " << divergentL2Result;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "L2NormDiff: " << resultL2ToBasicKernel;
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}

std::vector<std::string> L2NormTestBetweenKernels::buildBasicTestOutput()
{
    std::vector<std::string> output;
    std::ostringstream oss;

    output.push_back("L2 Norm Between Kernels Test");

    oss << "Basic Kernel: " << basicSimInfo->getKernelName();
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "Divergent Kernel: " << divergentSimInfo->getKernelName();
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "Viscosity: " << basicSimInfo->getViscosity();
    output.push_back(oss.str());
    oss.str(std::string());

    output.push_back(oss.str());

    oss << basicSimInfo->getSimulationName();
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "L: " << basicSimInfo->getLx() << basicSimInfo->getSimulationParameterString();
    output.push_back(oss.str());
    oss.str(std::string());

    output.push_back(oss.str());

    oss << "DataToCalculate: " << dataToCalculate;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "NormalizeData: " << normalizeData;
    output.push_back(oss.str());
    oss.str(std::string());

    oss << "TimeStep: " << timeStep;
    output.push_back(oss.str());
    oss.str(std::string());

    output.push_back(oss.str());

    return output;
}

std::vector<std::string> L2NormTestBetweenKernels::buildErrorTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    oss << "Error Message: " << basicPostProcessingStrategy->getErrorMessage(normalizeData);
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}

//! \}

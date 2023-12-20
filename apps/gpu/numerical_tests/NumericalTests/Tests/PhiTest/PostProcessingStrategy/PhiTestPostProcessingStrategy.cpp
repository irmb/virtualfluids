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
#include "PhiTestPostProcessingStrategy.h"

#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"
#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"

#include "Tests/PhiTest/PhiTestParameterStruct.h"

std::shared_ptr<PhiTestPostProcessingStrategy> PhiTestPostProcessingStrategy::getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests)
{
    return std::shared_ptr<PhiTestPostProcessingStrategy>(new PhiTestPostProcessingStrategy(simResult, analyticalResult, testPara, dataToCalcTests));
}

PhiTestPostProcessingStrategy::PhiTestPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests)
    : PostProcessingStrategyImp(simResult), analyticalResult(analyticalResult), dataToCalculate(dataToCalcTests)
{
    startTimeStepCalculation = testPara->startTimeStepCalculation;
    endTimeStepCalculation = testPara->endTimeStepCalculation;
    phiDiff.resize(dataToCalculate.size());

    isEvaluated = false;
    fftCalculator = FFTCalculator::getInstance();
}

std::vector<std::vector<double> > PhiTestPostProcessingStrategy::reduceDataToTimeSteps(std::vector<std::vector<double> > data)
{
    std::vector<int> timeStepsToDelete;

    for (int i = simResult->getTimeSteps().size() - 1; i >= 0; i--) {
        if (simResult->getTimeSteps().at(i) > endTimeStepCalculation)
            timeStepsToDelete.push_back(i);
        if (simResult->getTimeSteps().at(i) < startTimeStepCalculation)
            timeStepsToDelete.push_back(i);
    }

    for (int i = 0; i < timeStepsToDelete.size(); i++)
        data.erase(data.begin() + timeStepsToDelete.at(i));

    return data;
}

void PhiTestPostProcessingStrategy::evaluate()
{
    if (!isEvaluated) {
        bool transpose = false;
        int xNodes = simResult->getNumberOfXNodes();
        int zNodes = simResult->getNumberOfZNodes();
        int timeStepLength = simResult->getTimeStepLength();
        for (int i = 0; i < dataToCalculate.size(); i++) {
            std::vector<std::vector<double>> basicData;
            if (dataToCalculate.at(i) == "Vx")
                basicData = simResult->getVx();
            if (dataToCalculate.at(i) == "Vy")
                basicData = simResult->getVy();
            if (dataToCalculate.at(i) == "Vz") 
                basicData = simResult->getVz();
            if (dataToCalculate.at(i) == "Press")
                basicData = simResult->getPress();
            if (dataToCalculate.at(i) == "Rho")
                basicData = simResult->getRho();

            std::vector<std::vector<double>> dataForCalculation;
            dataForCalculation = reduceDataToTimeSteps(basicData);
            phiDiff.at(i) = fftCalculator->calcPhiDiff(dataForCalculation, transpose, xNodes, zNodes, timeStepLength);
        }
        isEvaluated = true;
    }
}

double PhiTestPostProcessingStrategy::getPhiDiff(std::string dataToCalc)
{
    for (int i = 0; i < dataToCalculate.size(); i++)
        if(dataToCalculate.at(i) == dataToCalc)
            return phiDiff.at(i);
}

//! \}

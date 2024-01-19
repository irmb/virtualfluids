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
#include "NyTestPostProcessingStrategy.h"

#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"
#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"

#include "Tests/NyTest/NyTestParameterStruct.h"


std::shared_ptr<NyTestPostProcessingStrategy> NyTestPostProcessingStrategy::getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<NyTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests)
{
    return std::shared_ptr<NyTestPostProcessingStrategy>(new NyTestPostProcessingStrategy(simResult, analyticalResult, testPara, dataToCalcTests));
}

NyTestPostProcessingStrategy::NyTestPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<NyTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests)
    : PostProcessingStrategyImp(simResult), analyticalResult(analyticalResult), dataToCalculate(dataToCalcTests)
{
    startTimeStepCalculation = testPara->startTimeStepCalculation;
    endTimeStepCalculation = testPara->endTimeStepCalculation;
    ny.resize(dataToCalculate.size());

    isEvaluated = false;
    fftCalculator = FFTCalculator::getInstance();
}

std::vector<std::vector<double> > NyTestPostProcessingStrategy::reduceDataToTimeSteps(std::vector<std::vector<double> > data, unsigned int startTimeStep, unsigned int endTimeStep)
{
    std::vector<int> timeStepsToDelete;

    for (int i = simResult->getTimeSteps().size() - 1; i >= 0; i--) {
        if (simResult->getTimeSteps().at(i) > endTimeStep)
            timeStepsToDelete.push_back(i);
        if (simResult->getTimeSteps().at(i) < startTimeStep)
            timeStepsToDelete.push_back(i);
    }

    for (int i = 0; i < timeStepsToDelete.size(); i++)
        data.erase(data.begin() + timeStepsToDelete.at(i));

    return data;
}

void NyTestPostProcessingStrategy::evaluate()
{
    if (!isEvaluated) {
        for (int i = 0; i < dataToCalculate.size(); i++) {
            if (dataToCalculate.at(i) == "Vx") {
                ny.at(i) = fftCalculator->calcNy(reduceDataToTimeSteps(simResult->getVx(), startTimeStepCalculation, endTimeStepCalculation), false, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
            }
            if (dataToCalculate.at(i) == "Vy") {
                ny.at(i) = fftCalculator->calcNy(reduceDataToTimeSteps(simResult->getVy(), startTimeStepCalculation, endTimeStepCalculation), false, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
            }
            if (dataToCalculate.at(i) == "Vz") {
                ny.at(i) = fftCalculator->calcNy(reduceDataToTimeSteps(simResult->getVz(), startTimeStepCalculation, endTimeStepCalculation), true, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
            }
            if (dataToCalculate.at(i) == "Press") {
                ny.at(i) = fftCalculator->calcNy(reduceDataToTimeSteps(simResult->getVy(), startTimeStepCalculation, endTimeStepCalculation), false, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
            }
            if (dataToCalculate.at(i) == "Rho") {
                ny.at(i) = fftCalculator->calcNy(reduceDataToTimeSteps(simResult->getVz(), startTimeStepCalculation, endTimeStepCalculation), true, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
            }
        }
        isEvaluated = true; 
    }
}

double NyTestPostProcessingStrategy::getNy(std::string dataToCalculate)
{
    for (int i = 0; i < ny.size(); i++) {
        if (dataToCalculate == this->dataToCalculate.at(i))
            return ny.at(i);
    }
}
    
//! \}

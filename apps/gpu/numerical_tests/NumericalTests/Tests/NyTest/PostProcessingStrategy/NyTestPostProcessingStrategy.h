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
#ifndef NY_TEST_POST_PROCESSING_STRATEGY_H
#define NY_TEST_POST_PROCESSING_STRATEGY_H

#include "Utilities/PostProcessingStrategy/PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class FFTCalculator;
struct NyTestParameterStruct;

class NyTestPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
    static std::shared_ptr<NyTestPostProcessingStrategy> getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<NyTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests);
    void evaluate();

    double getNy(std::string dataToCalculate);

private:
    NyTestPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<NyTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests);
    
    std::vector<std::vector<double> > reduceDataToTimeSteps(std::vector<std::vector<double> > data, unsigned int startTimeStep, unsigned int endTimeStep);

    std::shared_ptr<AnalyticalResults> analyticalResult;
    std::vector<std::string> dataToCalculate;
    std::shared_ptr<FFTCalculator> fftCalculator;
    unsigned int startTimeStepCalculation;
    unsigned int endTimeStepCalculation;
    std::vector<double> ny;
    bool isEvaluated;
};
#endif 
//! \}

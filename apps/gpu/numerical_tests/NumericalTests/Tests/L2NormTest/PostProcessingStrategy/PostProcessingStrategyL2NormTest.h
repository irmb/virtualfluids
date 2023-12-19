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
#ifndef L2NORM_POST_PROCESSING_STRATEGY_H
#define L2NORM_POST_PROCESSING_STRATEGY_H

#include "Utilities/PostProcessingStrategy/PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class L2NormCalculator;
class L2NormCalculatorFactory;
struct L2NormTestParameterStruct;

class L2NormPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
    static std::shared_ptr<L2NormPostProcessingStrategy> getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests);
    void evaluate();

    std::vector<double> getL2Norm(std::string dataToCalc, std::string normalizeData);

    std::string getErrorMessage(std::string aNormalizeData);

private:
    L2NormPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests);
    bool isEvaluated;

    std::shared_ptr<AnalyticalResults> analyticalResult;
    std::vector<std::shared_ptr<L2NormCalculator> > l2Normcalculator;
    
    std::vector<std::string> dataToCalculate;
    std::vector<std::string> normalizeData;
    unsigned int basicTimeStep;
    unsigned int divergentTimeStep;
    std::vector<std::vector<double>> l2NormBasic;
    std::vector<std::vector<double>> l2NormDivergent;
};
#endif
//! \}

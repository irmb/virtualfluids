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
#ifndef L2_NORM_TEST_H
#define L2_NORM_TEST_H

#include "Utilities/Test/TestImp.h"

#include <memory>

class L2NormCalculator;
class AnalyticalResults;
class L2NormPostProcessingStrategy;

struct L2NormTestParameterStruct;

class L2NormTest : public TestImp
{
public:
    static std::shared_ptr<L2NormTest> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate, double maxL2NormDiff, std::string normalizeData);

    void update();
    void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormPostProcessingStrategy> postProStrategy);

    void evaluate();
    std::string getLogFileOutput();
    std::string getErrorLogFileOutput();

private:
    L2NormTest(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate, double maxL2NormDiff, std::string normalizeData);
    std::vector<std::string> buildTestOutput();
    std::vector<std::string> buildBasicTestOutput();
    std::vector<std::string> buildErrorTestOutput();


    unsigned int basicTimeStep, divergentTimeStep;
    double resultBasicTimestep, resultDivergentTimeStep;
    std::string dataToCalculate;
    double diffL2Norm;
    double maxL2NormDiff;
    bool testPassed;

    std::string normalizeData;

    std::vector<std::shared_ptr<L2NormPostProcessingStrategy> > l2NormPostProStrategies;
};
#endif 
//! \}

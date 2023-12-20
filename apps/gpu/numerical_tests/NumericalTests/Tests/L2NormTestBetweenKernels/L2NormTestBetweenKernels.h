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
#ifndef L2_NORM_TEST_BETWEEN_KERNELS_H
#define L2_NORM_TEST_BETWEEN_KERNELS_H

#include "Utilities/Test/TestImp.h"

#include <memory>

class L2NormBetweenKernelPostProcessingStrategy;
class L2NormCalculator;
class L2NormCalculatorFactory;
class AnalyticalResults;
class SimulationResults;

class L2NormTestBetweenKernels : public TestImp
{
public:
    static std::shared_ptr<L2NormTestBetweenKernels> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep, std::string normalizeWith, std::shared_ptr<L2NormCalculatorFactory> factory);

    void update();
    void evaluate();
    std::string getLogFileOutput();
    std::string getErrorLogFileOutput();
    double getBasicL2Result();

    void setBasicSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy);
    void setDivergentKernelSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy);

private:
    L2NormTestBetweenKernels(std::shared_ptr<ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep, std::string normalizeWith, std::shared_ptr<L2NormCalculatorFactory> factory);
    int calcTimeStepInResults(unsigned int timeStep);
    std::vector<std::string> buildTestOutput();
    std::vector<std::string> buildBasicTestOutput();
    std::vector<std::string> buildErrorTestOutput();

    unsigned int timeStep;
    std::string dataToCalculate;
    std::shared_ptr<NumericalTestSimulation> basicSim;
    std::shared_ptr<SimulationInfo> basicSimInfo;
    std::shared_ptr<SimulationResults> basicSimResults;
    std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> basicPostProcessingStrategy;
    double basicL2Result;
    std::shared_ptr<NumericalTestSimulation> divergentSim;
    std::shared_ptr<SimulationInfo> divergentSimInfo;
    std::shared_ptr<SimulationResults> divergentSimResults;
    std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> divergentPostProcessingStrategy;
    double divergentL2Result;
    std::shared_ptr<L2NormCalculator> l2Normcalculator;
    std::string normalizeData;
    double resultL2ToBasicKernel;
};
#endif 
//! \}

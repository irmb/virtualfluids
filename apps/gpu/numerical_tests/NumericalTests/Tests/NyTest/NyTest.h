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
#ifndef NY_TEST_H
#define NY_TEST_H

#include "Utilities/Test/TestImp.h"

#include <memory>
#include <vector>
#include <iostream>

class FFTCalculator;
class NyTestPostProcessingStrategy;
struct NyTestParameterStruct;

class NyTest : public TestImp 
{
public:
    static std::shared_ptr<NyTest> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<NyTestParameterStruct> testPara, std::string dataToCalculate);
    
    void update();
    void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<NyTestPostProcessingStrategy> postProStrategy);
    
    void evaluate();

    std::string getDataToCalculate();
    std::vector<int> getLx();
    std::vector<double> getNy();
    std::vector<double> getNyDiff();
    double getOrderOfAccuracyNyDiff();

private:
    NyTest(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<NyTestParameterStruct> testPara, std::string dataToCalculate);
    double calcOrderOfAccuracy(std::vector<double> data);
    TestStatus checkTestPassed(double orderOfAccuracy);
    bool checkNy(std::vector<double> ny);
    std::vector<double> calcNyDiff(std::vector<double> ny);
    std::vector<std::string> buildTestOutput();
    std::vector<std::string> buildBasicTestOutput();
    std::vector<std::string> buildErrorTestOutput();
    std::vector<std::string> buildSimulationFailedTestOutput();


    unsigned int startStepCalculation, endStepCalculation;
    std::vector<double> lx;
    std::vector<double> ny, nyDiff;
    double orderOfAccuracy;
    double minOrderOfAccuracy;
    double viscosity;
    std::string dataToCalculate;

    std::vector<std::shared_ptr<NyTestPostProcessingStrategy> > postProStrategies;

};
#endif

//! \}

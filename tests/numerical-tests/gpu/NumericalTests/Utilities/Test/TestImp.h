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
#ifndef TEST_IMP_H
#define TEST_IMP_H

#include "Utilities/Test/Test.h"

#include <vector>
#include <memory>

class NumericalTestSimulation;
class SimulationResults;
class SimulationInfo;
class ColorConsoleOutput;
class PostProcessingStrategy;

class TestImp : public Test
{
public:
    void run() override;
    void update() override;
    TestStatus getTestStatus() override;
    void makeConsoleOutput() override;

    void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<PostProcessingStrategy> postProStrategy);
        
protected:
    TestImp(std::shared_ptr<ColorConsoleOutput> colorOutput);

    virtual void evaluate() = 0;
    virtual std::vector<std::string> buildTestOutput() = 0;
    virtual std::vector<std::string> buildBasicTestOutput() = 0;
    virtual std::vector<std::string> buildErrorTestOutput() = 0;
    std::vector<std::string> buildSimulationFailedTestOutput();
    
    bool CheckAllSimulationRun();

    std::vector<std::shared_ptr<NumericalTestSimulation> > simulations;
    std::vector<std::shared_ptr<PostProcessingStrategy> > postProStrategies;
    std::vector<std::shared_ptr<SimulationInfo> > simInfos;
    std::vector<bool> simulationRun;
    std::shared_ptr<ColorConsoleOutput> colorOutput;
    TestStatus testStatus;

private:
    TestImp() {};
};
#endif 
//! \}

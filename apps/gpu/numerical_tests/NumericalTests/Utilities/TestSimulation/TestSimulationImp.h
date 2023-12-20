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
#ifndef TEST_SIMULATION_IMP_H
#define TEST_SIMULATION_IMP_H

#include "TestSimulation.h"
#include "Utilities/NumericalTestSimulation/NumericalTestSimulation.h"

#include <functional>
#include <ctime>
#include <vector>

class ToVectorWriter;
class ColorConsoleOutput;
class SimulationInfo;
class AnalyticalResults;
class AnalyticalResults2DToVTKWriter;
class SimulationResults;
class TimeTracking;

struct TestSimulationDataStruct;

class TestSimulationImp : public TestSimulation, public NumericalTestSimulation
{
public:
    TestSimulationImp(std::function<void()> runSimulation, std::shared_ptr<TestSimulationDataStruct> testSimData,
                      std::shared_ptr<SimulationResults> simResult, std::shared_ptr<TimeTracking> timeTracking,
                      std::shared_ptr<ToVectorWriter> toVectorWriter,
                      std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter,
                      std::shared_ptr<ColorConsoleOutput> colorOutput);
    void run() override;

    std::shared_ptr<SimulationParameter> getSimulationParameter() override;
    std::shared_ptr<SimulationInfo> getSimulationInfo();
    std::shared_ptr<TimeTracking> getTimeTracking() override;

    SimulationStatus getSimulationStatus() override;

    void makeSimulationHeadOutput() override;
    void startPostProcessing() override;

    void setParameter(std::shared_ptr<Parameter> para) override;

    std::shared_ptr<SimulationResults> getSimulationResults();
    std::shared_ptr<AnalyticalResults> getAnalyticalResults();
    void registerSimulationObserver(std::shared_ptr<SimulationObserver> simObserver) override;
    std::vector<std::string> getDataToCalcTests();

private:
    void notifyObserver();

    void writeAnalyticalResultsToVTK();
    void checkSimulationResults();

    std::shared_ptr<SimulationParameter> simPara;
    std::shared_ptr<ToVectorWriter> toVectorWriter;
    std::shared_ptr<InitialCondition> initialCondition;
    std::shared_ptr<SimulationInfo> simInfo;
    std::shared_ptr<SimulationResults> simResult;
    std::shared_ptr<TimeTracking> timeTracking;

    std::shared_ptr<AnalyticalResults> analyticalResult;

    std::shared_ptr<ColorConsoleOutput> colorOutput;
    std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter;
    std::shared_ptr<Parameter> para;
    std::vector<std::shared_ptr<SimulationObserver>> simObserver;

    std::vector<std::string> dataToCalcTests;
    SimulationStatus status;

    std::function<void()> runSimulation;
};
#endif
//! \}

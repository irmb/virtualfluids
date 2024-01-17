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
#include <functional>

#include "TestSimulationImp.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities/DataWriter/AnalyticalResults2DToVTKWriter/AnalyticalResults2DToVTKWriter.h"
#include "Utilities/DataWriter/ToVectorWriter.h"
#include "Utilities/InitialCondition/InitialCondition.h"
#include "Utilities/KernelConfiguration/KernelConfiguration.h"
#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities/SimulationInfo/SimulationInfo.h"
#include "Utilities/Structs/TestSimulationDataStruct.h"
#include "Utilities/Test/SimulationObserver.h"
#include "Utilities/Time/TimeTracking.h"

TestSimulationImp::TestSimulationImp(std::function<void()> runSimulation,
                                     std::shared_ptr<TestSimulationDataStruct> testSimData,
                                     std::shared_ptr<SimulationResults> simResult,
                                     std::shared_ptr<TimeTracking> timeTracking,
                                     std::shared_ptr<ToVectorWriter> toVectorWriter,
                                     std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter,
                                     std::shared_ptr<ColorConsoleOutput> colorOutput)
{
    this->simPara = testSimData->simParameter;
    this->simInfo = testSimData->simInformation;
    this->analyticalResult = testSimData->analyticalResult;
    this->initialCondition = testSimData->initialCondition;

    this->timeTracking = timeTracking;

    this->simResult = simResult;
    this->toVectorWriter = toVectorWriter;

    this->anaResultWriter = anaResultWriter;
    this->colorOutput = colorOutput;

    this->simObserver.resize(0);
    this->dataToCalcTests = simInfo->getDataToCalcTests();
    this->status = initialized;
    this->runSimulation = runSimulation;
}

void TestSimulationImp::run()
{
    makeSimulationHeadOutput();
    timeTracking->setSimulationStartTime();
    runSimulation();
    timeTracking->setSimulationEndTime();
    startPostProcessing();
}

std::shared_ptr<SimulationParameter> TestSimulationImp::getSimulationParameter()
{
    return simPara;
}

std::shared_ptr<AnalyticalResults> TestSimulationImp::getAnalyticalResults()
{
    return analyticalResult;
}

std::shared_ptr<SimulationInfo> TestSimulationImp::getSimulationInfo()
{
    return simInfo;
}

std::shared_ptr<TimeTracking> TestSimulationImp::getTimeTracking()
{
    return timeTracking;
}

SimulationStatus TestSimulationImp::getSimulationStatus()
{
    return status;
}

void TestSimulationImp::registerSimulationObserver(std::shared_ptr<SimulationObserver> simObserver)
{
    this->simObserver.push_back(simObserver);
}

std::vector<std::string> TestSimulationImp::getDataToCalcTests()
{
    return dataToCalcTests;
}

void TestSimulationImp::notifyObserver()
{
    for (uint i = 0; i < simObserver.size(); i++)
        simObserver.at(i)->update();
}

void TestSimulationImp::writeAnalyticalResultsToVTK()
{
    if (!analyticalResult->isCalculated())
        analyticalResult->calc(simResult);

    anaResultWriter->writeAnalyticalResult(para, analyticalResult);
}

void TestSimulationImp::checkSimulationResults()
{
    bool dataOkay = simResult->checkYourData();
    if (!dataOkay)
        status = crashed;
}

void TestSimulationImp::makeSimulationHeadOutput()
{
    colorOutput->makeSimulationHeadOutput(simInfo);
}

void TestSimulationImp::startPostProcessing()
{
    status = executed;

    timeTracking->setResultCheckStartTime();
    checkSimulationResults();
    timeTracking->setResultCheckEndTime();

    timeTracking->setTestStartTime();
    notifyObserver();
    timeTracking->setTestEndTime();

    timeTracking->setAnalyticalResultWriteStartTime();
    writeAnalyticalResultsToVTK();
    timeTracking->setAnalyticalResultWriteEndTime();
}

void TestSimulationImp::setParameter(std::shared_ptr<Parameter> para)
{
    this->para = para;
}

std::shared_ptr<SimulationResults> TestSimulationImp::getSimulationResults()
{
    return simResult;
}
//! \}

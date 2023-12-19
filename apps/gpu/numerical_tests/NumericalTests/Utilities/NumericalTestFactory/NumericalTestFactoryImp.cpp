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
#include "NumericalTestFactoryImp.h"

#include "Utilities/Structs/ConfigDataStruct.h"
#include "Utilities/Structs/LogFileParameterStruct.h"
#include "Utilities/Structs/NumericalTestStruct.h"
#include "Utilities/Structs/SimulationDataStruct.h"
#include "Utilities/Structs/TestStruct.h"
#include "Utilities/Structs/TestSimulationDataStruct.h"
#include "Utilities/VirtualFluidSimulationFactory/VirtualFluidSimulationFactory.h"

#include "Simulations/TaylorGreenVortexUx/AnalyticalResults/AnalyticalResultsTaylorGreenVortexUx.h"
#include "Simulations/TaylorGreenVortexUx/InitialConditions/InitialConditionTaylorGreenVortexUx.h"
#include "Simulations/TaylorGreenVortexUx/LogFileInformation/LogFileInformationTaylorGreenVortexUx.h"
#include "Simulations/TaylorGreenVortexUx/SimulationInfo/SimulationInfoTaylorGreenVortexUx.h"
#include "Simulations/TaylorGreenVortexUx/SimulationParameter/SimulationParameterTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUz/SimulationParameter/SimulationParameterTaylorGreenVortexUz.h"
#include "Simulations/TaylorGreenVortexUz/LogFileInformation/LogFileInformationTaylorGreenVortexUz.h"
#include "Simulations/TaylorGreenVortexUz/SimulationInfo/SimulationInfoTaylorGreenVortexUz.h"
#include "Simulations/TaylorGreenVortexUz/AnalyticalResults/AnalyticalResultsTaylorGreenVortexUz.h"
#include "Simulations/TaylorGreenVortexUz/InitialConditions/InitialConditionTaylorGreenVortexUz.h"

#include "Simulations/ShearWave/SimulationParameter/ShearWaveSimulationParameter.h"
#include "Simulations/ShearWave/LogFileInformation/ShearWaveLogFileInformation.h"
#include "Simulations/ShearWave/SimulationInfo/ShearWaveSimulationInfo.h"
#include "Simulations/ShearWave/AnalyticalResults/ShearWaveAnalyticalResults.h"
#include "Simulations/ShearWave/InitialConditions/InitialConditionShearWave.h"

#include "Tests/NyTest/NyTest.h"
#include "Tests/NyTest/LogFileInformation/NyTestLogFileInformation.h"
#include "Tests/NyTest/PostProcessingStrategy/NyTestPostProcessingStrategy.h"
#include "Tests/PhiTest/PhiTest.h"
#include "Tests/PhiTest/LogFileInformation/PhiTestLogFileInformation.h"
#include "Tests/PhiTest/PostProcessingStrategy/PhiTestPostProcessingStrategy.h"
#include "Tests/L2NormTest/L2NormTest.h"
#include "Tests/L2NormTest/LogFileInformation/L2NormLogFileInformation.h"
#include "Tests/L2NormTest/PostProcessingStrategy/PostProcessingStrategyL2NormTest.h"
#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernels.h"
#include "Tests/L2NormTestBetweenKernels/PostProcessingStrategy/L2NormBetweenKernelPostProcessingStrategy.h"
#include "Tests/L2NormTestBetweenKernels/LogFileInformation/L2NormLogFileInformationBetweenKernels.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutputImp.h"
#include "Utilities/Calculator/L2NormCalculator/L2NormCalculatorFactory/L2NormCalculatorFactoryImp.h"
#include "Utilities/DataWriter/AnalyticalResults2DToVTKWriter/AnalyticalResults2DToVTKWriterImp.h"
#include "Utilities/DataWriter/Y2dSliceToResults/Y2dSliceToResults.h"

#include "Utilities/LogFileInformation/LogFileHead/LogFileHead.h"
#include "Utilities/LogFileInformation/BasicSimulationInfo/BasicSimulationInfo.h"
#include "Utilities/LogFileInformation/BasicTestLogFileInformation/BasicTestLogFileInformation.h"
#include "Utilities/LogFileInformation/LogFileInformationImp.h"
#include "Utilities/LogFileInformation/LogFileTimeInformation/LogFileTimeInformation.h"
#include "Utilities/LogFileWriter/LogFileWriterImp.h"
#include "Utilities/LogFileQueue/LogFileQueueImp.h"

#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities/TestQueue/TestQueueImp.h"
#include "Utilities/TestSimulation/TestSimulationImp.h"
#include "Utilities/Time/TimeImp.h"

#include <algorithm>


std::shared_ptr<NumericalTestFactoryImp> NumericalTestFactoryImp::getNewInstance(std::shared_ptr<ConfigDataStruct> configFileData)
{
    return std::shared_ptr<NumericalTestFactoryImp>(new NumericalTestFactoryImp(configFileData));
}

NumericalTestFactoryImp::NumericalTestFactoryImp(std::shared_ptr<ConfigDataStruct> configFileData)
{
    colorOutput = ColorConsoleOutputImp::getInstance();
    myTestQueue = TestQueueImp::getNewInstance(colorOutput);
    myLogFileWriterQueue = LogFileQueueImp::getNewInstance(configFileData->logFilePath);
    anaResultWriter = AnalyticalResults2DToVTKWriterImp::getInstance(configFileData->writeAnalyticalToVTK);
    l2NormCalculatorFactory = L2NormCalculatorFactoryImp::getInstance();
    l2NormTestsBetweenKernels.resize(0);
    init(configFileData);
}

std::vector<std::shared_ptr<TestSimulation> > NumericalTestFactoryImp::getTestSimulations()
{
    return myTestSimulations;
}

std::shared_ptr<TestQueue> NumericalTestFactoryImp::getTestQueue()
{
    return myTestQueue;
}

std::shared_ptr<LogFileQueue> NumericalTestFactoryImp::getLogFileQueue()
{
    return myLogFileWriterQueue;
}

void NumericalTestFactoryImp::init(std::shared_ptr<ConfigDataStruct> configFileData)
{
    simID = 1;
    numberOfSimulations = configFileData->numberOfSimulations;

    for (size_t i = 0; i < configFileData->kernelsToTest.size(); i++) {
        for (size_t j = 0; j < configFileData->viscosity.size(); j++) {
            for (size_t k = 0; k < configFileData->taylorGreenVortexUxParameter.size(); k++) {
                std::shared_ptr<SimulationDataStruct> simDataStruct = makeTaylorGreenUxSimulationData(configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->taylorGreenVortexUxParameter.at(k), configFileData->taylorGreenVortexUxGridInformation);
                if (simDataStruct->simGroupRun) {
                    std::shared_ptr<NumericalTestStruct> numericalTestStruct = makeNumericalTestStruct(configFileData, simDataStruct, configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->taylorGreenVortexUxParameter.at(k)->basicTimeStepLength);
                    addNumericalTestStruct(numericalTestStruct);
                }
            }

            for (size_t k = 0; k < configFileData->taylorGreenVortexUzParameter.size(); k++) {
                std::shared_ptr<SimulationDataStruct> simDataStruct = makeTaylorGreenUzSimulationData(configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->taylorGreenVortexUzParameter.at(k), configFileData->taylorGreenVortexUzGridInformation);
                if (simDataStruct->simGroupRun) {
                    std::shared_ptr<NumericalTestStruct> numericalTestStruct = makeNumericalTestStruct(configFileData, simDataStruct, configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->taylorGreenVortexUzParameter.at(k)->basicTimeStepLength);
                    addNumericalTestStruct(numericalTestStruct);
                }
            }

            for (size_t k = 0; k < configFileData->shearWaveParameter.size(); k++) {
                std::shared_ptr<SimulationDataStruct> simDataStruct = makeShearWaveSimulationData(configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->shearWaveParameter.at(k), configFileData->shearWaveGridInformation);
                if (simDataStruct->simGroupRun) {
                    std::shared_ptr<NumericalTestStruct> numericalTestStruct = makeNumericalTestStruct(configFileData, simDataStruct, configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->shearWaveParameter.at(k)->basicTimeStepLength);
                    addNumericalTestStruct(numericalTestStruct);
                }

            }
            
        }
    }
}

std::shared_ptr<NumericalTestStruct> NumericalTestFactoryImp::makeNumericalTestStruct(std::shared_ptr<ConfigDataStruct> configFileData, std::shared_ptr<SimulationDataStruct> simDataStruct, std::string kernel, double viscosity, int basicTimeStepLength)
{
    std::shared_ptr<NumericalTestStruct> numTestStruct = std::shared_ptr<NumericalTestStruct>(new NumericalTestStruct);

    std::vector<std::shared_ptr<TestSimulationImp> > testSim = makeTestSimulations(simDataStruct->testSimData, configFileData->vectorWriterInfo, configFileData->ySliceForCalculation);
    numTestStruct->testSimulations = testSim;
    std::shared_ptr<BasicTestLogFileInformation> basicTestLogFileInfo = BasicTestLogFileInformation::getNewInstance();
    std::vector<std::shared_ptr<TestLogFileInformation> > testLogFileInfo;
    
    std::shared_ptr<TestStruct> phiTestStruct = makePhiTestsStructs(configFileData->phiTestParameter, testSim, viscosity);
    initTestStruct(phiTestStruct, numTestStruct, testLogFileInfo, basicTestLogFileInfo);

    std::shared_ptr<TestStruct> nyTestStruct = makeNyTestsStructs(configFileData->nyTestParameter, testSim, viscosity);
    initTestStruct(nyTestStruct, numTestStruct, testLogFileInfo, basicTestLogFileInfo);
        
    std::shared_ptr<TestStruct> l2NormTestSruct = makeL2NormTestsStructs(configFileData->l2NormTestParameter, testSim);
    initTestStruct(l2NormTestSruct, numTestStruct, testLogFileInfo, basicTestLogFileInfo);

    std::shared_ptr<TestStruct> l2NormTestBetweenKernelStruct = makeL2NormTestsBetweenKernelsStructs(configFileData->l2NormTestBetweenKernelsParameter, testSim, kernel);
    initTestStruct(l2NormTestBetweenKernelStruct, numTestStruct, testLogFileInfo, basicTestLogFileInfo);

    std::vector<std::shared_ptr<SimulationInfo> > simInfo;
    for (size_t i = 0; i < simDataStruct->testSimData.size(); i++)
        simInfo.push_back(simDataStruct->testSimData.at(i)->simInformation);

    std::shared_ptr<LogFileWriter> logFileWriter = makeLogFileWriter(testLogFileInfo, simDataStruct->logFileInformation, simInfo, kernel, viscosity, basicTimeStepLength, configFileData->logFilePara, basicTestLogFileInfo);
    numTestStruct->logFileWriter = logFileWriter;

    return numTestStruct;
}

void NumericalTestFactoryImp::addNumericalTestStruct(std::shared_ptr<NumericalTestStruct> numericalTestStruct)
{
    for (size_t i = 0; i < numericalTestStruct->testSimulations.size(); i++)
        myTestSimulations.push_back(numericalTestStruct->testSimulations.at(i));

    for (size_t i = 0; i < numericalTestStruct->tests.size(); i++)
        myTestQueue->addTest(numericalTestStruct->tests.at(i));

    myLogFileWriterQueue->addLogFileWriter(numericalTestStruct->logFileWriter);
}

std::shared_ptr<SimulationDataStruct> NumericalTestFactoryImp::makeTaylorGreenUxSimulationData(std::string kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
    std::shared_ptr<SimulationDataStruct> simDataStruct = std::shared_ptr<SimulationDataStruct>(new SimulationDataStruct);

    if (gridInfoStruct.size() > 0) {
        for (size_t i = 0; i < gridInfoStruct.size(); i++) {
            std::shared_ptr<TestSimulationDataStruct> aTestSimData = std::shared_ptr<TestSimulationDataStruct>(new TestSimulationDataStruct);
            aTestSimData->simParameter = SimulationParameterTaylorGreenUx::getNewInstance(kernel, viscosity, simParaStruct, gridInfoStruct.at(i));
            aTestSimData->initialCondition = InitialConditionTaylorGreenUx::getNewInstance(simParaStruct, gridInfoStruct.at(i));
            aTestSimData->simInformation = SimulationInfoTaylorGreenUx::getNewInstance(simID, kernel, viscosity, simParaStruct, gridInfoStruct.at(i), numberOfSimulations);
            simID++;
            aTestSimData->analyticalResult = AnalyticalResultsTaylorGreenUx::getNewInstance(viscosity, simParaStruct);
            simDataStruct->testSimData.push_back(aTestSimData);
        }
        simDataStruct->logFileInformation = LogFileInformationTaylorGreenUx::getNewInstance(simParaStruct, gridInfoStruct);
        simDataStruct->simGroupRun = true;
    }
    else {
        simDataStruct->simGroupRun = false;
    }
    return simDataStruct;
}

std::shared_ptr<SimulationDataStruct> NumericalTestFactoryImp::makeTaylorGreenUzSimulationData(std::string kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
    std::shared_ptr<SimulationDataStruct> simDataStruct = std::shared_ptr<SimulationDataStruct>(new SimulationDataStruct);
    if (gridInfoStruct.size() > 0) {
        for (size_t i = 0; i < gridInfoStruct.size(); i++) {
            std::shared_ptr<TestSimulationDataStruct> aTestSimData = std::shared_ptr<TestSimulationDataStruct>(new TestSimulationDataStruct);
            aTestSimData->simParameter = SimulationParameterTaylorGreenUz::getNewInstance(kernel, viscosity, simParaStruct, gridInfoStruct.at(i));
            aTestSimData->initialCondition = InitialConditionTaylorGreenUz::getNewInstance(simParaStruct, gridInfoStruct.at(i));
            aTestSimData->simInformation = SimulationInfoTaylorGreenUz::getNewInstance(simID, kernel, viscosity, simParaStruct, gridInfoStruct.at(i), numberOfSimulations);
            simID++;
            aTestSimData->analyticalResult = AnalyticalResultsTaylorGreenUz::getNewInstance(viscosity, simParaStruct);
            simDataStruct->testSimData.push_back(aTestSimData);
        }
        simDataStruct->logFileInformation = LogFileInformationTaylorGreenUz::getNewInstance(simParaStruct, gridInfoStruct);
        simDataStruct->simGroupRun = true;
    }
    else {
        simDataStruct->simGroupRun = false;
    }
    return simDataStruct;
}

std::shared_ptr<SimulationDataStruct> NumericalTestFactoryImp::makeShearWaveSimulationData(std::string kernel, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
    std::shared_ptr<SimulationDataStruct> simDataStruct = std::shared_ptr<SimulationDataStruct>(new SimulationDataStruct);
    if (gridInfoStruct.size() > 0) {
        for (size_t i = 0; i < gridInfoStruct.size(); i++) {
            std::shared_ptr<TestSimulationDataStruct> aTestSimData = std::shared_ptr<TestSimulationDataStruct>(new TestSimulationDataStruct);
            aTestSimData->simParameter = ShearWaveSimulationParameter::getNewInstance(kernel, viscosity, simParaStruct, gridInfoStruct.at(i));
            aTestSimData->initialCondition = InitialConditionShearWave::getNewInstance(simParaStruct, gridInfoStruct.at(i));
            aTestSimData->simInformation = ShearWaveSimulationInfo::getNewInstance(simID, kernel, viscosity, simParaStruct, gridInfoStruct.at(i), numberOfSimulations);
            simID++;
            aTestSimData->analyticalResult = ShearWaveAnalyticalResults::getNewInstance(viscosity, simParaStruct);
            simDataStruct->testSimData.push_back(aTestSimData);
        }
        simDataStruct->logFileInformation = ShearWaveInformation::getNewInstance(simParaStruct, gridInfoStruct);
        simDataStruct->simGroupRun = true;
    }        
    else {
        simDataStruct->simGroupRun = false;
    }
    return simDataStruct;
}

std::vector<std::shared_ptr<TestSimulationImp> > NumericalTestFactoryImp::makeTestSimulations(std::vector<std::shared_ptr<TestSimulationDataStruct> > testSimDataStruct, std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int ySliceForCalculation)
{
    std::vector<std::shared_ptr<TestSimulationImp> > testSimulations;
    for (size_t i = 0; i < testSimDataStruct.size(); i++) {
        std::shared_ptr<TimeImp> time = TimeImp::getNewInstance();
        testSimDataStruct.at(i)->simInformation->setTimeInfo(time);
        std::shared_ptr<SimulationResults> simResult = SimulationResults::getNewInstance(testSimDataStruct.at(i)->simParameter);
        std::shared_ptr<ToVectorWriter> toVectorWriter = Y2dSliceToResults::getNewInstance(vectorWriterInfo, testSimDataStruct.at(i)->simParameter->getTimeStepLength(), simResult, ySliceForCalculation);
        

        auto currentTestSimData = testSimDataStruct.at(i);
        auto para = vf::gpu::tests::makeParameter(currentTestSimData->simParameter);
        currentTestSimData->initialCondition->setParameter(para);
        auto vfsim = vf::gpu::tests::makeVirtualFluidSimulation(para, currentTestSimData->initialCondition, toVectorWriter);

        auto testSim = std::make_shared<TestSimulationImp>(vfsim, currentTestSimData, simResult, time, toVectorWriter, anaResultWriter, colorOutput);
        testSim->setParameter(para);

        testSimulations.push_back(testSim);
    }

    return testSimulations;
}

std::shared_ptr<TestStruct> NumericalTestFactoryImp::makePhiTestsStructs(std::shared_ptr<PhiTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp>> testSimumlations, double viscosity)
{
    std::shared_ptr<TestStruct> testStruct = std::shared_ptr<TestStruct>(new TestStruct);

    if (testParameter->basicTestParameter->runTest && testSimumlations.size() > 1) {
        std::shared_ptr<PhiTestLogFileInformation> testLogFileInfo = PhiTestLogFileInformation::getNewInstance(testParameter);
        
        std::vector<std::shared_ptr<PhiTestPostProcessingStrategy> > postProcessingStrategies;
        for (size_t i = 0; i < testSimumlations.size(); i++)
            postProcessingStrategies.push_back(PhiTestPostProcessingStrategy::getNewInstance(testSimumlations.at(i)->getSimulationResults(), testSimumlations.at(i)->getAnalyticalResults(), testParameter, testSimumlations.at(i)->getDataToCalcTests()));

        for (size_t i = 0; i < testSimumlations.at(0)->getDataToCalcTests().size(); i++) {
            std::vector<std::shared_ptr<PhiTest> > phiTests = makePhiTests(testParameter, testSimumlations, postProcessingStrategies, viscosity, testSimumlations.at(0)->getDataToCalcTests().at(i));
            testLogFileInfo->addTestGroup(phiTests);
            for (size_t j = 0; j < phiTests.size(); j++)
                testStruct->tests.push_back(phiTests.at(j));
        }
        testStruct->logFileInfo = testLogFileInfo;
        testStruct->testName = "PhiTest";
    }

    return testStruct;
}

std::vector<std::shared_ptr<PhiTest>> NumericalTestFactoryImp::makePhiTests(std::shared_ptr<PhiTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp>> testSim, std::vector<std::shared_ptr<PhiTestPostProcessingStrategy>> phiPostProStrategy, double viscosity, std::string dataToCalculate)
{
    std::vector<std::shared_ptr<PhiTest> > phiTests;
    for (size_t i = 1; i < testSim.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            std::shared_ptr<PhiTest> test = PhiTest::getNewInstance(colorOutput, viscosity, testParameter, dataToCalculate);
            test->addSimulation(testSim.at(j), testSim.at(j)->getSimulationInfo(), phiPostProStrategy.at(j));
            test->addSimulation(testSim.at(i), testSim.at(i)->getSimulationInfo(), phiPostProStrategy.at(i));

            testSim.at(j)->registerSimulationObserver(test);
            testSim.at(i)->registerSimulationObserver(test);

            phiTests.push_back(test);
        }
    }
    return phiTests;
}

std::shared_ptr<TestStruct> NumericalTestFactoryImp::makeNyTestsStructs(std::shared_ptr<NyTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp>> testSimumlations, double viscosity)
{
    std::shared_ptr<TestStruct> testStruct = std::shared_ptr<TestStruct>(new TestStruct);

    if (testParameter->basicTestParameter->runTest && testSimumlations.size() > 1) {
        std::shared_ptr<NyTestLogFileInformation> testLogFileInfo = NyTestLogFileInformation::getNewInstance(testParameter);

        std::vector<std::shared_ptr<NyTestPostProcessingStrategy> > postProcessingStrategies;
        for (size_t i = 0; i < testSimumlations.size(); i++)
            postProcessingStrategies.push_back(NyTestPostProcessingStrategy::getNewInstance(testSimumlations.at(i)->getSimulationResults(), testSimumlations.at(i)->getAnalyticalResults(), testParameter, testSimumlations.at(i)->getDataToCalcTests()));

        for (size_t i = 0; i < testSimumlations.at(0)->getDataToCalcTests().size(); i++) {
            std::vector<std::shared_ptr<NyTest> > nyTests = makeNyTests(testParameter, testSimumlations, postProcessingStrategies, viscosity, testSimumlations.at(0)->getDataToCalcTests().at(i));
            testLogFileInfo->addTestGroup(nyTests);
            for (size_t j = 0; j < nyTests.size(); j++)
                testStruct->tests.push_back(nyTests.at(j));
        }
        testStruct->logFileInfo = testLogFileInfo;
        testStruct->testName = "NyTest";
    }

    return testStruct;
}

std::vector<std::shared_ptr<NyTest>> NumericalTestFactoryImp::makeNyTests(std::shared_ptr<NyTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp>> testSim, std::vector<std::shared_ptr<NyTestPostProcessingStrategy>> nuPostProStrategy, double viscosity, std::string dataToCalculate)
{
    std::vector<std::shared_ptr<NyTest> > nyTests;
    for (size_t i = 1; i < testSim.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            std::shared_ptr<NyTest> test = NyTest::getNewInstance(colorOutput, viscosity, testParameter, dataToCalculate);
            test->addSimulation(testSim.at(j), testSim.at(j)->getSimulationInfo(), nuPostProStrategy.at(j));
            test->addSimulation(testSim.at(i), testSim.at(i)->getSimulationInfo(), nuPostProStrategy.at(i));

            testSim.at(j)->registerSimulationObserver(test);
            testSim.at(i)->registerSimulationObserver(test);

            nyTests.push_back(test);
        }
    }
    return nyTests;
}

std::shared_ptr<TestStruct> NumericalTestFactoryImp::makeL2NormTestsStructs(std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSimumlations)
{
    std::shared_ptr<TestStruct> testStruct = std::shared_ptr<TestStruct> (new TestStruct);

    if (testParameter->basicTestParameter->runTest) {
        std::vector<std::shared_ptr<L2NormPostProcessingStrategy> >  postProcessingStrategies;
        for (size_t i = 0; i < testSimumlations.size(); i++)
            postProcessingStrategies.push_back(L2NormPostProcessingStrategy::getNewInstance(testSimumlations.at(i)->getSimulationResults(), testSimumlations.at(i)->getAnalyticalResults(), testParameter, l2NormCalculatorFactory, testSimumlations.at(i)->getDataToCalcTests()));

        std::vector<std::shared_ptr<L2NormTest> > tests = makeL2NormTests(testSimumlations, postProcessingStrategies, testParameter);
        std::shared_ptr<L2NormInformation> testLogFileInfo = L2NormInformation::getNewInstance(tests, testParameter, testSimumlations.at(0)->getDataToCalcTests());

        for(size_t i = 0; i < tests.size(); i++)
            testStruct->tests.push_back(tests.at(i));
        testStruct->logFileInfo = testLogFileInfo;
        testStruct->testName = "L2NormTest";
    }
    return testStruct;
}

std::vector<std::shared_ptr<L2NormTest> > NumericalTestFactoryImp::makeL2NormTests(std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormPostProcessingStrategy> > postProStrategy, std::shared_ptr<L2NormTestParameterStruct> testParameter)
{
    std::vector<std::shared_ptr<L2NormTest> > l2Tests;
    for (size_t k = 0; k < testParameter->normalizeData.size(); k++) {
        for (size_t i = 0; i < testSim.size(); i++) {
            for (size_t j = 0; j < testSim.at(i)->getDataToCalcTests().size(); j++) {
                std::shared_ptr<L2NormTest> test = L2NormTest::getNewInstance(colorOutput, testParameter, testSim.at(i)->getDataToCalcTests().at(j), testParameter->maxDiff.at(k), testParameter->normalizeData.at(k));
                test->addSimulation(testSim.at(i), testSim.at(i)->getSimulationInfo(), postProStrategy.at(i));
                testSim.at(i)->registerSimulationObserver(test);
                l2Tests.push_back(test);
            }
        }
    }
    return l2Tests;
}

std::shared_ptr<TestStruct> NumericalTestFactoryImp::makeL2NormTestsBetweenKernelsStructs(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::string kernelName)
{
    std::shared_ptr<TestStruct> testStruct = std::shared_ptr<TestStruct>(new TestStruct);
    testStruct->testName = "L2NormTestBetweenKernel";

    if (testPara->basicTestParameter->runTest) {

        std::vector<std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> > postProcessingStrategies;
        for (size_t i = 0; i < testSim.size(); i++)
            postProcessingStrategies.push_back(L2NormBetweenKernelPostProcessingStrategy::getNewInstance(testSim.at(i)->getSimulationResults(), testSim.at(i)->getAnalyticalResults(), testPara, l2NormCalculatorFactory, testSim.at(i)->getDataToCalcTests()));

        if (kernelName == testPara->basicKernel) {
            std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > > tests = makeL2NormTestsBetweenKernels(testPara, testSim, postProcessingStrategies);
            
            if (l2NormTestsBetweenKernels.size() == 0) {
                l2NormTestsBetweenKernels = tests;
            }
            else {
                for (size_t i = 0; i < tests.size(); i++)
                    for (size_t j = 0; j < tests.at(i).size(); j++)
                        l2NormTestsBetweenKernels.at(i).push_back(tests.at(i).at(j));
            }

        }else{
            std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests = linkL2NormTestsBetweenKernels(testPara, testSim, postProcessingStrategies);
            for (size_t i = 0; i < tests.size(); i++)
                testStruct->tests.push_back(tests.at(i));
            testStruct->logFileInfo = L2NormBetweenKernelsInformation::getNewInstance(tests, testPara, testSim.at(0)->getDataToCalcTests());
        }
    }
    return testStruct;    
}

std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > >  NumericalTestFactoryImp::makeL2NormTestsBetweenKernels(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> > postProcessingStrategies)
{
    std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > > testsForAllKernels;

    std::vector<std::shared_ptr<L2NormTestBetweenKernels> > testForOneKernel;

    for (size_t l = 0; l < testPara->kernelsToTest.size() - 1; l++) {
        for (size_t k = 0; k < testSim.size(); k++) {
            for(size_t j = 0; j < testSim.at(k)->getDataToCalcTests().size(); j++){
                for (size_t m = 0; m < testPara->normalizeData.size(); m++) {
                    for (size_t i = 0; i < testPara->timeSteps.size(); i++) {
                        std::shared_ptr<L2NormTestBetweenKernels> aTest = L2NormTestBetweenKernels::getNewInstance(colorOutput, testSim.at(k)->getDataToCalcTests().at(j), testPara->timeSteps.at(i), testPara->normalizeData.at(m), l2NormCalculatorFactory);
                        aTest->setBasicSimulation(testSim.at(k), testSim.at(k)->getSimulationInfo(), postProcessingStrategies.at(k));
                        testSim.at(k)->registerSimulationObserver(aTest);
                        testForOneKernel.push_back(aTest);
                    }
                }
            }
        }
        testsForAllKernels.push_back(testForOneKernel);
        testForOneKernel.resize(0);
    }
        
    return testsForAllKernels;
}

std::vector<std::shared_ptr<L2NormTestBetweenKernels> > NumericalTestFactoryImp::linkL2NormTestsBetweenKernels(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> > postProcessingStrategies)
{
    std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests;

    if (testSim.size() > 0)
        if (l2NormTestsBetweenKernels.at(0).size() == 0)
            l2NormTestsBetweenKernels.erase(l2NormTestsBetweenKernels.begin());

    for (size_t k = 0; k < testSim.size(); k++) {
        for (size_t j = 0; j < testSim.at(k)->getDataToCalcTests().size(); j++) {
            for (size_t m = 0; m < testPara->normalizeData.size(); m++) {
                for (size_t i = 0; i < testPara->timeSteps.size(); i++) {
                    std::shared_ptr<L2NormTestBetweenKernels> aTest = l2NormTestsBetweenKernels.at(0).at(0);
                    l2NormTestsBetweenKernels.at(0).erase(l2NormTestsBetweenKernels.at(0).begin());
                    aTest->setDivergentKernelSimulation(testSim.at(k), testSim.at(k)->getSimulationInfo(), postProcessingStrategies.at(k));
                    testSim.at(k)->registerSimulationObserver(aTest);
                    tests.push_back(aTest);
                }
            }
        }
    }
    return tests;
}

void NumericalTestFactoryImp::initTestStruct(std::shared_ptr<TestStruct> testStruct, std::shared_ptr<NumericalTestStruct> numericalTestStruct, std::vector<std::shared_ptr<TestLogFileInformation> > &testLogFileInfo, std::shared_ptr<BasicTestLogFileInformation> basicTestLogFileInfo)
{
    for (size_t i = 0; i < testStruct->tests.size(); i++)
        numericalTestStruct->tests.push_back(testStruct->tests.at(i));
    if (testStruct->tests.size() > 0) {
        testLogFileInfo.push_back(testStruct->logFileInfo);
        basicTestLogFileInfo->addTest(testStruct->testName, true);
    }
    else {
        basicTestLogFileInfo->addTest(testStruct->testName, false);
    }
}

std::shared_ptr<LogFileWriter> NumericalTestFactoryImp::makeLogFileWriter(std::vector<std::shared_ptr<TestLogFileInformation> > testLogFiles, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::vector<std::shared_ptr<SimulationInfo> > simInfo, std::string kernel, double viscosity, int basicTimeStepLength, std::shared_ptr<LogFileParameterStruct> logFilePara, std::shared_ptr<BasicTestLogFileInformation> basicTestLogFileInfo)
{
    std::shared_ptr<LogFileHead> logFileHead = LogFileHead::getNewInstance(logFilePara->devices);
    std::shared_ptr<BasicSimulationInfo> basicSimInfo = BasicSimulationInfo::getNewInstance(logFilePara->numberOfTimeSteps, viscosity, basicTimeStepLength, kernel);

    std::shared_ptr<LogFileTimeInformation> logFileTimeInfo = LogFileTimeInformation::getNewInstance(simInfo, logFilePara->writeAnalyticalToVTK);

    std::shared_ptr<LogFileWriterImp> logFileWriter = LogFileWriterImp::getNewInstance(logFileHead, basicSimInfo, basicTestLogFileInfo, testLogFiles, logFileTimeInfo, simLogInfo, kernel, viscosity);

    return logFileWriter;
}
//! \}

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
#ifndef NUMERICAL_TEST_FACTORY_IMP_H
#define NUMERICAL_TEST_FACTORY_IMP_H

#include "NumericalTestFactory.h"

struct ConfigDataStruct;
struct GridInformationStruct;
struct L2NormTestParameterStruct;
struct L2NormTestBetweenKernelsParameterStruct;
struct LogFileParameterStruct;
struct NumericalTestStruct;
struct PhiTestParameterStruct;
struct NyTestParameterStruct;
struct ShearWaveParameterStruct;
struct SimulationDataStruct;
struct TestStruct;
struct TaylorGreenVortexUxParameterStruct;
struct TaylorGreenVortexUzParameterStruct;
struct TestSimulationDataStruct;
struct VectorWriterInformationStruct;


class AnalyticalResults2DToVTKWriter;
class BasicTestLogFileInformation;
class ColorConsoleOutput;
class L2NormCalculatorFactory;
class L2NormTest;
class L2NormPostProcessingStrategy;
class L2NormBetweenKernelPostProcessingStrategy;
class L2NormTestBetweenKernels;
class LogFileTimeInformation;
class LogFileQueueImp;
class LogFileWriter;
class PhiTest;
class PhiTestPostProcessingStrategy;
class NyTest;
class NyTestPostProcessingStrategy;
class SimulationInfo;
class SimulationLogFileInformation;
class SimulationParameter;
class SimulationResults;
class TestQueueImp;
class TestSimulation;
class TestSimulationImp;
class TestLogFileInformation;

class NumericalTestFactoryImp : public NumericalTestFactory
{
public:
    static std::shared_ptr<NumericalTestFactoryImp> getNewInstance(std::shared_ptr<ConfigDataStruct> configFileData);

    std::vector<std::shared_ptr<TestSimulation> > getTestSimulations();
    std::shared_ptr<TestQueue> getTestQueue();
    std::shared_ptr<LogFileQueue> getLogFileQueue();

private:
    NumericalTestFactoryImp() {};
    NumericalTestFactoryImp(std::shared_ptr<ConfigDataStruct> configFileData);

    void init(std::shared_ptr<ConfigDataStruct> configFileData);

    std::shared_ptr<NumericalTestStruct> makeNumericalTestStruct(std::shared_ptr<ConfigDataStruct> configFileData, std::shared_ptr<SimulationDataStruct> simDataStruct, std::string kernel, double viscosity, int basicTimeStepLength);
    void addNumericalTestStruct(std::shared_ptr<NumericalTestStruct> numericalTestStruct);

    std::shared_ptr<SimulationDataStruct> makeTaylorGreenUxSimulationData(std::string kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);
    std::shared_ptr<SimulationDataStruct> makeTaylorGreenUzSimulationData(std::string kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);
    std::shared_ptr<SimulationDataStruct> makeShearWaveSimulationData(std::string kernel, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);

    std::vector<std::shared_ptr<TestSimulationImp> > makeTestSimulations(std::vector<std::shared_ptr<TestSimulationDataStruct> > testSimDataStruct, std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int ySliceForCalculation);

    std::shared_ptr<TestStruct> makePhiTestsStructs(std::shared_ptr<PhiTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSimumlations, double viscosity);
    std::shared_ptr<TestStruct> makeNyTestsStructs(std::shared_ptr<NyTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSimumlations, double viscosity);
    std::shared_ptr<TestStruct> makeL2NormTestsStructs(std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSimumlations);
    std::shared_ptr<TestStruct> makeL2NormTestsBetweenKernelsStructs(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::string kernel);
    
    std::vector<std::shared_ptr<PhiTest> > makePhiTests(std::shared_ptr<PhiTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<PhiTestPostProcessingStrategy> > phiAndNuPostProStrategy, double viscosity, std::string dataToCalculate);
    std::vector<std::shared_ptr<NyTest> > makeNyTests(std::shared_ptr<NyTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<NyTestPostProcessingStrategy> > phiAndNuPostProStrategy, double viscosity, std::string dataToCalculate);
    std::vector<std::shared_ptr<L2NormTest> > makeL2NormTests(std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormPostProcessingStrategy> > postProStrategy, std::shared_ptr<L2NormTestParameterStruct> testParameter);
    std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > > makeL2NormTestsBetweenKernels(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> > postProcessingStrategies);
    std::vector<std::shared_ptr<L2NormTestBetweenKernels> > linkL2NormTestsBetweenKernels(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> > postProcessingStrategies);

    std::shared_ptr<LogFileWriter> makeLogFileWriter(std::vector<std::shared_ptr<TestLogFileInformation> > testLogFiles, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::vector<std::shared_ptr<SimulationInfo> > simInfo, std::string kernel, double viscosity, int basicTimeStepLength, std::shared_ptr<LogFileParameterStruct> logFilePara, std::shared_ptr<BasicTestLogFileInformation> basicTestLogFileInfo);

    void initTestStruct(std::shared_ptr<TestStruct> testStruct, std::shared_ptr<NumericalTestStruct> numericalTestStruct, std::vector<std::shared_ptr<TestLogFileInformation> > &testLogFileInfo, std::shared_ptr<BasicTestLogFileInformation> basicTestLogFileInfo);

    std::vector<std::shared_ptr<TestSimulation> > myTestSimulations;
    std::shared_ptr<TestQueueImp> myTestQueue;
    std::shared_ptr<LogFileQueueImp> myLogFileWriterQueue;
    std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > > l2NormTestsBetweenKernels;
    std::shared_ptr<ColorConsoleOutput> colorOutput;
    std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter;
    std::shared_ptr<L2NormCalculatorFactory> l2NormCalculatorFactory;

    int simID;
    int numberOfSimulations;
    // int simPerKernel, numberOfTestGroupsBetweenKernels, numberOfTestsForOneSimulation, numberOfTestsBetweenKernels;
    // int posBasicSimulationForL2Test, posDivergentSimulationForL2Test;
};
#endif

//! \}

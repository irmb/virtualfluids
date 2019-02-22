#include "NumericalTestFactoryImp.h"

#include "Utilities/Structs/ConfigDataStruct.h"
#include "Utilities/Structs/LogFileParameterStruct.h"
#include "Utilities/Structs/NumericalTestStruct.h"
#include "Utilities/Structs/SimulationDataStruct.h"
#include "Utilities/Structs/TestSimulationDataStruct.h"

#include "Simulations/TaylorGreenVortexUx/AnalyticalResults/AnalyticalResultsTaylorGreenVortexUx.h"
#include "Simulations/TaylorGreenVortexUx/InitialConditions/InitialConditionTaylorGreenVortexUx.h"
#include "Simulations/TaylorGreenVortexUx/LogFileInformation/LogFileInformationTaylorGreenVortexUx.h"
#include "Simulations/TaylorGreenVortexUx/SimulationInfo/SimulationInfoTaylorGreenVortexUx.h"
#include "Simulations/TaylorGreenVortexUx/SimulationParameter/SimulationParameterTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUz/SimulationParameter/SimulationParameterTaylorGreenVortexUz.h"
#include "Simulations/TaylorGreenVortexUz/LogFileInformation/LogFileInformationTaylorGreenvortexUz.h"
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
#include "Tests/NyTest/NyTestStruct.h"
#include "Tests/PhiTest/PhiTest.h"
#include "Tests/PhiTest/LogFileInformation/PhiTestLogFileInformation.h"
#include "Tests/PhiTest/PostProcessingStrategy/PhiTestPostProcessingStrategy.h"
#include "Tests/PhiTest/PhiTestStruct.h"
#include "Tests/L2NormTest/L2NormTest.h"
#include "Tests/L2NormTest/LogFileInformation/L2NormLogFileInformation.h"
#include "Tests/L2NormTest/PostProcessingStrategy/PostProcessingStrategyL2NormTest.h"
#include "Tests/L2NormTest/L2NormTestStruct.h"
#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernels.h"
#include "Tests/L2NormTestBetweenKernels/PostProcessingStrategy/L2NormBetweenKernelPostProcessingStrategy.h"
#include "Tests/L2NormTestBetweenKernels/LogFileInformation/L2NormLogFileInformationBetweenKernels.h"
#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernelsStruct.h"

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

	for (int i = 0; i < configFileData->kernelsToTest.size(); i++) {
		for (int j = 0; j < configFileData->viscosity.size(); j++) {
			for (int k = 0; k < configFileData->taylorGreenVortexUxParameter.size(); k++) {
				std::shared_ptr<SimulationDataStruct> simDataStruct = makeTaylorGreenUxSimulationData(configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->taylorGreenVortexUxParameter.at(k), configFileData->taylorGreenVortexUxGridInformation);
				if (simDataStruct->simGroupRun) {
					std::shared_ptr<NumericalTestStruct> numericalTestStruct = makeNumericalTestStruct(configFileData, simDataStruct, configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->taylorGreenVortexUxParameter.at(k)->basicTimeStepLength);
					addNumericalTestStruct(numericalTestStruct);
				}
			}

			for (int k = 0; k < configFileData->taylorGreenVortexUzParameter.size(); k++) {
				std::shared_ptr<SimulationDataStruct> simDataStruct = makeTaylorGreenUzSimulationData(configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->taylorGreenVortexUzParameter.at(k), configFileData->taylorGreenVortexUzGridInformation);
				if (simDataStruct->simGroupRun) {
					std::shared_ptr<NumericalTestStruct> numericalTestStruct = makeNumericalTestStruct(configFileData, simDataStruct, configFileData->kernelsToTest.at(i), configFileData->viscosity.at(j), configFileData->taylorGreenVortexUzParameter.at(k)->basicTimeStepLength);
					addNumericalTestStruct(numericalTestStruct);
				}
			}

			for (int k = 0; k < configFileData->shearWaveParameter.size(); k++) {
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
	
	std::shared_ptr<PhiTestStruct> phiTestStruct = makePhiTestsStructs(configFileData->phiTestParameter, testSim, viscosity);
	for (int i = 0; i < phiTestStruct->tests.size(); i++)
		numTestStruct->tests.push_back(phiTestStruct->tests.at(i));
	if (phiTestStruct->tests.size() > 0) {
		testLogFileInfo.push_back(phiTestStruct->logFileInfo);
		basicTestLogFileInfo->addTest("PhiTest", true);
	}
	else {
		basicTestLogFileInfo->addTest("PhiTest", false);
	}

	std::shared_ptr<NyTestStruct> nyTestStruct = makeNyTestsStructs(configFileData->nyTestParameter, testSim, viscosity);
	for (int i = 0; i < nyTestStruct->tests.size(); i++)
		numTestStruct->tests.push_back(nyTestStruct->tests.at(i));
	if (nyTestStruct->tests.size() > 0) {
		testLogFileInfo.push_back(nyTestStruct->logFileInfo);
		basicTestLogFileInfo->addTest("NyTest", true);
	}
	else {
		basicTestLogFileInfo->addTest("NyTest", false);
	}
		

	std::shared_ptr<L2NormTestStruct> l2NormTestSruct = makeL2NormTestsStructs(configFileData->l2NormTestParameter, testSim);
	for (int i = 0; i < l2NormTestSruct->tests.size(); i++)
		numTestStruct->tests.push_back(l2NormTestSruct->tests.at(i));
	if (l2NormTestSruct->tests.size() > 0) {
		testLogFileInfo.push_back(l2NormTestSruct->logFileInfo);
		basicTestLogFileInfo->addTest("L2NormTest", true);
	}
	else {
		basicTestLogFileInfo->addTest("L2NormTest", false);
	}

	std::shared_ptr<L2NormTestBetweenKernelsStruct> l2NormTestBetweenKernelStruct = makeL2NormTestsBetweenKernelsStructs(configFileData->l2NormTestBetweenKernelsParameter, testSim, kernel);
	for (int i = 0; i < l2NormTestBetweenKernelStruct->tests.size(); i++)
		numTestStruct->tests.push_back(l2NormTestBetweenKernelStruct->tests.at(i));
	if (l2NormTestBetweenKernelStruct->tests.size() > 0) {
		testLogFileInfo.push_back(l2NormTestBetweenKernelStruct->logFileInfo);
		basicTestLogFileInfo->addTest("L2NormTestBetweenKernel", true);
	}
	else {
		basicTestLogFileInfo->addTest("L2NormTestBetweenKernel", false);
	}

	std::vector<std::shared_ptr<SimulationInfo> > simInfo;
	for (int i = 0; i < simDataStruct->testSimData.size(); i++)
		simInfo.push_back(simDataStruct->testSimData.at(i)->simInformation);

	std::shared_ptr<LogFileWriter> logFileWriter = makeLogFileWriter(testLogFileInfo, simDataStruct->logFileInformation, simInfo, kernel, viscosity, basicTimeStepLength, configFileData->logFilePara, basicTestLogFileInfo);
	numTestStruct->logFileWriter = logFileWriter;

	return numTestStruct;
}

void NumericalTestFactoryImp::addNumericalTestStruct(std::shared_ptr<NumericalTestStruct> numericalTestStruct)
{
	for (int i = 0; i < numericalTestStruct->testSimulations.size(); i++)
		myTestSimulations.push_back(numericalTestStruct->testSimulations.at(i));

	for (int i = 0; i < numericalTestStruct->tests.size(); i++)
		myTestQueue->addTest(numericalTestStruct->tests.at(i));

	myLogFileWriterQueue->addLogFileWriter(numericalTestStruct->logFileWriter);
}

std::shared_ptr<SimulationDataStruct> NumericalTestFactoryImp::makeTaylorGreenUxSimulationData(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
	std::shared_ptr<SimulationDataStruct> simDataStruct = std::shared_ptr<SimulationDataStruct>(new SimulationDataStruct);

	if (gridInfoStruct.size() > 0) {
		for (int i = 0; i < gridInfoStruct.size(); i++) {
			std::shared_ptr<TestSimulationDataStruct> aTestSimData = std::shared_ptr<TestSimulationDataStruct>(new TestSimulationDataStruct);
			aTestSimData->simParameter = SimulationParameterTaylorGreenUx::getNewInstance(kernelName, viscosity, simParaStruct, gridInfoStruct.at(i));
			aTestSimData->initialCondition = InitialConditionTaylorGreenUx::getNewInstance(simParaStruct, gridInfoStruct.at(i));
			aTestSimData->simInformation = SimulationInfoTaylorGreenUx::getNewInstance(simID, kernelName, viscosity, simParaStruct, gridInfoStruct.at(i), numberOfSimulations);
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

std::shared_ptr<SimulationDataStruct> NumericalTestFactoryImp::makeTaylorGreenUzSimulationData(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
	std::shared_ptr<SimulationDataStruct> simDataStruct = std::shared_ptr<SimulationDataStruct>(new SimulationDataStruct);
	if (gridInfoStruct.size() > 0) {
		for (int i = 0; i < gridInfoStruct.size(); i++) {
			std::shared_ptr<TestSimulationDataStruct> aTestSimData = std::shared_ptr<TestSimulationDataStruct>(new TestSimulationDataStruct);
			aTestSimData->simParameter = SimulationParameterTaylorGreenUz::getNewInstance(kernelName, viscosity, simParaStruct, gridInfoStruct.at(i));
			aTestSimData->initialCondition = InitialConditionTaylorGreenUz::getNewInstance(simParaStruct, gridInfoStruct.at(i));
			aTestSimData->simInformation = SimulationInfoTaylorGreenUz::getNewInstance(simID, kernelName, viscosity, simParaStruct, gridInfoStruct.at(i), numberOfSimulations);
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

std::shared_ptr<SimulationDataStruct> NumericalTestFactoryImp::makeShearWaveSimulationData(std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
	std::shared_ptr<SimulationDataStruct> simDataStruct = std::shared_ptr<SimulationDataStruct>(new SimulationDataStruct);
	if (gridInfoStruct.size() > 0) {
		for (int i = 0; i < gridInfoStruct.size(); i++) {
			std::shared_ptr<TestSimulationDataStruct> aTestSimData = std::shared_ptr<TestSimulationDataStruct>(new TestSimulationDataStruct);
			aTestSimData->simParameter = ShearWaveSimulationParameter::getNewInstance(kernelName, viscosity, simParaStruct, gridInfoStruct.at(i));
			aTestSimData->initialCondition = InitialConditionShearWave::getNewInstance(simParaStruct, gridInfoStruct.at(i));
			aTestSimData->simInformation = ShearWaveSimulationInfo::getNewInstance(simID, kernelName, viscosity, simParaStruct, gridInfoStruct.at(i), numberOfSimulations);
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
	std::vector<std::shared_ptr<TestSimulationImp> > testSimumlations;
	for (int i = 0; i < testSimDataStruct.size(); i++) {
		std::shared_ptr<TimeImp> time = TimeImp::getNewInstance();
		testSimDataStruct.at(i)->simInformation->setTimeInfo(time);
		std::shared_ptr<SimulationResults> simResult = SimulationResults::getNewInstance(testSimDataStruct.at(i)->simParameter);
		std::shared_ptr<ToVectorWriter> toVectorWriter = Y2dSliceToResults::getNewInstance(vectorWriterInfo, testSimDataStruct.at(i)->simParameter->getTimeStepLength(), simResult, ySliceForCalculation);
		
		testSimumlations.push_back(TestSimulationImp::getNewInsance(testSimDataStruct.at(i), simResult, time, toVectorWriter, anaResultWriter, colorOutput));
	}

	return testSimumlations;
}

std::shared_ptr<PhiTestStruct> NumericalTestFactoryImp::makePhiTestsStructs(std::shared_ptr<PhiTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp>> testSimumlations, double viscosity)
{
	std::shared_ptr<PhiTestStruct> testStruct = std::shared_ptr<PhiTestStruct>(new PhiTestStruct);

	if (testParameter->basicTestParameter->runTest && testSimumlations.size() > 1) {
		testStruct->logFileInfo = PhiTestLogFileInformation::getNewInstance(testParameter);


		std::vector<std::shared_ptr<PhiTestPostProcessingStrategy> > postProcessingStrategies;
		for (int i = 0; i < testSimumlations.size(); i++)
			postProcessingStrategies.push_back(PhiTestPostProcessingStrategy::getNewInstance(testSimumlations.at(i)->getSimulationResults(), testSimumlations.at(i)->getAnalyticalResults(), testParameter, testSimumlations.at(i)->getDataToCalcTests()));

		for (int i = 0; i < testSimumlations.at(0)->getDataToCalcTests().size(); i++) {
			std::vector<std::shared_ptr<PhiTest> > phiTests = makePhiTests(testParameter, testSimumlations, postProcessingStrategies, viscosity, testSimumlations.at(0)->getDataToCalcTests().at(i));
			testStruct->logFileInfo->addTestGroup(phiTests);
			for (int j = 0; j < phiTests.size(); j++)
				testStruct->tests.push_back(phiTests.at(j));
		}
	}

	return testStruct;
}

std::shared_ptr<NyTestStruct> NumericalTestFactoryImp::makeNyTestsStructs(std::shared_ptr<NyTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp>> testSimumlations, double viscosity)
{
	std::shared_ptr<NyTestStruct> testStruct = std::shared_ptr<NyTestStruct>(new NyTestStruct);

	if (testParameter->basicTestParameter->runTest && testSimumlations.size() > 1) {
		testStruct->logFileInfo = NyTestLogFileInformation::getNewInstance(testParameter);


		std::vector<std::shared_ptr<NyTestPostProcessingStrategy> > postProcessingStrategies;
		for (int i = 0; i < testSimumlations.size(); i++)
			postProcessingStrategies.push_back(NyTestPostProcessingStrategy::getNewInstance(testSimumlations.at(i)->getSimulationResults(), testSimumlations.at(i)->getAnalyticalResults(), testParameter, testSimumlations.at(i)->getDataToCalcTests()));

		for (int i = 0; i < testSimumlations.at(0)->getDataToCalcTests().size(); i++) {
			std::vector<std::shared_ptr<NyTest> > nyTests = makeNyTests(testParameter, testSimumlations, postProcessingStrategies, viscosity, testSimumlations.at(0)->getDataToCalcTests().at(i));
			testStruct->logFileInfo->addTestGroup(nyTests);
			for (int j = 0; j < nyTests.size(); j++)
				testStruct->tests.push_back(nyTests.at(j));
		}
	}

	return testStruct;
}

std::shared_ptr<L2NormTestStruct> NumericalTestFactoryImp::makeL2NormTestsStructs(std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSimumlations)
{
	std::shared_ptr<L2NormTestStruct> testStruct = std::shared_ptr<L2NormTestStruct> (new L2NormTestStruct);

	if (testParameter->basicTestParameter->runTest) {
		std::vector<std::shared_ptr<L2NormPostProcessingStrategy> >  postProcessingStrategies;
		for (int i = 0; i < testSimumlations.size(); i++)
			postProcessingStrategies.push_back(L2NormPostProcessingStrategy::getNewInstance(testSimumlations.at(i)->getSimulationResults(), testSimumlations.at(i)->getAnalyticalResults(), testParameter, l2NormCalculatorFactory, testSimumlations.at(i)->getDataToCalcTests()));

		testStruct->tests = makeL2NormTests(testSimumlations, postProcessingStrategies, testParameter);
		testStruct->logFileInfo = L2NormInformation::getNewInstance(testStruct->tests, testParameter, testSimumlations.at(0)->getDataToCalcTests());
	}
	return testStruct;
}

std::vector<std::shared_ptr<L2NormTest> > NumericalTestFactoryImp::makeL2NormTests(std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormPostProcessingStrategy> > postProStrategy, std::shared_ptr<L2NormTestParameterStruct> testParameter)
{
	std::vector<std::shared_ptr<L2NormTest> > l2Tests;
	for (int k = 0; k < testParameter->normalizeData.size(); k++) {
		for (int i = 0; i < testSim.size(); i++) {
			for (int j = 0; j < testSim.at(i)->getDataToCalcTests().size(); j++) {
				std::shared_ptr<L2NormTest> test = L2NormTest::getNewInstance(colorOutput, testParameter, testSim.at(i)->getDataToCalcTests().at(j), testParameter->maxDiff.at(k), testParameter->normalizeData.at(k));
				test->addSimulation(testSim.at(i), testSim.at(i)->getSimulationInfo(), postProStrategy.at(i));
				testSim.at(i)->registerSimulationObserver(test);
				l2Tests.push_back(test);
			}
		}
	}
	return l2Tests;
}

std::shared_ptr<L2NormTestBetweenKernelsStruct> NumericalTestFactoryImp::makeL2NormTestsBetweenKernelsStructs(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::string kernelName)
{
	std::shared_ptr<L2NormTestBetweenKernelsStruct> testStruct = std::shared_ptr<L2NormTestBetweenKernelsStruct>(new L2NormTestBetweenKernelsStruct);

	if (testPara->basicTestParameter->runTest) {

		std::vector<std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> > postProcessingStrategies;
		for (int i = 0; i < testSim.size(); i++)
			postProcessingStrategies.push_back(L2NormBetweenKernelPostProcessingStrategy::getNewInstance(testSim.at(i)->getSimulationResults(), testSim.at(i)->getAnalyticalResults(), testPara, l2NormCalculatorFactory, testSim.at(i)->getDataToCalcTests()));

		if (kernelName == testPara->basicKernel) {
			std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > > tests = makeL2NormTestsBetweenKernels(testPara, testSim, postProcessingStrategies);
			
			if (l2NormTestsBetweenKernels.size() == 0) {
				l2NormTestsBetweenKernels = tests;
			}
			else {
				for (int i = 0; i < tests.size(); i++)
					for (int j = 0; j < tests.at(i).size(); j++)
						l2NormTestsBetweenKernels.at(i).push_back(tests.at(i).at(j));
			}

		}else{
			std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests = linkL2NormTestsBetweenKernels(testPara, testSim, postProcessingStrategies);
			testStruct->tests = tests;
			testStruct->logFileInfo = L2NormBetweenKernelsInformation::getNewInstance(tests, testPara, testSim.at(0)->getDataToCalcTests());
		}
	}
	return testStruct;	
}

std::vector<std::shared_ptr<PhiTest>> NumericalTestFactoryImp::makePhiTests(std::shared_ptr<PhiTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp>> testSim, std::vector<std::shared_ptr<PhiTestPostProcessingStrategy>> phiPostProStrategy, double viscosity, std::string dataToCalculate)
{
	std::vector<std::shared_ptr<PhiTest> > phiTests;
	for (int i = 1; i < testSim.size(); i++) {
		for (int j = 0; j < i; j++) {
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

std::vector<std::shared_ptr<NyTest>> NumericalTestFactoryImp::makeNyTests(std::shared_ptr<NyTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp>> testSim, std::vector<std::shared_ptr<NyTestPostProcessingStrategy>> nuPostProStrategy, double viscosity, std::string dataToCalculate)
{
	std::vector<std::shared_ptr<NyTest> > nyTests;
	for (int i = 1; i < testSim.size(); i++) {
		for (int j = 0; j < i; j++) {
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


std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > >  NumericalTestFactoryImp::makeL2NormTestsBetweenKernels(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> > postProcessingStrategies)
{
	std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > > testsForAllKernels;

	std::vector<std::shared_ptr<L2NormTestBetweenKernels> > testForOneKernel;

	for (int l = 0; l < testPara->kernelsToTest.size() - 1; l++) {
		for (int k = 0; k < testSim.size(); k++) {
			for(int j = 0; j < testSim.at(k)->getDataToCalcTests().size(); j++){
				for (int m = 0; m < testPara->normalizeData.size(); m++) {
					for (int i = 0; i < testPara->timeSteps.size(); i++) {
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

	for (int k = 0; k < testSim.size(); k++) {
		for (int j = 0; j < testSim.at(k)->getDataToCalcTests().size(); j++) {
			for (int m = 0; m < testPara->normalizeData.size(); m++) {
				for (int i = 0; i < testPara->timeSteps.size(); i++) {
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

std::shared_ptr<LogFileWriter> NumericalTestFactoryImp::makeLogFileWriter(std::vector<std::shared_ptr<TestLogFileInformation> > testLogFiles, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::vector<std::shared_ptr<SimulationInfo> > simInfo, std::string kernelName, double viscosity, int basicTimeStepLength, std::shared_ptr<LogFileParameterStruct> logFilePara, std::shared_ptr<BasicTestLogFileInformation> basicTestLogFileInfo)
{
	std::shared_ptr<LogFileHead> logFileHead = LogFileHead::getNewInstance(logFilePara->devices);
	std::shared_ptr<BasicSimulationInfo> basicSimInfo = BasicSimulationInfo::getNewInstance(logFilePara->numberOfTimeSteps, viscosity, basicTimeStepLength, kernelName);

	std::shared_ptr<LogFileTimeInformation> logFileTimeInfo = LogFileTimeInformation::getNewInstance(simInfo, logFilePara->writeAnalyticalToVTK);

	std::shared_ptr<LogFileWriterImp> logFileWriter = LogFileWriterImp::getNewInstance(logFileHead, basicSimInfo, basicTestLogFileInfo, testLogFiles, logFileTimeInfo, simLogInfo, kernelName, viscosity);

	return logFileWriter;
}
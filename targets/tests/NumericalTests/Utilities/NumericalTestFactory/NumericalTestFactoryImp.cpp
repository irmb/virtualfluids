#include "NumericalTestFactoryImp.h"

#include "Simulations/TaylorGreenVortexUx/SimulationParameter/SimulationParameterTaylorGreenVortexUx.h"
#include "Simulations/TaylorGreenVortexUx/LogFileInformation/LogFileInformationTaylorGreenVortexUx.h"
#include "Simulations\TaylorGreenVortexUx\SimulationInfo\SimulationInfoTaylorGreenVortexUx.h"
#include "Simulations\TaylorGreenVortexUx\AnalyticalResults\AnalyticalResultsTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUz/SimulationParameter/SimulationParameterTaylorGreenVortexUz.h"
#include "Simulations/TaylorGreenVortexUz/LogFileInformation/LogFileInformationTaylorGreenvortexUz.h"
#include "Simulations\TaylorGreenVortexUz\SimulationInfo\SimulationInfoTaylorGreenVortexUz.h"
#include "Simulations\TaylorGreenVortexUz\AnalyticalResults\AnalyticalResultsTaylorGreenVortexUz.h"

#include "Simulations/ShearWave/SimulationParameter/ShearWaveSimulationParameter.h"
#include "Simulations/ShearWave/LogFileInformation/ShearWaveLogFileInformation.h"
#include "Simulations\ShearWave\SimulationInfo\ShearWaveSimulationInfo.h"
#include "Simulations\ShearWave\AnalyticalResults\ShearWaveAnalyticalResults.h"

#include "Tests/PhiAndNuTest/PhiAndNuTest.h"
#include "Tests\PhiAndNuTest\LogFileInformation\PhiAndNuLogFileInformation.h"
#include "Tests\PhiAndNuTest\PostProcessingStrategy\PostProcessingStrategyPhiAndNuTest.h"
#include "Tests\L2NormTest\L2NormTest.h"
#include "Tests\L2NormTest\LogFileInformation\L2NormLogFileInformation.h"
#include "Tests\L2NormTest\PostProcessingStrategy\PostProcessingStrategyL2NormTest.h"
#include "Tests\L2NormTestBetweenKernels\L2NormTestBetweenKernels.h"
#include "Tests\L2NormTestBetweenKernels\PostProcessingStrategy\L2NormBetweenKernelPostProcessingStrategy.h"
#include "Tests\L2NormTestBetweenKernels\LogFileInformation\L2NormLogFileInformationBetweenKernels.h"

#include "Utilities\ColorConsoleOutput\ColorConsoleOutputImp.h"
#include "Utilities\DataWriter\AnalyticalResults2DToVTKWriter\AnalyticalResults2DToVTKWriterImp.h"
#include "Utilities\LogFileInformation\LogFileTimeInformation\LogFileTimeInformation.h"
#include "Utilities\LogFileWriter\LogFileWriterImp.h"
#include "Utilities\LogFileQueue\LogFileQueueImp.h"
#include "Utilities\PostProcessingResults\PostProcessingResultsImp.h"
#include "Utilities\Results\SimulationResults\SimulationResults.h"
#include "Utilities\TestQueue\TestQueueImp.h"
#include "Utilities/TestSimulation/TestSimulationImp.h"

#include <algorithm>

std::shared_ptr<NumericalTestFactoryImp> NumericalTestFactoryImp::getNewInstance(std::shared_ptr<ConfigDataStruct> configFileData)
{
	return std::shared_ptr<NumericalTestFactoryImp>(new NumericalTestFactoryImp(configFileData));
}

NumericalTestFactoryImp::NumericalTestFactoryImp(std::shared_ptr<ConfigDataStruct> configFileData) : cfd(configFileData)
{
	colorOutput = ColorConsoleOutputImp::getInstance();
	testQueue = TestQueueImp::getNewInstance(colorOutput);
	logFileWriterQueue = LogFileQueueImp::getNewInstance(cfd->logFilePath);
	anaResultWriter = AnalyticalResults2DToVTKWriterImp::getInstance();
	init();
}

std::vector<std::shared_ptr<TestSimulation>> NumericalTestFactoryImp::getTestSimulations()
{
	return testSimulations;
}

std::shared_ptr<TestQueue> NumericalTestFactoryImp::getTestQueue()
{
	return testQueue;
}

std::shared_ptr<LogFileQueue> NumericalTestFactoryImp::getLogFileQueue()
{
	return logFileWriterQueue;
}

void NumericalTestFactoryImp::init()
{
	calcNumberOfSimulations();
	simID = 1;
	if (cfd->l2NormBetweenKernelTest && cfd->kernelsToTest.size() > 1)
	{
		sortKernels();
		simPerKernel = numberOfSimulations / cfd->kernelsToTest.size();
		numberOfTestGroupsBetweenKernels = (cfd->kernelsToTest.size() - 1) * simPerKernel;
		numberOfTestsForOneSimulation = cfd->dataToCalcL2NormBetweenKernel.size() * cfd->timeStepsL2NormBetweenKernel.size();
		numberOfTestsBetweenKernels = numberOfTestGroupsBetweenKernels * numberOfTestsForOneSimulation;
		posBasicSimulationForL2Test = 0;
		posDivergentSimulationForL2Test = 0;
		for (int i = 0; i < numberOfTestGroupsBetweenKernels; i++) {
			for (int j = 0; j < cfd->dataToCalcL2NormBetweenKernel.size(); j++) {
				for (int k = 0; k < cfd->timeStepsL2NormBetweenKernel.size(); k++) {
					l2KernelTests.push_back(L2NormTestBetweenKernels::getNewInstance(colorOutput, cfd->dataToCalcL2NormBetweenKernel.at(j), cfd->timeStepsL2NormBetweenKernel.at(k)));
				}
			}
		}
	}

	for (int i = 0; i < cfd->kernelsToTest.size(); i++) {
		for (int j = 0; j < cfd->viscosity.size(); j++) {
			for (int k = 0; k < cfd->u0TGVux.size(); k++)
				if (shouldSimulationGroupRun(cfd->tgvUx))
					makeTaylorGreenUxSimulations(cfd->kernelsToTest.at(i), cfd->viscosity.at(j), cfd->u0TGVux.at(k), cfd->amplitudeTGVux.at(k), cfd->basisTimeStepLengthTGVux.at(k));

			for (int k = 0; k < cfd->v0TGVuz.size(); k++)
				if (shouldSimulationGroupRun(cfd->tgvUz))
					makeTaylorGreenUzSimulations(cfd->kernelsToTest.at(i), cfd->viscosity.at(j), cfd->v0TGVuz.at(k), cfd->amplitudeTGVuz.at(k), cfd->basisTimeStepLengthTGVuz.at(k));

			for (int k = 0; k < cfd->u0SW.size(); k++)
				if (shouldSimulationGroupRun(cfd->sw))
					makeShearWaveSimulations(cfd->kernelsToTest.at(i), cfd->viscosity.at(j), cfd->u0SW.at(k), cfd->v0SW.at(k), cfd->basisTimeStepLengthSW.at(k));
			
		}
	}
}

void NumericalTestFactoryImp::makeTaylorGreenUxSimulations(std::string kernelName, double viscosity, double ux, double amplitude, int basicTimeStepLength)
{
	std::vector< std::shared_ptr< SimulationParameter>> simParaTGV;
	std::vector< std::shared_ptr< SimulationInfo>> simInfoTGV;
	std::vector< std::shared_ptr< AnalyticalResults>> analyResultTGV;

	for (int i = 0; i < cfd->tgvUx.size(); i++) {
		if (cfd->tgvUx.at(i)) {
			simParaTGV.push_back(SimulationParameterTaylorGreenUx::getNewInstance(kernelName, ux, amplitude, viscosity, cfd->rho0, cfd->lx.at(i), cfd->lz.at(i), cfd->l0TGVux, cfd->numberOfTimeSteps, basicTimeStepLength, calcStartStepForToVectorWriter(), cfd->ySliceForCalculation, cfd->grids.at(i), cfd->maxLevel, cfd->numberOfGridLevels, cfd->writeFiles, cfd->startStepFileWriter, cfd->filePath, cfd->devices));
			simInfoTGV.push_back(SimulationInfoTaylorGreenUx::getNewInstance(ux, amplitude, cfd->l0TGVux, cfd->lx.at(i), viscosity, kernelName, numberOfSimulations));
			analyResultTGV.push_back(AnalyticalResultsTaylorGreenUx::getNewInstance(viscosity, ux, amplitude, cfd->l0TGVux, cfd->rho0));
		}
	}
	std::shared_ptr<LogFileInformationTaylorGreenUx> tgInfo = LogFileInformationTaylorGreenUx::getNewInstance(ux, amplitude, cfd->tgvUx, cfd->lx, cfd->l0TGVux);

	makePeriodicBoundaryConditionSimulationAndTests(simParaTGV, simInfoTGV, analyResultTGV, tgInfo, kernelName, cfd->tgvUx, viscosity, basicTimeStepLength);
}

void NumericalTestFactoryImp::makeTaylorGreenUzSimulations(std::string kernelName, double viscosity, double uz, double amplitude, int basicTimeStepLength)
{
	std::vector< std::shared_ptr< SimulationParameter>> simParaTGV;
	std::vector< std::shared_ptr< SimulationInfo>> simInfoTGV;
	std::vector< std::shared_ptr< AnalyticalResults>> analyResultTGV;

	for (int i = 0; i < cfd->tgvUz.size(); i++) {
		if (cfd->tgvUz.at(i)) {
			simParaTGV.push_back(SimulationParameterTaylorGreenUz::getNewInstance(kernelName, uz, amplitude, viscosity, cfd->rho0, cfd->lx.at(i), cfd->lz.at(i), cfd->l0TGVuz, cfd->numberOfTimeSteps, basicTimeStepLength, calcStartStepForToVectorWriter(), cfd->ySliceForCalculation, cfd->grids.at(i), cfd->maxLevel, cfd->numberOfGridLevels, cfd->writeFiles, cfd->startStepFileWriter, cfd->filePath, cfd->devices));
			simInfoTGV.push_back(SimulationInfoTaylorGreenUz::getNewInstance(uz, amplitude, cfd->l0TGVuz, cfd->lz.at(i), viscosity, kernelName, numberOfSimulations));
			analyResultTGV.push_back(AnalyticalResultsTaylorGreenUz::getNewInstance(viscosity, uz, amplitude, cfd->l0TGVuz, cfd->rho0));
		}
	}
	std::shared_ptr< LogFileInformationTaylorGreenUz> tgInfo = LogFileInformationTaylorGreenUz::getNewInstance(uz, amplitude, cfd->tgvUz, cfd->lz, cfd->l0TGVuz);

	makePeriodicBoundaryConditionSimulationAndTests(simParaTGV, simInfoTGV, analyResultTGV, tgInfo, kernelName, cfd->tgvUz, viscosity, basicTimeStepLength);
}

void NumericalTestFactoryImp::makeShearWaveSimulations(std::string kernelName, double viscosity, double u0, double v0, int basicTimeStepLength)
{
	std::vector< std::shared_ptr< SimulationParameter>> simParaSW;
	std::vector< std::shared_ptr< SimulationInfo>> simInfoSW;
	std::vector< std::shared_ptr< AnalyticalResults>> analyResultSW;

	for (int i = 0; i < cfd->sw.size(); i++)
		if (cfd->sw.at(i)) {
			simParaSW.push_back(ShearWaveSimulationParameter::getNewInstance(kernelName, u0, v0, viscosity, cfd->rho0, cfd->lx.at(i), cfd->lz.at(i), cfd->l0SW, cfd->numberOfTimeSteps, basicTimeStepLength, calcStartStepForToVectorWriter(), cfd->ySliceForCalculation, cfd->grids.at(i), cfd->maxLevel, cfd->numberOfGridLevels, cfd->writeFiles, cfd->startStepFileWriter, cfd->filePath, cfd->devices));
			simInfoSW.push_back(ShearWaveSimulationInfo::getNewInstance(u0, v0, cfd->l0SW, cfd->lx.at(i), viscosity, kernelName, numberOfSimulations));
			analyResultSW.push_back(ShearWaveAnalyticalResults::getNewInstance(viscosity, u0, v0, cfd->l0SW, cfd->rho0));
		}

	std::shared_ptr< ShearWaveInformation> swInfo = ShearWaveInformation::getNewInstance(u0, v0, cfd->sw, cfd->lx, cfd->l0SW);

	makePeriodicBoundaryConditionSimulationAndTests(simParaSW, simInfoSW, analyResultSW, swInfo, kernelName, cfd->sw, viscosity, basicTimeStepLength);
}

void NumericalTestFactoryImp::makePeriodicBoundaryConditionSimulationAndTests(std::vector<std::shared_ptr<SimulationParameter>> simPara, std::vector<std::shared_ptr<SimulationInfo>> simInfo, std::vector<std::shared_ptr<AnalyticalResults>> analyResult, std::shared_ptr<SimulationLogFileInformation> simlogFileInfo, std::string kernelName, std::vector<bool> simulationsRun, double viscosity, int basicTimeStepLength)
{
	std::vector< std::shared_ptr< SimulationResults>> simResults;
	for (int i = 0; i < simPara.size(); i++) {
		std::shared_ptr< SimulationResults> simResult = SimulationResults::getNewInstance(simPara.at(i)->getLx(), 1, simPara.at(i)->getLz(), simPara.at(i)->getTimeStepLength());
		simResults.push_back(simResult);
	}

	std::vector< std::shared_ptr< TestSimulation>> testSim = buildTestSimulation(simPara, simInfo, simResults, analyResult);

	std::vector< std::shared_ptr< TestLogFileInformation>> testLogFileInfo;

	if (cfd->nuAndPhiTest && checkNuAndPhiTestCouldRun(simulationsRun)) {
		std::vector< std::shared_ptr< PhiAndNuTestPostProcessingStrategy>> phiAndNuPostProStrategy;
		for(int i = 0; i < testSim.size();i++)
			phiAndNuPostProStrategy.push_back(PhiAndNuTestPostProcessingStrategy::getNewInstance(simResults.at(i), analyResult.at(i), cfd->dataToCalcPhiAndNuTest, cfd->startTimeStepCalculationPhiNu, cfd->endTimeStepCalculationPhiNu));
		std::shared_ptr< PhiAndNuInformation> phiNuLogFileInfo = PhiAndNuInformation::getNewInstance(cfd->startTimeStepCalculationPhiNu, cfd->endTimeStepCalculationPhiNu);
		for (int i = 0; i < cfd->dataToCalcPhiAndNuTest.size(); i++) {
			std::vector< std::shared_ptr< PhiAndNuTest>> phiAndNuTests = makePhiAndNuTests(testSim, simInfo, phiAndNuPostProStrategy, viscosity, cfd->dataToCalcPhiAndNuTest.at(i));
			phiNuLogFileInfo->addTestGroup(phiAndNuTests);
		}
		testLogFileInfo.push_back(phiNuLogFileInfo);
	}

	if (cfd->l2NormTest) {
		std::vector< std::shared_ptr< L2NormTest>> l2NormTests = makeL2NormTests(testSim, simInfo, simResults, analyResult);
		std::shared_ptr< L2NormInformation> l2NormLogFileInfo = L2NormInformation::getNewInstance(l2NormTests, cfd->basicTimeStepL2Norm, cfd->divergentTimeStepL2Norm);
		testLogFileInfo.push_back(l2NormLogFileInfo);
	}

	if (cfd->l2NormBetweenKernelTest){
		std::vector< std::shared_ptr< L2NormTestBetweenKernels>> tests = makeL2NormTestsBetweenKernels(testSim, simInfo, simResults, analyResult);
		if (tests.size() > 0){
			std::shared_ptr<L2NormBetweenKernelsInformation> l2NormBetweenKernelLogFileInfo = L2NormBetweenKernelsInformation::getNewInstance(tests, cfd->basicKernelL2NormTest, cfd->timeStepsL2NormBetweenKernel, cfd->dataToCalcL2NormBetweenKernel);
			testLogFileInfo.push_back(l2NormBetweenKernelLogFileInfo);
		}
	}

	for (int i = 0; i < testSim.size(); i++)
		testSimulations.push_back(testSim.at(i));

	std::shared_ptr< LogFileTimeInformation> logFileTimeInfo = LogFileTimeInformation::getNewInstance(testSim, cfd->writeFiles);

	makeLogFileWriter(testLogFileInfo, logFileTimeInfo, simlogFileInfo, kernelName, viscosity, basicTimeStepLength);
}

std::vector<std::shared_ptr<TestSimulation>> NumericalTestFactoryImp::buildTestSimulation(std::vector< std::shared_ptr< SimulationParameter>> simPara, std::vector< std::shared_ptr< SimulationInfo>> simInfo, std::vector< std::shared_ptr< SimulationResults>> simResults, std::vector< std::shared_ptr< AnalyticalResults>> analyResult)
{
	std::vector< std::shared_ptr< TestSimulation>> testSim;
	testSim.resize(0);

	for (int i = 0; i < simPara.size(); i++) {
		testSim.push_back(TestSimulationImp::getNewInsance(simID, simPara.at(i), simInfo.at(i), colorOutput, simResults.at(i), analyResult.at(i), anaResultWriter, cfd->writeAnalyticalToVTK));
		simID++;
	}
	return testSim;
}

std::vector<std::shared_ptr<PhiAndNuTest>> NumericalTestFactoryImp::makePhiAndNuTests(std::vector<std::shared_ptr<TestSimulation>> testSim, std::vector<std::shared_ptr<SimulationInfo>> simInfo, std::vector< std::shared_ptr< PhiAndNuTestPostProcessingStrategy>> phiAndNuPostProStrategy, double viscosity, std::string dataToCalculate)
{
	std::vector< std::shared_ptr< PhiAndNuTest>> phiAndNuTests;

	for (int i = 1; i < testSim.size(); i++) {
		for (int j = 0; j < i; j++) {
				std::shared_ptr< PhiAndNuTest> test = PhiAndNuTest::getNewInstance(colorOutput, dataToCalculate, cfd->minOrderOfAccuracy, viscosity, cfd->startTimeStepCalculationPhiNu, cfd->endTimeStepCalculationPhiNu);
				test->addSimulation(testSim.at(j), simInfo.at(j), phiAndNuPostProStrategy.at(j));
				test->addSimulation(testSim.at(i), simInfo.at(i), phiAndNuPostProStrategy.at(i));

				testSim.at(j)->registerSimulationObserver(test);
				testSim.at(i)->registerSimulationObserver(test);

				phiAndNuTests.push_back(test);
				testQueue->addTest(test);
		}
	}
	return phiAndNuTests;
}

std::vector<std::shared_ptr<L2NormTest>> NumericalTestFactoryImp::makeL2NormTests(std::vector<std::shared_ptr<TestSimulation>> testSim, std::vector<std::shared_ptr<SimulationInfo>> simInfo, std::vector< std::shared_ptr< SimulationResults>> simResults, std::vector<std::shared_ptr<AnalyticalResults>> analyResult)
{
	std::vector<std::shared_ptr<L2NormTest>> l2Tests;
	for (int i = 0; i < testSim.size(); i++) {
		std::shared_ptr< L2NormPostProcessingStrategy> l2NormPostProStrategy = L2NormPostProcessingStrategy::getNewInstance(simResults.at(i), analyResult.at(i), cfd->dataToCalcL2Test, cfd->basicTimeStepL2Norm, cfd->divergentTimeStepL2Norm);
		for (int j = 0; j < cfd->dataToCalcL2Test.size(); j++) {
			std::shared_ptr<L2NormTest> test = L2NormTest::getNewInstance(colorOutput, cfd->dataToCalcL2Test.at(j), cfd->maxL2NormDiff, cfd->basicTimeStepL2Norm, cfd->divergentTimeStepL2Norm);
			test->addSimulation(testSim.at(i), simInfo.at(i), l2NormPostProStrategy);
			testSim.at(i)->registerSimulationObserver(test);
			l2Tests.push_back(test);
			testQueue->addTest(test);
		}
		
	}
	return l2Tests;
}

std::vector<std::shared_ptr<L2NormTestBetweenKernels>> NumericalTestFactoryImp::makeL2NormTestsBetweenKernels(std::vector<std::shared_ptr<TestSimulation>> testSim, std::vector<std::shared_ptr<SimulationInfo>> simInfo, std::vector<std::shared_ptr<SimulationResults>> simResults, std::vector<std::shared_ptr<AnalyticalResults>> analyResult)
{
	std::vector< std::shared_ptr< L2NormTestBetweenKernels>> tests;
	for (int i = 0; i < testSim.size(); i++)
	{
		std::shared_ptr< L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy = L2NormBetweenKernelPostProcessingStrategy::getNewInstance(simResults.at(i), analyResult.at(i), cfd->timeStepsL2NormBetweenKernel, cfd->dataToCalcL2NormBetweenKernel);
		if (simID - 1 <= simPerKernel)
		{
			for (int j = 0; j < numberOfTestsForOneSimulation; j++)
			{
				for (int k = 0; k < cfd->kernelsToTest.size() - 1; k++)
				{
					l2KernelTests.at(posBasicSimulationForL2Test + simPerKernel * k)->setBasicSimulation(testSim.at(i), simInfo.at(i), postProcessingStrategy);
					testSim.at(i)->registerSimulationObserver(l2KernelTests.at(posBasicSimulationForL2Test + simPerKernel * k));
				}
				posBasicSimulationForL2Test++;
			}
		}
		else
		{
			for (int j = 0; j < numberOfTestsForOneSimulation; j++)
			{
				l2KernelTests.at(posDivergentSimulationForL2Test)->setDivergentKernelSimulation(testSim.at(i), simInfo.at(i), postProcessingStrategy);
				tests.push_back(l2KernelTests.at(posDivergentSimulationForL2Test));
				testQueue->addTest(l2KernelTests.at(posDivergentSimulationForL2Test));
				testSim.at(i)->registerSimulationObserver(l2KernelTests.at(posDivergentSimulationForL2Test));
				posDivergentSimulationForL2Test++;
			}
		}
	}
	return tests;
}

void NumericalTestFactoryImp::makeLogFileWriter(std::vector<std::shared_ptr<TestLogFileInformation>> testLogFiles, std::shared_ptr<LogFileTimeInformation> logFileTimeInfo, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::string kernelName, double viscosity, int basicTimeStepLength)
{
	std::shared_ptr< LogFileWriterImp> logFileWriter = LogFileWriterImp::getNewInstance(testLogFiles, logFileTimeInfo, simLogInfo, kernelName, viscosity, cfd->devices, cfd->numberOfTimeSteps, calcStartStepForToVectorWriter(), basicTimeStepLength);
	logFileWriterQueue->addLogFileWriter(logFileWriter);
}

void NumericalTestFactoryImp::sortKernels()
{
	bool basicKernelInKernelList = false;
	for (int i = 0; i < cfd->kernelsToTest.size(); i++) {
		if (cfd->kernelsToTest.at(i) == cfd->basicKernelL2NormTest)
			basicKernelInKernelList = true;
	}
	if (!basicKernelInKernelList)
		cfd->kernelsToTest.push_back(cfd->basicKernelL2NormTest);

	std::vector< std::string> kernels = cfd->kernelsToTest;

	while (kernels.at(0)!= cfd->basicKernelL2NormTest){
		kernels.push_back(kernels.at(0));
		std::vector< std::string>::iterator it = kernels.begin();
		kernels.erase(it);
	}
	cfd->kernelsToTest = kernels;
}

bool NumericalTestFactoryImp::shouldSimulationGroupRun(std::vector<bool> test)
{
	for (int i = 0; i < test.size(); i++) {
		if (test.at(i))
			return true;
	}
	return false;
}

unsigned int NumericalTestFactoryImp::calcStartStepForToVectorWriter()
{
	std::vector< unsigned int> startStepsTests;
	startStepsTests.push_back(cfd->basicTimeStepL2Norm);
	startStepsTests.push_back(cfd->startTimeStepCalculationPhiNu);

	std::sort(startStepsTests.begin(), startStepsTests.end());

	return startStepsTests.at(0);
}

bool NumericalTestFactoryImp::checkNuAndPhiTestCouldRun(std::vector<bool> test)
{
	int numberOfTestInGroup = 0;
	for (int i = 0; i < test.size(); i++) {
		if (test.at(i))
			numberOfTestInGroup++;
	}
	return numberOfTestInGroup > 1;
}

void NumericalTestFactoryImp::calcNumberOfSimulations()
{
	int counter = 0;

	int tgvCounterU0 = 0;
	for (int i = 0; i < cfd->tgvUx.size(); i++)
		if (cfd->tgvUx.at(i))
			tgvCounterU0++;
	tgvCounterU0 *= cfd->u0TGVux.size();
	counter += tgvCounterU0;

	int tgvCounterV0 = 0;
	for (int i = 0; i < cfd->tgvUz.size(); i++)
		if (cfd->tgvUz.at(i))
			tgvCounterV0++;
	tgvCounterV0 *= cfd->v0TGVuz.size();
	counter += tgvCounterV0;


	int swCounter = 0;
	for (int i = 0; i < cfd->sw.size(); i++)
		if (cfd->sw.at(i))
			swCounter++;
	swCounter *= cfd->u0SW.size();
	counter += swCounter;

	counter *= cfd->viscosity.size();
	counter *= cfd->kernelsToTest.size();

	numberOfSimulations = counter;
}
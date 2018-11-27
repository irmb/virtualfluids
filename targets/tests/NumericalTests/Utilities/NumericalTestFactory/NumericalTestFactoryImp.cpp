#include "NumericalTestFactoryImp.h"

#include "Simulation/TaylorGreenVortex/SimulationParameter/TaylorGreenSimulationParameter.h"
#include "Simulation/TaylorGreenVortex/LogFileInformation/TaylorGreenLogFileInformation.h"
#include "Simulation\TaylorGreenVortex\SimulationInfo\TaylorGreenVortexSimulationInfo.h"
#include "Simulation\TaylorGreenVortex\AnalyticalResults\TaylorGreenVortexAnalyticalResults.h"

#include "Simulation/ShearWave/SimulationParameter/ShearWaveSimulationParameter.h"
#include "Simulation/ShearWave/LogFileInformation/ShearWaveLogFileInformation.h"
#include "Simulation\ShearWave\SimulationInfo\ShearWaveSimulationInfo.h"
#include "Simulation\ShearWave\AnalyticalResults\ShearWaveAnalyticalResults.h"

#include "Tests/PhiAndNuTest/PhiAndNuTest.h"
#include "Tests\PhiAndNuTest\LogFileInformation\PhiAndNuLogFileInformation.h"
#include "Tests\L2NormTest\L2NormTest.h"
#include "Tests\L2NormTest\LogFileInformation\L2NormLogFileInformation.h"

#include "Utilities\ColorConsoleOutput\ColorConsoleOutputImp.h"
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
	for (int i = 0; i < cfd->kernelsToTest.size(); i++) {
		for (int j = 0; j < cfd->viscosity.size(); j++) {
			for (int k = 0; k < cfd->u0TGV.size(); k++) {
				if (shouldSimulationGroupRun(cfd->tgv)) {
					makeTaylorGreenSimulations(cfd->kernelsToTest.at(i), cfd->viscosity.at(j), cfd->u0TGV.at(k), cfd->amplitudeTGV.at(k));
				}
			}
			for (int k = 0; k < cfd->u0SW.size(); k++) {
				if (shouldSimulationGroupRun(cfd->sw))
					makeShearWaveSimulations(cfd->kernelsToTest.at(i), cfd->viscosity.at(j), cfd->u0SW.at(k), cfd->v0SW.at(k));
			}
		}
	}
}

void NumericalTestFactoryImp::makeTaylorGreenSimulations(std::string kernelName, double viscosity, double u0, double amplitude)
{
	std::vector< std::shared_ptr< SimulationParameter>> simParaTGV;
	std::vector< std::shared_ptr< SimulationInfo>> simInfoTGV;
	std::vector< std::shared_ptr< AnalyticalResults>> analyResultTGV;

	for (int i = 0; i < cfd->tgv.size(); i++) {
		if (cfd->tgv.at(i)) {
			simParaTGV.push_back(TaylorGreenSimulationParameter::getNewInstance(kernelName, u0, amplitude, viscosity, cfd->rho0, cfd->lx.at(i), cfd->lz.at(i), cfd->l0, cfd->numberOfTimeSteps, cfd->basisTimeStepLength, calcStartStepForToVectorWriter(), cfd->ySliceForCalculation, cfd->grids.at(i), cfd->maxLevel, cfd->numberOfGridLevels, cfd->writeFiles, cfd->startStepFileWriter, cfd->filePath, cfd->devices));
			simInfoTGV.push_back(TaylorGreenVortexSimulationInfo::getNewInstance(u0, amplitude, cfd->l0, cfd->lx.at(i), viscosity, kernelName, numberOfSimulations));
			analyResultTGV.push_back(TaylorGreenAnalyticalResults::getNewInstance(viscosity, u0, amplitude, cfd->l0, cfd->rho0));
		}
	}
	std::shared_ptr< TaylorGreenInformation> tgInfo = TaylorGreenInformation::getNewInstance(u0, amplitude, cfd->tgv, cfd->lx, cfd->l0);

	makePeriodicBoundaryConditionSimulationAndTests(simParaTGV, simInfoTGV, analyResultTGV, tgInfo, kernelName, cfd->tgv, viscosity, cfd->nuAndPhiTestTGV, cfd->l2NormTestTGV);
}

void NumericalTestFactoryImp::makeShearWaveSimulations(std::string kernelName, double viscosity, double u0, double v0)
{
	std::vector< std::shared_ptr< SimulationParameter>> simParaSW;
	std::vector< std::shared_ptr< SimulationInfo>> simInfoSW;
	std::vector< std::shared_ptr< AnalyticalResults>> analyResultSW;

	for (int i = 0; i < cfd->sw.size(); i++)
		if (cfd->sw.at(i)) {
			simParaSW.push_back(ShearWaveSimulationParameter::getNewInstance(kernelName, u0, v0, viscosity, cfd->rho0, cfd->lx.at(i), cfd->lz.at(i), cfd->l0, cfd->numberOfTimeSteps, cfd->basisTimeStepLength, calcStartStepForToVectorWriter(), cfd->ySliceForCalculation, cfd->grids.at(i), cfd->maxLevel, cfd->numberOfGridLevels, cfd->writeFiles, cfd->startStepFileWriter, cfd->filePath, cfd->devices));
			simInfoSW.push_back(ShearWaveSimulationInfo::getNewInstance(u0, v0, cfd->l0, cfd->lx.at(i), viscosity, kernelName, numberOfSimulations));
			analyResultSW.push_back(ShearWaveAnalyticalResults::getNewInstance(viscosity, u0, v0, cfd->l0, cfd->rho0));
		}

	std::shared_ptr< ShearWaveInformation> swInfo = ShearWaveInformation::getNewInstance(u0, v0, cfd->sw, cfd->lx, cfd->l0);

	makePeriodicBoundaryConditionSimulationAndTests(simParaSW, simInfoSW, analyResultSW, swInfo, kernelName, cfd->sw, viscosity, cfd->nuAndPhiTestSW, cfd->l2NormTestSW);
}

void NumericalTestFactoryImp::makePeriodicBoundaryConditionSimulationAndTests(std::vector<std::shared_ptr<SimulationParameter>> simPara, std::vector<std::shared_ptr<SimulationInfo>> simInfo, std::vector<std::shared_ptr<AnalyticalResults>> analyResult, std::shared_ptr<SimulationLogFileInformation> simlogFileInfo, std::string kernelName, std::vector<bool> simulationsRun, double viscosity, bool nuAndPhiTest, bool l2NormTest)
{
	std::vector< std::shared_ptr< SimulationResults>> simResults;
	std::vector< std::shared_ptr< PostProcessingResults>> postProResults;
	for (int i = 0; i < simPara.size(); i++) {
		std::shared_ptr< SimulationResults> simResult = SimulationResults::getNewInstance(simPara.at(i)->getLx(), 1, simPara.at(i)->getLz(), simPara.at(i)->getTimeStepLength());
		simResults.push_back(simResult);
		postProResults.push_back(PostProcessingResultsImp::getNewInstance(simResult, analyResult.at(i), cfd->dataToCalcPhiAndNuTest, cfd->dataToCalcL2Test, cfd->startTimeStepCalculationPhiNu, cfd->endTimeStepCalculationPhiNu, cfd->basicTimeStepL2Norm, cfd->divergentTimeStepL2Norm));
	}

	std::vector< std::shared_ptr< TestSimulation>> testSim = buildTestSimulation(simPara, simInfo, simResults);

	std::vector< std::shared_ptr< TestLogFileInformation>> testLogFileInfo;

	if (nuAndPhiTest && checkNuAndPhiTestCouldRun(simulationsRun)) {
		std::vector< std::shared_ptr< PhiAndNuTest>> phiAndNuTests = makePhiAndNuTests(testSim, simInfo, postProResults, viscosity);
		std::shared_ptr< PhiAndNuInformation> phiNuLogFileInfo = PhiAndNuInformation::getNewInstance(phiAndNuTests, cfd->startTimeStepCalculationPhiNu, cfd->endTimeStepCalculationPhiNu);
		testLogFileInfo.push_back(phiNuLogFileInfo);
	}

	if (l2NormTest) {
		std::vector< std::shared_ptr< L2NormTest>> l2NormTests = makeL2NormTests(testSim, simInfo, postProResults);
		std::shared_ptr< L2NormInformation> l2NormLogFileInfo = L2NormInformation::getNewInstance(l2NormTests, cfd->basicTimeStepL2Norm, cfd->divergentTimeStepL2Norm);
		testLogFileInfo.push_back(l2NormLogFileInfo);
	}

	for (int i = 0; i < testSim.size(); i++)
		testSimulations.push_back(testSim.at(i));

	std::shared_ptr< LogFileTimeInformation> logFileTimeInfo = LogFileTimeInformation::getNewInstance(testSim, cfd->writeFiles);

	makeLogFileWriter(testLogFileInfo, logFileTimeInfo, simlogFileInfo, kernelName, viscosity);
}

std::vector<std::shared_ptr<TestSimulation>> NumericalTestFactoryImp::buildTestSimulation(std::vector< std::shared_ptr< SimulationParameter>> simPara, std::vector< std::shared_ptr< SimulationInfo>> simInfo, std::vector< std::shared_ptr< SimulationResults>> simResults)
{
	std::vector< std::shared_ptr< TestSimulation>> testSim;
	testSim.resize(0);

	for (int i = 0; i < simPara.size(); i++) {
		testSim.push_back(TestSimulationImp::getNewInsance(simID, simPara.at(i), simInfo.at(i), colorOutput, simResults.at(i)));
		simID++;
	}
	return testSim;
}

std::vector<std::shared_ptr<PhiAndNuTest>> NumericalTestFactoryImp::makePhiAndNuTests(std::vector<std::shared_ptr<TestSimulation>> testSim, std::vector<std::shared_ptr<SimulationInfo>> simInfo, std::vector< std::shared_ptr< PostProcessingResults>> postProResults, double viscosity)
{
	std::vector< std::shared_ptr< PhiAndNuTest>> phiAndNuTests;

	for (int i = 1; i < testSim.size(); i++) {
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < cfd->dataToCalcPhiAndNuTest.size(); k++) {
				std::shared_ptr< PhiAndNuTest> test = PhiAndNuTest::getNewInstance(colorOutput, cfd->dataToCalcPhiAndNuTest.at(k), cfd->minOrderOfAccuracy, viscosity, cfd->startTimeStepCalculationPhiNu, cfd->endTimeStepCalculationPhiNu);
				test->addSimulation(testSim.at(j), simInfo.at(j), postProResults.at(j));
				test->addSimulation(testSim.at(i), simInfo.at(i), postProResults.at(i));

				testSim.at(j)->registerSimulationObserver(test);
				testSim.at(i)->registerSimulationObserver(test);

				phiAndNuTests.push_back(test);
				testQueue->addTest(test);
			}
		}
	}
	return phiAndNuTests;
}

std::vector<std::shared_ptr<L2NormTest>> NumericalTestFactoryImp::makeL2NormTests(std::vector<std::shared_ptr<TestSimulation>> testSim, std::vector<std::shared_ptr<SimulationInfo>> simInfo, std::vector< std::shared_ptr< PostProcessingResults>> postProResults)
{
	std::vector<std::shared_ptr<L2NormTest>> l2Tests;
	for (int i = 0; i < testSim.size(); i++) {
		for (int j = 0; j < cfd->dataToCalcL2Test.size(); j++) {
			std::shared_ptr<L2NormTest> test = L2NormTest::getNewInstance(postProResults.at(i), colorOutput, cfd->dataToCalcL2Test.at(j), cfd->maxL2NormDiff, cfd->basicTimeStepL2Norm, cfd->divergentTimeStepL2Norm);
			test->addSimulation(testSim.at(i), simInfo.at(i), postProResults.at(i));
			testSim.at(i)->registerSimulationObserver(test);
			l2Tests.push_back(test);
			testQueue->addTest(test);
		}
		
	}
	return l2Tests;
}

void NumericalTestFactoryImp::makeLogFileWriter(std::vector<std::shared_ptr<TestLogFileInformation>> testLogFiles, std::shared_ptr<LogFileTimeInformation> logFileTimeInfo, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::string kernelName, double viscosity)
{
	std::shared_ptr< LogFileWriterImp> logFileWriter = LogFileWriterImp::getNewInstance(testLogFiles, logFileTimeInfo, simLogInfo, kernelName, viscosity, cfd->devices, cfd->numberOfTimeSteps, cfd->basisTimeStepLength, calcStartStepForToVectorWriter());
	logFileWriterQueue->addLogFileWriter(logFileWriter);
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

	int tgvCounter = 0;
	for (int i = 0; i < cfd->tgv.size(); i++)
		if (cfd->tgv.at(i))
			tgvCounter++;
	tgvCounter *= cfd->u0TGV.size();
	counter += tgvCounter;

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

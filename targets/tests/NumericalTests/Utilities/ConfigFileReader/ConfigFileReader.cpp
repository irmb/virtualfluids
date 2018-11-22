#include "ConfigFileReader.h"

#include <fstream>
#include <iostream>

#include "utilities/input/Input.h"
#include "utilities/StringUtil/StringUtil.h"

#include "Utilities/TestSimulation/TestSimulationImp.h"
#include "Utilities\TestQueue\TestQueueImp.h"

#include "Utilities\LogFileWriter\LogFileWriter.h"
#include "Utilities\LogFileQueue\LogFileQueueImp.h"

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

#include "Utilities/LogFileInformation/LogFileInformation.h"
#include "Utilities/LogFileInformation/BasicSimulationInfo/BasicSimulationInfo.h"
#include "Utilities\LogFileInformation\LogFileTimeInformation\LogFileTimeInformation.h"

#include "Utilities\ColorConsoleOutput\ColorConsoleOutputImp.h"

std::shared_ptr<ConfigFileReader> ConfigFileReader::getNewInstance()
{
	return std::shared_ptr<ConfigFileReader>(new ConfigFileReader());
}

ConfigFileReader::ConfigFileReader()
{
	logInfo.resize(0);

	testSimulation.resize(0);

	lx.resize(5);
	lx.at(0) = 32.0;
	lx.at(1) = 64.0;
	lx.at(2) = 128.0;
	lx.at(3) = 256.0;
	lx.at(4) = 512.0;

	lz.resize(5);
	lz.at(0) = lx.at(0) * 3.0 / 2.0;
	lz.at(1) = lx.at(1) * 3.0 / 2.0;
	lz.at(2) = lx.at(2) * 3.0 / 2.0;
	lz.at(3) = lx.at(3) * 3.0 / 2.0;
	lz.at(4) = lx.at(4) * 3.0 / 2.0;

	l0 = 32.0;
	rho0 = 1.0;


	colorOutput = ColorConsoleOutputImp::getInstance();
	testQueue = TestQueueImp::getNewInstance(colorOutput);
}

void ConfigFileReader::readConfigFile(const std::string aFilePath)
{
	std::ifstream stream;
	stream.open(aFilePath.c_str(), std::ios::in);
	if (stream.fail())
		throw "can not open config file!\n";

	std::unique_ptr<input::Input> input = input::Input::makeInput(stream, "config");

	devices = StringUtil::toIntVector(input->getValue("Devices"));
	kernelsToTest = StringUtil::toStringVector(input->getValue("KernelsToTest"));

	viscosity = StringUtil::toDoubleVector(input->getValue("Viscosity"));

	minOrderOfAccuracy = StringUtil::toDouble(input->getValue("MinOrderOfAccuracy"));
	dataToCalcPhiAndNuTest = StringUtil::toString(input->getValue("DataToCalcPhiAndNuTest"));
	startStepCalculationPhiNu = StringUtil::toInt(input->getValue("StartTimeStepCalculation_PhiNu"));
	endStepCalculationPhiNu = StringUtil::toInt(input->getValue("EndTimeStepCalculation_PhiNu"));

	basicTimeStepL2Norm = StringUtil::toInt(input->getValue("BasicTimeStep_L2"));
	divergentDataL2Norm = StringUtil::toInt(input->getValue("DivergentDataTimeStep_L2"));

	amplitudeTGV = StringUtil::toDoubleVector(input->getValue("Amplitude_TGV"));
	u0TGV = StringUtil::toDoubleVector(input->getValue("u0_TGV"));
	nuAndPhiTestTGV = StringUtil::toBool(input->getValue("PhiAndNuTest_TGV"));
	l2NormTestTGV = StringUtil::toBool(input->getValue("L2NormTest_TGV"));

	v0SW = StringUtil::toDoubleVector(input->getValue("v0_SW"));
	u0SW = StringUtil::toDoubleVector(input->getValue("u0_SW"));
	nuAndPhiTestSW = StringUtil::toBool(input->getValue("PhiAndNuTest_SW"));
	l2NormTestSW = StringUtil::toBool(input->getValue("L2NormTest_SW"));

	numberOfTimeSteps = StringUtil::toInt(input->getValue("NumberOfTimeSteps"));
	basisTimeStepLength = StringUtil::toInt(input->getValue("BasisTimeStepLength"));
	

	grids.resize(5);
	grids.at(0) = input->getValue("GridPath32");
	grids.at(1) = input->getValue("GridPath64");
	grids.at(2) = input->getValue("GridPath128");
	grids.at(3) = input->getValue("GridPath256");
	grids.at(4) = input->getValue("GridPath512");

	numberOfGridLevels = StringUtil::toInt(input->getValue("NumberOfGridLevels"));
	maxLevel = numberOfGridLevels - 1;

	ySliceForCalculation = StringUtil::toInt(input->getValue("ySliceForCalculation"));

	writeFiles = StringUtil::toBool(input->getValue("WriteFiles"));
	filePath = input->getValue("PathForFileWriting");
	startStepFileWriter = StringUtil::toInt(input->getValue("StartStepFileWriter"));
	logFilePath = input->getValue("PathLogFile");;

	tgv.resize(5);
	tgv.at(0) = StringUtil::toBool(input->getValue("TaylorGreenVortex32"));
	tgv.at(1) = StringUtil::toBool(input->getValue("TaylorGreenVortex64"));
	tgv.at(2) = StringUtil::toBool(input->getValue("TaylorGreenVortex128"));
	tgv.at(3) = StringUtil::toBool(input->getValue("TaylorGreenVortex256"));
	tgv.at(4) = StringUtil::toBool(input->getValue("TaylorGreenVortex512"));

	sw.resize(5);
	sw.at(0) = StringUtil::toBool(input->getValue("ShearWave32"));
	sw.at(1) = StringUtil::toBool(input->getValue("ShearWave64"));
	sw.at(2) = StringUtil::toBool(input->getValue("ShearWave128"));
	sw.at(3) = StringUtil::toBool(input->getValue("ShearWave256"));
	sw.at(4) = StringUtil::toBool(input->getValue("ShearWave512"));

	stream.close();

	checkConfigFileData();
	logFileWriterQueue = LogFileQueueImp::getNewInstance(logFilePath);

	init();
}

void ConfigFileReader::checkConfigFileData()
{
	if (u0TGV.size() != amplitudeTGV.size()) {
		std::cout << "Length u0_TGV is unequal to Lenght Amplitude_TGV!" << std::endl << std::flush;
		exit(1);
	}
		
	if (u0SW.size() != v0SW.size()) {
		std::cout << "Length u0_SW is unequal to Lenght v0_SW!" << std::endl << std::flush;
		exit(1);
	}	
}

void ConfigFileReader::init()
{
	for (int i = 0; i < kernelsToTest.size(); i++) {
		simID = 1;
		for (int j = 0; j < viscosity.size(); j++) {
			for (int k = 0; k < u0TGV.size(); k++) {
				if (shouldSimulationGroupRun(tgv))
					makeTaylorGreenSimulations(kernelsToTest.at(i), viscosity.at(j), u0TGV.at(k), amplitudeTGV.at(k));
			}
			for (int k = 0; k < u0SW.size(); k++) {
				if (shouldSimulationGroupRun(sw))
					makeShearWaveSimulations(kernelsToTest.at(i), viscosity.at(j), u0SW.at(k), v0SW.at(k));
			}
		}
	}
	
}

std::vector< std::shared_ptr< TestSimulation>> ConfigFileReader::buildTestSimulation(std::vector< std::shared_ptr< SimulationParameter>> simPara, std::vector< std::shared_ptr< SimulationInfo>> simInfo)
{
	std::vector< std::shared_ptr< TestSimulation>> testSim;
	testSim.resize(0);

	for (int i = 0; i < simPara.size(); i++) {
		testSim.push_back(TestSimulationImp::getNewInsance(simID, simPara.at(i), simInfo.at(i), colorOutput));
		simID++;
	}
	return testSim;
}

void ConfigFileReader::makeTaylorGreenSimulations(std::string kernelName, double viscosity, double u0, double amplitude)
{
	std::vector< std::shared_ptr< SimulationParameter>> simParaTGV;
	simParaTGV.resize(0);
	std::vector< std::shared_ptr< SimulationInfo>> simInfoTGV;
	simInfoTGV.resize(0);
	std::vector< std::shared_ptr< AnalyticalResults>> analyResultTGV;
	analyResultTGV.resize(0);

	for (int i = 0; i < tgv.size(); i++)
		if (tgv.at(i)) {
			simParaTGV.push_back(TaylorGreenSimulationParameter::getNewInstance(kernelName, u0, amplitude, viscosity, rho0, lx.at(i), lz.at(i), l0, numberOfTimeSteps, basisTimeStepLength, calcStartStepForToVectorWriter(), ySliceForCalculation, grids.at(i), maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, filePath, devices));
			simInfoTGV.push_back(TaylorGreenVortexSimulationInfo::getNewInstance(u0, amplitude, l0, lx.at(i), viscosity, kernelName, "TaylorGreenVortex"));
			analyResultTGV.push_back(TaylorGreenAnalyticalResults::getNewInstance(viscosity, u0, amplitude, l0, rho0));
		}

	std::vector< std::shared_ptr< TestSimulation>> testSimTGV = buildTestSimulation(simParaTGV, simInfoTGV);

	std::vector< std::shared_ptr< TestLogFileInformation>> testLogFileInfo;

	if (nuAndPhiTestTGV) {
		std::vector< std::shared_ptr< PhiAndNuTest>> phiAndNuTests = makePhiAndNuTests(testSimTGV, simInfoTGV, viscosity);
		std::shared_ptr< PhiAndNuInformation> phiNuLogFileInfo = PhiAndNuInformation::getNewInstance(phiAndNuTests);
		testLogFileInfo.push_back(phiNuLogFileInfo);
	}

	if (l2NormTestTGV) {
		std::vector< std::shared_ptr< L2NormTest>> l2NormTests = makeL2NormTests(testSimTGV, simInfoTGV, analyResultTGV);
	}
	
		

	for (int i = 0; i < testSimTGV.size(); i++)
		testSimulation.push_back(testSimTGV.at(i));

	std::shared_ptr< TaylorGreenInformation> tgInfo = TaylorGreenInformation::getNewInstance(u0, amplitude, tgv, lx, l0);
	std::shared_ptr< LogFileTimeInformation> logFileTimeInfo = LogFileTimeInformation::getNewInstance(testSimTGV, writeFiles);

	makeLogFileWriter(testLogFileInfo, logFileTimeInfo, tgInfo, kernelName, viscosity);
}

void ConfigFileReader::makeShearWaveSimulations(std::string kernelName, double viscosity, double u0, double v0)
{
	std::vector< std::shared_ptr< SimulationParameter>> simParaSW;
	simParaSW.resize(0);
	std::vector< std::shared_ptr< SimulationInfo>> simInfoSW;
	simInfoSW.resize(0);

	for (int i = 0; i < tgv.size(); i++)
		if (tgv.at(i)) {
			simParaSW.push_back(ShearWaveSimulationParameter::getNewInstance(kernelName, u0, v0, viscosity, rho0, lx.at(i), lz.at(i), l0, numberOfTimeSteps, basisTimeStepLength, calcStartStepForToVectorWriter(), ySliceForCalculation, grids.at(i), maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, filePath, devices));
			simInfoSW.push_back(ShearWaveSimulationInfo::getNewInstance(u0, v0, l0, lx.at(i), viscosity, kernelName, "ShearWave"));
		}

	std::vector< std::shared_ptr< TestSimulation>> testSimSW = buildTestSimulation(simParaSW, simInfoSW);

	std::vector< std::shared_ptr< TestLogFileInformation>> testLogFileInfo;

	if (nuAndPhiTestSW) {
		std::vector< std::shared_ptr< PhiAndNuTest>> phiAndNuTests = makePhiAndNuTests(testSimSW, simInfoSW, viscosity);
		std::shared_ptr< PhiAndNuInformation> phiNuLogFileInfo = PhiAndNuInformation::getNewInstance(phiAndNuTests);
		testLogFileInfo.push_back(phiNuLogFileInfo);
	}
		

	for (int i = 0; i < testSimSW.size(); i++)
		testSimulation.push_back(testSimSW.at(i));

	std::shared_ptr< ShearWaveInformation> swInfo = ShearWaveInformation::getNewInstance(u0, v0, sw, lx, l0);
	std::shared_ptr< LogFileTimeInformation> logFileTimeInfo = LogFileTimeInformation::getNewInstance(testSimSW, writeFiles);
	makeLogFileWriter(testLogFileInfo, logFileTimeInfo, swInfo, kernelName, viscosity);
}

std::vector< std::shared_ptr< PhiAndNuTest>> ConfigFileReader::makePhiAndNuTests(std::vector< std::shared_ptr< TestSimulation>> testSim, std::vector< std::shared_ptr< SimulationInfo>> simInfo, double viscosity)
{
	std::vector< std::shared_ptr< PhiAndNuTest>> phiAndNuTests;

	for (int i = 1; i < testSim.size(); i++) {
		for (int j = 0; j < i; j++) {
			std::shared_ptr< PhiAndNuTest> test = PhiAndNuTest::getNewInstance(colorOutput, dataToCalcPhiAndNuTest, minOrderOfAccuracy, viscosity, startStepCalculationPhiNu, endStepCalculationPhiNu);
			test->addSimulation(testSim.at(j), simInfo.at(j));
			test->addSimulation(testSim.at(i), simInfo.at(i));

			testSim.at(j)->registerSimulationObserver(test);
			testSim.at(i)->registerSimulationObserver(test);

			phiAndNuTests.push_back(test);
			testQueue->addTest(test);
		}
	}
	return phiAndNuTests;
}

std::vector<std::shared_ptr<L2NormTest>> ConfigFileReader::makeL2NormTests(std::vector<std::shared_ptr<TestSimulation>> testSim, std::vector< std::shared_ptr< SimulationInfo>> simInfo, std::vector<std::shared_ptr< AnalyticalResults>> analyticalResults)
{
	std::vector<std::shared_ptr<L2NormTest>> l2Tests;
	for (int i = 0; i < testSim.size(); i++) {
		std::shared_ptr<L2NormTest> test = L2NormTest::getNewInstance(analyticalResults.at(i), colorOutput);
		test->addSimulation(testSim.at(i), simInfo.at(i));
		testSim.at(i)->registerSimulationObserver(test);
		l2Tests.push_back(test);
		testQueue->addTest(test);
	}

	return l2Tests;
}

bool ConfigFileReader::shouldSimulationGroupRun(std::vector<bool> test)
{
	for (int i = 0; i < test.size(); i++) {
		if (test.at(i))
			return true;
	}
	return false;
}

unsigned int ConfigFileReader::calcStartStepForToVectorWriter()
{
	std::vector< unsigned int> startStepsTests;
	startStepsTests.push_back(basicTimeStepL2Norm);
	startStepsTests.push_back(startStepCalculationPhiNu);

	std::sort(startStepsTests.begin(), startStepsTests.end());

	return startStepsTests.at(0);
}

void ConfigFileReader::makeLogFileWriter(std::vector< std::shared_ptr< TestLogFileInformation>> testLogFiles, std::shared_ptr< LogFileTimeInformation> logFileTimeInfo, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::string kernelName, double viscosity)
{
	std::shared_ptr< LogFileWriter> logFileWriter = LogFileWriter::getNewInstance(testLogFiles, logFileTimeInfo, simLogInfo, kernelName, viscosity, devices, numberOfTimeSteps, basisTimeStepLength, calcStartStepForToVectorWriter());
	logFileWriterQueue->addLogFileWriter(logFileWriter);
}

std::vector<std::shared_ptr<TestSimulation>> ConfigFileReader::getTestSimulations()
{
	return testSimulation;
}

std::shared_ptr<TestQueue> ConfigFileReader::getTestQueue()
{
	return testQueue;
}

std::shared_ptr<LogFileQueue> ConfigFileReader::getLogFileQueue()
{
	return logFileWriterQueue;
}
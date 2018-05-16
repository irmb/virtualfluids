#include "reader.h"

#include <fstream>
#include <iostream>

#include "utilities/input/Input.h"
#include "utilities/StringUtil/StringUtil.h"

#include "Utilities/TestInformation/TestInformationImp.h"
#include "Utilities/TestParameter/TaylorGreenTestParameter/TaylorGreenTestParameter.h"
#include "Utilities/TestParameter/ShearWaveTestParameter/ShearWaveTestParameter.h"

#include "Tests/PhiAndNuTest/PhiAndNuTest.h"
#include "Utilities/TestCout/TestCoutImp.h"
#include "Utilities/SimulationInfo/SimulationInfoImp.h"

#include "Utilities/LogFileInformation/LogFileInformation.h"
#include "Utilities/LogFileInformation/LogFileInformationOutput/LogFileInformationOutput.h"
#include "Utilities/LogFileInformation/BasicSimulationInfo/BasicSimulationInfo.h"
#include "Utilities/LogFileInformation/TaylorGreenInformation/TaylorGreenInformation.h"
#include "Utilities/LogFileInformation/ShearWaveInformation/ShearWaveInformation.h"
#include "Utilities/LogFileInformation/SimulationTimeInformation/SimulationTimeInformation.h"

bool Reader::testShouldRun(std::vector<bool> test)
{
	for (int i = 0; i < test.size(); i++) {
		if (test.at(i))
			return true;
	}
	return false;
}

std::shared_ptr<Reader> Reader::getNewInstance(const std::string aFilePath)
{
	return std::shared_ptr<Reader>(new Reader(aFilePath));
}

std::shared_ptr<TestInformation> Reader::getTestInformation()
{
	return testInfo;
}

std::vector<std::shared_ptr<TestParameter>> Reader::getTestParameter()
{
	return testParameter;
}

Reader::Reader(const std::string aFilePath)
{
	std::ifstream stream;
	stream.open(aFilePath.c_str(), std::ios::in);
	std::unique_ptr<input::Input> input = input::Input::makeInput(stream, "config");

	devices = StringUtil::toVector(input->getValue("Devices"));

	viscosity = StringUtil::toDouble(input->getValue("Viscosity"));
	minOrderOfAccuracy = StringUtil::toDouble(input->getValue("MinOrderOfAccuracy"));

	amplitudeTGV = StringUtil::toDouble(input->getValue("Amplitude_TGV"));
	u0TGV = StringUtil::toDouble(input->getValue("u0_TGV"));
	v0SW = StringUtil::toDouble(input->getValue("v0_SW"));
	u0SW = StringUtil::toDouble(input->getValue("u0_SW"));
	
	numberOfTimeSteps = StringUtil::toInt(input->getValue("NumberOfTimeSteps"));
	basisTimeStepLength = StringUtil::toInt(input->getValue("BasisTimeStepLength"));
	startStepCalculation = StringUtil::toInt(input->getValue("StartStepCalculation"));

	grids.resize(5);
	grids.at(0) = input->getValue("GridPath32");
	grids.at(1) = input->getValue("GridPath64");
	grids.at(2) = input->getValue("GridPath128");
	grids.at(3) = input->getValue("GridPath256");
	grids.at(4) = input->getValue("GridPath512");

	ySliceForCalculation = StringUtil::toInt(input->getValue("ySliceForCalculation"));
	
	writeFiles = StringUtil::toBool(input->getValue("WriteFiles"));
	filePath = input->getValue("PathForFileWriting");
	startStepFileWriter = StringUtil::toInt(input->getValue("StartStepFileWriter"));
	logFilePath= input->getValue("PathLogFile");;

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

	l.resize(5);
	l.at(0) = 32.0;
	l.at(1) = 64.0;
	l.at(2) = 128.0;
	l.at(3) = 256.0;
	l.at(4) = 512.0;

	stream.close();

	testResults.resize(0);
	simInfo.resize(0);
	testOutput = TestCoutImp::getNewInstance();
	makeTestParameter();
	makeTestInformation();
}

void Reader::makeTestInformation()
{
	testInfo = TestInformationImp::getNewInstance();

	makeSimulationInfo();
	testInfo->setSimulationInfo(simInfo);

	testInfo->setLogFilePath(logFilePath);
	makeLogFileInformation();
	testInfo->setLogFileInformation(logInfo);
}

void Reader::makeTestParameter()
{
	if (testShouldRun(tgv)) {
		std::shared_ptr< PhiAndNuTest> tgvTestResults = PhiAndNuTest::getNewInstance("TaylorGreenVortex", minOrderOfAccuracy, testOutput);
		testResults.push_back(tgvTestResults);
		for (int i = 0; i < tgv.size(); i++) {
			if (tgv.at(i)) {
				testParameter.push_back(TaylorGreenTestParameter::getNewInstance(u0TGV, amplitudeTGV, viscosity, l.at(i), numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, grids.at(i), writeFiles, startStepFileWriter, filePath, tgvTestResults, devices));
			}
		}
	}

	if (testShouldRun(sw)) {
		std::shared_ptr< PhiAndNuTest> swTestResults = PhiAndNuTest::getNewInstance("ShearWave", minOrderOfAccuracy, testOutput);
		testResults.push_back(swTestResults);

		for (int i = 0; i < sw.size(); i++) {
			if (sw.at(i)) {
				testParameter.push_back(ShearWaveTestParameter::getNewInstance(u0SW, v0SW, viscosity, l.at(i), numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, grids.at(i), writeFiles, startStepFileWriter, filePath, swTestResults, devices));
			}
		}
	}
}

void Reader::makeSimulationInfo()
{
	for (int i = 0; i < tgv.size(); i++) {
		if (tgv.at(i)) {
			simInfo.push_back(SimulationInfoImp::getNewInstance(testOutput, "TaylorGreenVortex", l.at(i)));
		}
	}
	for (int i = 0; i < sw.size(); i++) {
		if (sw.at(i)) {
			simInfo.push_back(SimulationInfoImp::getNewInstance(testOutput, "ShearWave", l.at(i)));
		}
	}
}

void Reader::makeLogFileInformation()
{
	
	logInfo.push_back(LogFileInformationOutput::getNewInstance(devices));
	logInfo.push_back(BasicSimulationInfo::getNewInstance(numberOfTimeSteps, basisTimeStepLength, startStepCalculation, viscosity));

	if (testShouldRun(tgv))
		logInfo.push_back(TaylorGreenInformation::getNewInstance(u0TGV, amplitudeTGV));
	if (testShouldRun(sw))
		logInfo.push_back(ShearWaveInformation::getNewInstance(u0SW, v0SW));

	logInfo.push_back(SimulationTimeInformation::getNewInstance(simInfo, writeFiles));

	for (int i = 0; i < testResults.size(); i++) {
		logInfo.push_back(testResults.at(i));
	}
}

std::vector< std::shared_ptr< TestResults> > Reader::getAllTestResults()
{
	std::vector< std::shared_ptr< TestResults> > allTestResults;
	for (int i = 0; i < testResults.size(); i++) {
		allTestResults.push_back(testResults.at(i));
	}
	return allTestResults;
}

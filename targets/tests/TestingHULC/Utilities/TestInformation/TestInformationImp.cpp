#include "TestInformationImp.h"

#include "Utilities/SimulationInfo/SimulationInfo.h"
#include "Utilities/LogFileWriter/LogFileWriter.h"
#include "Utilities/LogFileInformation/LogFileInformation.h"
#include "Utilities/TestResults/TestResults.h"
#include "Utilities/TestCout/TestCout.h"

#include <iomanip>

std::shared_ptr<TestInformationImp> TestInformationImp::getNewInstance()
{
	return std::shared_ptr<TestInformationImp>(new TestInformationImp());
}

void TestInformationImp::setSimulationInfo(std::vector<std::shared_ptr<SimulationInfo>> simInfo)
{
	this->simInfos = simInfo;
}

void TestInformationImp::makeSimulationHeadOutput(int i)
{
	simInfos.at(i)->makeSimulationHeadOutput();
}

void TestInformationImp::setSimulationStartTime(int i)
{
	simInfos.at(i)->setStartTime();
}

void TestInformationImp::setSimulationEndTime(int i)
{
	simInfos.at(i)->setEndTime();
}

void TestInformationImp::writeLogFile()
{
	std::shared_ptr<LogFileWriter> logFile = LogFileWriter::getNewInstance(logFilePath);

	for (int i = 0; i < logInfos.size(); i++) {
		logFile->makeOutput(logInfos.at(i)->getOutput());
	}
}

void TestInformationImp::makeFinalTestOutput()
{
	colorOutput->makeFinalTestOutput(getNumberOfPassedTests(),getNumberOfTests());
	for (int i = 0; i < testResults.size(); i++) {
		testResults.at(i)->makeFinalOutput();
	}
	colorOutput->makeFinalTestOutput(getNumberOfPassedTests(), getNumberOfTests());
}

void TestInformationImp::setLogFilePath(std::string aLogFilePath)
{
	this->logFilePath = aLogFilePath;
}

void TestInformationImp::setLogFileInformation(std::vector<std::shared_ptr<LogFileInformation>> logInfo)
{
	this->logInfos = logInfo;
}

void TestInformationImp::setTestResults(std::vector<std::shared_ptr<TestResults>> testResults)
{
	this->testResults = testResults;
}

void TestInformationImp::setColorOutput(std::shared_ptr<TestCout> colorOutput)
{
	this->colorOutput = colorOutput;
}

void TestInformationImp::makeHashLine()
{
	oss << "#################################################" << std::endl;
}

void TestInformationImp::makeCenterHead(std::string output)
{
	makeHashLine();
	oss << "#" << std::setfill(' ') << std::right << std::setw(24 + output.length() / 2) << output << std::setw(24 - output.length() / 2) << "#" << std::endl;
	makeHashLine();
}

int TestInformationImp::getNumberOfTests()
{
	int numberOfTests = 0;
	for (int i = 0; i < testResults.size(); i++) {
		numberOfTests += testResults.at(i)->getNumberOfTests();
	}
	return numberOfTests;
}

int TestInformationImp::getNumberOfPassedTests()
{
	int numberOfPassedTests = 0;
	for (int i = 0; i < testResults.size(); i++) {
		numberOfPassedTests += testResults.at(i)->getNumberOfPassedTests();
	}
	return numberOfPassedTests;
}

TestInformationImp::TestInformationImp()
{
	simInfos.resize(0);
	logInfos.resize(0);
}
#include "EvaluationParameter.h"

std::shared_ptr<EvaluationParameter> EvaluationParameter::getNewInstance(std::string testNames, int numberOfEqualTests, int lx, std::string dataToCalculate, std::string logFilePath, double minOrderOfAccuracy, bool writeFiles, double viscosity)
{
	return std::shared_ptr<EvaluationParameter>(new EvaluationParameter(testNames, numberOfEqualTests, lx, dataToCalculate, logFilePath, minOrderOfAccuracy, writeFiles, viscosity));
}

EvaluationParameter::EvaluationParameter(std::string testName, int numberOfEqualTests, int lx, std::string dataToCalculate, std::string logFilePath, double minOrderOfAccuracy, bool writeFiles, double viscosity)
	: lx(lx),  testName(testName), numberOfEqualTests(numberOfEqualTests), dataToCalculate(dataToCalculate), logFilePath(logFilePath), minOrderOfAccuracy(minOrderOfAccuracy), writeFiles(writeFiles), viscosity(viscosity)
{

}

std::string EvaluationParameter::getTestName()
{
	return testName;
}

std::string EvaluationParameter::getDataToCalculate()
{
	return dataToCalculate;
}

int EvaluationParameter::getLx()
{
	return lx;
}

int EvaluationParameter::getNumberOfEqualTests()
{
	return numberOfEqualTests;
}

void EvaluationParameter::setStartTime()
{
	startTime = time(NULL);
}

void EvaluationParameter::setEndTime()
{
	endTime = time(NULL);
}

double EvaluationParameter::getTestTime()
{
	return difftime(endTime, startTime);
}

std::string EvaluationParameter::getLogFilePath()
{
	return logFilePath;
}

double EvaluationParameter::getMinOrderOfAccuracy()
{
	return minOrderOfAccuracy;
}

bool EvaluationParameter::getWriteFiles()
{
	return writeFiles;
}

double EvaluationParameter::getViscosity()
{
	return viscosity;
}

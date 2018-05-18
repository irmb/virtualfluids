#include "TestParameterImp.h"

#include "Utilities/TestResults/TestResults.h"

TestParameterImp::TestParameterImp(
	real viscosity, real lx, real lz, real l0,
	unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, 
	unsigned int startStepCalculation, unsigned int ySliceForCalculation, 
	std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
	bool writeFiles, unsigned int startStepFileWriter, 
	std::shared_ptr<TestResults> testResults,
	std::vector<int> devices)
		:viscosity(viscosity), lx(lx), l0(l0), lz(lz),
		numberOfTimeSteps(numberOfTimeSteps), basisTimeStepLength(basisTimeStepLength), 
		startStepCalculation(startStepCalculation), ySliceForCalculation(ySliceForCalculation), 
		gridPath(gridPath), maxLevel(maxLevel), numberOfGridLevels(numberOfGridLevels),
		writeFiles(writeFiles), startStepFileWriter(startStepFileWriter), 
		testResults(testResults), devices(devices)
{
	timeStepLength = basisTimeStepLength*(lx / l0)*(lx / l0);
	startTimeCalculation = timeStepLength * startStepCalculation;
	startTimeDataWriter = timeStepLength * startStepFileWriter;
	endTime = timeStepLength * numberOfTimeSteps;

}

double TestParameterImp::getViscosity()
{
	return viscosity;
}

std::string TestParameterImp::getGridPath()
{
	return gridPath;
}

std::string TestParameterImp::getFilePath()
{
	return filePath;
}

unsigned int TestParameterImp::getNumberOfGridLevels()
{
	return numberOfGridLevels;
}

unsigned int TestParameterImp::getEndTime()
{
	return endTime;
}

unsigned int TestParameterImp::getTimeStepLength()
{
	return timeStepLength;
}

unsigned int TestParameterImp::getLx()
{
	return lx;
}

unsigned int TestParameterImp::getLz()
{
	return lz;
}

unsigned int TestParameterImp::getYSliceForCalculation()
{
	return ySliceForCalculation;
}

unsigned int TestParameterImp::getStartTimeCalculation()
{
	return startTimeCalculation;
}

bool TestParameterImp::getWriteFiles()
{
	return writeFiles;
}

unsigned int TestParameterImp::getStartTimeDataWriter()
{
	return startTimeDataWriter;
}

std::vector<int> TestParameterImp::getDevices()
{
	return devices;
}

std::shared_ptr<InitialCondition> TestParameterImp::getInitialCondition()
{
	return initialCondition;
}

std::shared_ptr<Calculator> TestParameterImp::getCalculator()
{
	return calculator;
}

std::shared_ptr<TestResults> TestParameterImp::getTestResults()
{
	return testResults;
}

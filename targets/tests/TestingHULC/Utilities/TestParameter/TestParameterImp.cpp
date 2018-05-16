#include "TestParameterImp.h"

#include "Utilities/TestResults/TestResults.h"

TestParameterImp::TestParameterImp(
	real viscosity, unsigned int lx, 
	unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, 
	unsigned int startStepCalculation, unsigned int ySliceForCalculation, 
	std::string gridPath, 
	bool writeFiles, unsigned int startStepFileWriter, 
	std::shared_ptr<TestResults> testResults,
	std::vector<int> devices)
		:viscosity(viscosity), lx(lx),
		numberOfTimeSteps(numberOfTimeSteps), basisTimeStepLength(basisTimeStepLength), 
		startStepCalculation(startStepCalculation), ySliceForCalculation(ySliceForCalculation), 
		gridPath(gridPath), 
		writeFiles(writeFiles), startStepFileWriter(startStepFileWriter), 
		testResults(testResults), devices(devices)
{
	maxLevel = 0;
	numberOfGridLevels = 1;
	l0 = 32;
	rho0 = 1.0;

	lz = 3 * lx / 2;
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

double TestParameterImp::getVelocity()
{
	return velocity;
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

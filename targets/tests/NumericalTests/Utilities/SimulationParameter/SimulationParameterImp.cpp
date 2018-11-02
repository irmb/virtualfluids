#include "SimulationParameterImp.h"

#include "Utilities/TestResults/TestResults.h"

SimulationParameterImp::SimulationParameterImp(
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

double SimulationParameterImp::getViscosity()
{
	return viscosity;
}

std::string SimulationParameterImp::getGridPath()
{
	return gridPath;
}

std::string SimulationParameterImp::getFilePath()
{
	return filePath;
}

unsigned int SimulationParameterImp::getNumberOfGridLevels()
{
	return numberOfGridLevels;
}

unsigned int SimulationParameterImp::getEndTime()
{
	return endTime;
}

unsigned int SimulationParameterImp::getTimeStepLength()
{
	return timeStepLength;
}

unsigned int SimulationParameterImp::getLx()
{
	return lx;
}

unsigned int SimulationParameterImp::getLz()
{
	return lz;
}

unsigned int SimulationParameterImp::getYSliceForCalculation()
{
	return ySliceForCalculation;
}

unsigned int SimulationParameterImp::getStartTimeCalculation()
{
	return startTimeCalculation;
}

bool SimulationParameterImp::getWriteFiles()
{
	return writeFiles;
}

unsigned int SimulationParameterImp::getStartTimeDataWriter()
{
	return startTimeDataWriter;
}

std::vector<int> SimulationParameterImp::getDevices()
{
	return devices;
}

std::shared_ptr<InitialCondition> SimulationParameterImp::getInitialCondition()
{
	return initialCondition;
}

std::shared_ptr<Calculator> SimulationParameterImp::getCalculator()
{
	return calculator;
}

std::shared_ptr<TestResults> SimulationParameterImp::getTestResults()
{
	return testResults;
}
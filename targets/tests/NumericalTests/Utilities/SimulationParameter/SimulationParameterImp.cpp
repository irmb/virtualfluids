#include "SimulationParameterImp.h"

#include <experimental/filesystem>

SimulationParameterImp::SimulationParameterImp(std::string simName, real viscosity, real lx, real lz, real l0,
	unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, 
	unsigned int startStepCalculation, unsigned int ySliceForCalculation, 
	std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
	bool writeFiles, unsigned int startStepFileWriter,
	std::vector<int> devices)
		:simName(simName), viscosity(viscosity), lx(lx), l0(l0), lz(lz),
		numberOfTimeSteps(numberOfTimeSteps), basisTimeStepLength(basisTimeStepLength), 
		startStepCalculation(startStepCalculation), ySliceForCalculation(ySliceForCalculation), 
		gridPath(gridPath), maxLevel(maxLevel), numberOfGridLevels(numberOfGridLevels),
		writeFiles(writeFiles), startStepFileWriter(startStepFileWriter), devices(devices)
{
	timeStepLength = basisTimeStepLength*(lx / l0)*(lx / l0);
	startTimeCalculation = timeStepLength * startStepCalculation;
	startTimeDataWriter = timeStepLength * startStepFileWriter;
	endTime = timeStepLength * numberOfTimeSteps;

}

void SimulationParameterImp::generateFilePath(std::string filePath)
{
	std::experimental::filesystem::path dir(filePath);
	if (!(std::experimental::filesystem::exists(dir)))
		std::experimental::filesystem::create_directories(dir);
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

std::shared_ptr<KernelConfiguration> SimulationParameterImp::getKernelConfiguration()
{
	return kernelConfig;
}
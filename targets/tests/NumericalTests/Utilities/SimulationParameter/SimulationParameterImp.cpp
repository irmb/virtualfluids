#include "SimulationParameterImp.h"

#include "Utilities\KernelConfiguration\KernelConfigurationImp.h"

#include "Utilities\Structs\BasicSimulationParameterStruct.h"
#include "Utilities\Structs\GridInformationStruct.h"

#include <experimental/filesystem>

SimulationParameterImp::SimulationParameterImp(std::string kernelName, double viscosity, std::shared_ptr<BasicSimulationParameterStruct> basicSimPara, std::shared_ptr<GridInformationStruct> gridInfo)
	: viscosity(viscosity)
{
	kernelConfig = KernelConfigurationImp::getNewInstance(kernelName);

	devices = basicSimPara->devices;
	numberOfTimeSteps = basicSimPara->numberOfTimeSteps;
	
	gridPath = gridInfo->gridPath;
	lx = gridInfo->lx;
	lz = gridInfo->lz;
	maxLevel = gridInfo->maxLevel;
	numberOfGridLevels = gridInfo->numberOfGridLevels;
}

void SimulationParameterImp::generateFileDirectionInMyStystem(std::string filePath)
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
	return timeStepLength * numberOfTimeSteps;
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
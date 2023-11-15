#include "SimulationParameterImp.h"

#include "Utilities/KernelConfiguration/KernelConfigurationImp.h"

#include "Utilities/Structs/BasicSimulationParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#include <filesystem>

SimulationParameterImp::SimulationParameterImp(std::string kernel, double viscosity, std::shared_ptr<BasicSimulationParameterStruct> basicSimPara, std::shared_ptr<GridInformationStruct> gridInfo)
    : viscosity(viscosity)
{
    kernelConfig = KernelConfigurationImp::getNewInstance(kernel);

    devices = basicSimPara->devices;
    numberOfTimeSteps = basicSimPara->numberOfTimeSteps;
    
    gridPath = gridInfo->gridPath;
    lx = gridInfo->lx;
    lz = gridInfo->lz;
    l0 = basicSimPara->l0;
    maxLevel = gridInfo->maxLevel;
    numberOfGridLevels = gridInfo->numberOfGridLevels;
}

void SimulationParameterImp::generateFileDirectionInMyStystem(std::string filePath)
{
    std::filesystem::path dir(filePath);
    if (!(std::filesystem::exists(dir)))
        std::filesystem::create_directories(dir);
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

unsigned int SimulationParameterImp::getL0()
{
    return l0;
}

std::vector<unsigned int> SimulationParameterImp::getDevices()
{
    return devices;
}

double SimulationParameterImp::getMaxVelocity()
{
    return maxVelocity;
}

std::shared_ptr<KernelConfiguration> SimulationParameterImp::getKernelConfiguration()
{
    return kernelConfig;
}
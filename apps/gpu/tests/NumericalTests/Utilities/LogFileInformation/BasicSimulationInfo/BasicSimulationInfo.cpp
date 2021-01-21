#include "BasicSimulationInfo.h"

#include "VirtualFluids_GPU/Kernel/Utilities/Mapper/KernelMapper/KernelMapper.h"

std::shared_ptr<BasicSimulationInfo> BasicSimulationInfo::getNewInstance(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, KernelType kernel)
{
	return std::shared_ptr<BasicSimulationInfo>(new BasicSimulationInfo(numberOfTimeSteps, viscosity, basicTimeStepLength, kernel));
}

std::string BasicSimulationInfo::getOutput()
{
	makeCenterHead("Basic Simulation Information");
	oss << "Kernel=" << kernelName << std::endl;
	oss << "NumberOfTimeSteps=" << numberOfTimeSteps << std::endl;
	oss << "Viscosity=" << viscosity << std::endl;
	oss << "BasisTimeStepLength=" << basicTimeStepLength << std::endl;
	oss << std::endl;
	return oss.str();
}

BasicSimulationInfo::BasicSimulationInfo(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, KernelType kernel)
  : numberOfTimeSteps(numberOfTimeSteps), viscosity(viscosity), basicTimeStepLength(basicTimeStepLength)
{
	std::shared_ptr<KernelMapper> myKernelMapper = KernelMapper::getInstance();
	kernelName = myKernelMapper->getString(kernel);
}
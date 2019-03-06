#include "BasicSimulationInfo.h"

std::shared_ptr<BasicSimulationInfo> BasicSimulationInfo::getNewInstance(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, std::string kernelName)
{
	return std::shared_ptr<BasicSimulationInfo>(new BasicSimulationInfo(numberOfTimeSteps, viscosity, basicTimeStepLength, kernelName));
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

BasicSimulationInfo::BasicSimulationInfo(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, std::string kernelName):numberOfTimeSteps(numberOfTimeSteps), viscosity(viscosity), basicTimeStepLength(basicTimeStepLength), kernelName(kernelName)
{
	
}
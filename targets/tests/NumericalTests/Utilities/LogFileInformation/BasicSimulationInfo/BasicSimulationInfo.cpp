#include "BasicSimulationInfo.h"

std::shared_ptr<BasicSimulationInfo> BasicSimulationInfo::getNewInstance(int numberOfTimeSteps, int basisTimeStepLength, double viscosity)
{
	return std::shared_ptr<BasicSimulationInfo>(new BasicSimulationInfo(numberOfTimeSteps, basisTimeStepLength, viscosity));
}

std::string BasicSimulationInfo::getOutput()
{
	makeCenterHead("Basic Simulation Information");
	oss << "NumberOfTimeSteps: " << numberOfTimeSteps << std::endl;
	oss << "BasisTimeStepLength: " << basisTimeStepLength << std::endl;
	oss << "Viscosity: " << viscosity << std::endl;
	oss << std::endl;
	return oss.str();
}

BasicSimulationInfo::BasicSimulationInfo(int numberOfTimeSteps, int basisTimeStepLength, double viscosity):numberOfTimeSteps(numberOfTimeSteps), basisTimeStepLength(basisTimeStepLength), viscosity(viscosity)
{
	
}
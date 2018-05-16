#include "BasicSimulationInfo.h"

std::shared_ptr<LogFileInformation> BasicSimulationInfo::getNewInstance(int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation, double viscosity)
{
	return std::shared_ptr<LogFileInformation>(new BasicSimulationInfo(numberOfTimeSteps, basisTimeStepLength, startStepCalculation, viscosity));
}

std::string BasicSimulationInfo::getOutput()
{
	makeCenterHead("Basic Information");
	oss << "NumberOfTimeSteps: " << numberOfTimeSteps << std::endl;
	oss << "BasisTimeStepLength: " << basisTimeStepLength << std::endl;
	oss << "StartStepCalculation: " << startStepCalculation << std::endl;
	oss << "Viscosity: " << viscosity << std::endl;
	oss << std::endl;
	return oss.str();
}

BasicSimulationInfo::BasicSimulationInfo(int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation, double viscosity):numberOfTimeSteps(numberOfTimeSteps), basisTimeStepLength(basisTimeStepLength), startStepCalculation(startStepCalculation), viscosity(viscosity)
{
	
}
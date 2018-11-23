#include "SimulationInfoImp.h"

std::string SimulationInfoImp::getKernelName()
{
	return kernelName;
}

double SimulationInfoImp::getViscosity()
{
	return viscosity;
}

std::string SimulationInfoImp::getSimulationName()
{
	return simulationName;
}

std::string SimulationInfoImp::getSimulationParameterString()
{
	return simulationParameterString;
}

int SimulationInfoImp::getLx()
{
	return lx;
}

int SimulationInfoImp::getNumberOfSimulations()
{
	return numberOfSimulations;
}

int SimulationInfoImp::getSimulationID()
{
	return simID;
}

void SimulationInfoImp::setSimulationID(int simID)
{
	this->simID = simID;
}

SimulationInfoImp::SimulationInfoImp(int lx, double viscosity, std::string kernelName, int numberOfSimulations) : lx(lx), viscosity(viscosity), kernelName(kernelName), numberOfSimulations(numberOfSimulations)
{
}

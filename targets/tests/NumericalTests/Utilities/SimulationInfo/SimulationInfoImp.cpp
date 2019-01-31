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

SimulationInfoImp::SimulationInfoImp(int simID, std::string kernelName, double viscosity, int lx, int numberOfSimulations, std::string simulationName)
	: simID(simID), lx(lx), viscosity(viscosity), kernelName(kernelName), numberOfSimulations(numberOfSimulations), simulationName(simulationName)
{
}

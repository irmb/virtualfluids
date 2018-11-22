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

SimulationInfoImp::SimulationInfoImp(int lx, double viscosity, std::string kernelName) : lx(lx), viscosity(viscosity), kernelName(kernelName)
{
}

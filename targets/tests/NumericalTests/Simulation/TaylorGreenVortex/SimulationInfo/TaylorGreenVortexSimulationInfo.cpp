#include "TaylorGreenVortexSimulationInfo.h"

#include <sstream>

std::shared_ptr<SimulationInfo> TaylorGreenVortexSimulationInfo::getNewInstance(double u0, double amplitude, int l0, int lx, double viscosity, std::string kernelName, std::string simulationName)
{
	return std::shared_ptr<SimulationInfo>(new TaylorGreenVortexSimulationInfo(u0, amplitude, l0, lx, viscosity, kernelName, simulationName));
}

TaylorGreenVortexSimulationInfo::TaylorGreenVortexSimulationInfo(double u0, double amplitude, int l0, int lx, double viscosity, std::string kernelName, std::string simulationName): SimulationInfoImp(lx, viscosity, kernelName, simulationName)
{
	std::ostringstream oss;
	oss << " u0: " << u0 / (lx / l0) << " Amplitude: " << amplitude / (lx / l0);
	this->simulationParameterString = oss.str();

}
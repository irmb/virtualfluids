#include "TaylorGreenVortexSimulationInfo.h"

#include <sstream>

std::shared_ptr<SimulationInfo> TaylorGreenVortexSimulationInfo::getNewInstance(double u0, double amplitude, int l0, int lx, double viscosity, std::string kernelName)
{
	return std::shared_ptr<SimulationInfo>(new TaylorGreenVortexSimulationInfo(u0, amplitude, l0, lx, viscosity, kernelName));
}

TaylorGreenVortexSimulationInfo::TaylorGreenVortexSimulationInfo(double u0, double amplitude, int l0, int lx, double viscosity, std::string kernelName) : SimulationInfoImp(lx, viscosity, kernelName)
{
	std::ostringstream oss;
	oss << " u0: " << u0 / (lx / l0) << " Amplitude: " << amplitude / (lx / l0);
	this->simulationParameterString = oss.str();

	simulationName = "TaylorGreenVortex";
}
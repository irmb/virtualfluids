#include "SimulationInfoTaylorGreenVortexUx.h"

#include <sstream>

std::shared_ptr<SimulationInfo> SimulationInfoTaylorGreenUx::getNewInstance(double ux, double amplitude, int l0, int lx, double viscosity, std::string kernelName, int numberOfSimulations)
{
	return std::shared_ptr<SimulationInfo>(new SimulationInfoTaylorGreenUx(ux, amplitude, l0, lx, viscosity, kernelName, numberOfSimulations));
}

SimulationInfoTaylorGreenUx::SimulationInfoTaylorGreenUx(double ux, double amplitude, int l0, int lx, double viscosity, std::string kernelName, int numberOfSimulations) : SimulationInfoImp(lx, viscosity, kernelName, numberOfSimulations)
{
	std::ostringstream oss;
	oss << " ux: " << ux / (lx / l0) << " Amplitude: " << amplitude / (lx / l0);
	this->simulationParameterString = oss.str();

	simulationName = "TaylorGreenVortex Ux";
}
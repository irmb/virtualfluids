#include "SimulationInfoTaylorGreenVortexUz.h"

#include <sstream>

std::shared_ptr<SimulationInfoTaylorGreenUz> SimulationInfoTaylorGreenUz::getNewInstance(double uz, double amplitude, int l0, int lz, double viscosity, std::string kernelName, int numberOfSimulations)
{
	return std::shared_ptr<SimulationInfoTaylorGreenUz>(new SimulationInfoTaylorGreenUz(uz, amplitude, l0, lz, viscosity, kernelName, numberOfSimulations));
}

SimulationInfoTaylorGreenUz::SimulationInfoTaylorGreenUz(double uz, double amplitude, int l0, int lz, double viscosity, std::string kernelName, int numberOfSimulations) : SimulationInfoImp(lz, viscosity, kernelName, numberOfSimulations)
{
	std::ostringstream oss;
	oss << " uz: " << uz / (lz / l0) << " Amplitude: " << amplitude / (lz / l0);
	this->simulationParameterString = oss.str();

	simulationName = "TaylorGreenVortex Uz";
}
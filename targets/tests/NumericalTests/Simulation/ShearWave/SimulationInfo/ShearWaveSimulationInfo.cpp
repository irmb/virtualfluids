#include "ShearWaveSimulationInfo.h"

#include <sstream>

std::shared_ptr<SimulationInfo> ShearWaveSimulationInfo::getNewInstance(double u0, double v0, int l0, int lx, double viscosity, std::string kernelName)
{
	return std::shared_ptr<SimulationInfo>(new ShearWaveSimulationInfo(u0, v0, l0, lx, viscosity, kernelName));
}

ShearWaveSimulationInfo::ShearWaveSimulationInfo(double u0, double v0, int l0, int lx, double viscosity, std::string kernelName) : SimulationInfoImp(lx, viscosity, kernelName)
{
	std::ostringstream oss;
	oss << " u0: " << u0 / (lx / l0) << " v0: " << v0 / (lx / l0);
	this->simulationParameterString = oss.str();

	simulationName = "ShearWave";
}
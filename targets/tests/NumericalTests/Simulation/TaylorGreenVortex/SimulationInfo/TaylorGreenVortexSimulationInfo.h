#ifndef TAYLORGREENVORTEX_SIMULATION_INFO_H
#define TAYLORGREENVORTEX_SIMULATION_INFO_H

#include "Utilities\SimulationInfo\SimulationInfoImp.h"

#include <memory>

class TaylorGreenVortexSimulationInfo : public SimulationInfoImp
{
public:
	static std::shared_ptr< SimulationInfo> getNewInstance(double u0, double amplitude, int l0, int lx, double viscosity, std::string kernelName, std::string simulationName);

private:
	TaylorGreenVortexSimulationInfo() {};
	TaylorGreenVortexSimulationInfo(double u0, double amplitude, int l0, int lx, double viscosity, std::string kernelName, std::string simulationName);
};
#endif 
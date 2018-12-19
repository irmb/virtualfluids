#ifndef SIMULATION_INFO_TAYLORGREENVORTEX_UX_H
#define SIMULATION_INFO_TAYLORGREENVORTEX_UX_H

#include "Utilities\SimulationInfo\SimulationInfoImp.h"

#include <memory>

class SimulationInfoTaylorGreenUx : public SimulationInfoImp
{
public:
	static std::shared_ptr< SimulationInfo> getNewInstance(double ux, double amplitude, int l0, int lx, double viscosity, std::string kernelName, int numberOfSimulations);

private:
	SimulationInfoTaylorGreenUx() {};
	SimulationInfoTaylorGreenUx(double ux, double amplitude, int l0, int lx, double viscosity, std::string kernelName, int numberOfSimulations);
	
};
#endif 
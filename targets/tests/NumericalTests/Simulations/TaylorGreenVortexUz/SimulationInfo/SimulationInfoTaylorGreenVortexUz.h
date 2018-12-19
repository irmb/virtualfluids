#ifndef SIMULATION_INFO_TAYLORGREENVORTEX_UZ_H
#define SIMULATION_INFO_TAYLORGREENVORTEX_UZ_H

#include "Utilities\SimulationInfo\SimulationInfoImp.h"

#include <memory>

class SimulationInfoTaylorGreenUz : public SimulationInfoImp
{
public:
	static std::shared_ptr< SimulationInfoTaylorGreenUz> getNewInstance(double uz, double amplitude, int l0, int lz, double viscosity, std::string kernelName, int numberOfSimulations);

private:
	SimulationInfoTaylorGreenUz() {};
	SimulationInfoTaylorGreenUz(double uz, double amplitude, int l0, int lz, double viscosity, std::string kernelName, int numberOfSimulations);
	
};
#endif 
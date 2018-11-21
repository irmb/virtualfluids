#ifndef SHEARWAVE_SIMULATION_INFO_H
#define SHEARWAVE_SIMULATION_INFO_H

#include "Utilities\SimulationInfo\SimulationInfoImp.h"

#include <memory>

class ShearWaveSimulationInfo : public SimulationInfoImp
{
public:
	static std::shared_ptr< SimulationInfo> getNewInstance(double u0, double v0, int l0, int lx, double viscosity, std::string kernelName, std::string simulationName);

private:
	ShearWaveSimulationInfo() {};
	ShearWaveSimulationInfo(double u0, double v0, int l0, int lx, double viscosity, std::string kernelName, std::string simulationName);
};
#endif 
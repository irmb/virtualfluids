#ifndef SHEAR_WAVE_SIMULATION_PARAMETER_H
#define SHEAR_WAVE_SIMULATION_PARAMETER_H

#include "Utilities/SimulationParameter/SimulationParameterImp.h"

struct ShearWaveParameterStruct;

class ShearWaveSimulationParameter : public SimulationParameterImp
{
public:
	static std::shared_ptr<SimulationParameter> getNewInstance(std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo);
	double getMaxVelocity();

protected:
	ShearWaveSimulationParameter() {};
	ShearWaveSimulationParameter(std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo);

private:
	real ux, uz, rho0;
};

#endif

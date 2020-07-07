#ifndef SHEAR_WAVE_SIMULATION_PARAMETER_H
#define SHEAR_WAVE_SIMULATION_PARAMETER_H

#include "Utilities/SimulationParameter/SimulationParameterImp.h"

struct ShearWaveParameterStruct;

class ShearWaveSimulationParameter : public SimulationParameterImp
{
public:
	static std::shared_ptr<SimulationParameter> getNewInstance(KernelType kernel, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo);

protected:
	ShearWaveSimulationParameter() {};
	ShearWaveSimulationParameter(KernelType kernel, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo);

};

#endif

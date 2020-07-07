#ifndef SIMULATION_PARAMETER_TaylorGreen_Ux_H
#define SIMULATION_PARAMETER_TaylorGreen_Ux_H

#include "Utilities/SimulationParameter/SimulationParameterImp.h"

#include <string>
#include <memory>

struct TaylorGreenVortexUxParameterStruct;

class SimulationParameterTaylorGreenUx : public SimulationParameterImp
{
public:
	static std::shared_ptr<SimulationParameter> getNewInstance(KernelType kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo);

	
protected:
	SimulationParameterTaylorGreenUx(KernelType kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo);

};
#endif 

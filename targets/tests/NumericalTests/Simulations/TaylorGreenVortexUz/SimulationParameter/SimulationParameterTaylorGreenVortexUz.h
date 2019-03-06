#ifndef SIMULATION_PARAMETER_TAYLORGREENVORTEX_UZ_H
#define SIMULATION_PARAMETER_TAYLORGREENVORTEX_UZ_H

#include "Utilities/SimulationParameter/SimulationParameterImp.h"

#include <string>
#include <memory>

struct TaylorGreenVortexUzParameterStruct;
struct GridInformationStruct;

class SimulationParameterTaylorGreenUz : public SimulationParameterImp
{
public:
	static std::shared_ptr<SimulationParameterTaylorGreenUz> getNewInstance(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo);
	
protected:
	SimulationParameterTaylorGreenUz(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo);

};
#endif 

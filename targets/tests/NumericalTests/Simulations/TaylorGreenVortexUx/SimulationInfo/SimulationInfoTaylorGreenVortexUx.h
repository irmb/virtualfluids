#ifndef SIMULATION_INFO_TAYLORGREENVORTEX_UX_H
#define SIMULATION_INFO_TAYLORGREENVORTEX_UX_H

#include "Utilities\SimulationInfo\SimulationInfoImp.h"

#include <memory>

class TimeInfo;

struct TaylorGreenVortexUxParameterStruct;
struct GridInformationStruct;

class SimulationInfoTaylorGreenUx : public SimulationInfoImp
{
public:
	static std::shared_ptr<SimulationInfoTaylorGreenUx> getNewInstance(int simID, std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations);

private:
	SimulationInfoTaylorGreenUx() {};
	SimulationInfoTaylorGreenUx(int simID, std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations);
	
};
#endif 
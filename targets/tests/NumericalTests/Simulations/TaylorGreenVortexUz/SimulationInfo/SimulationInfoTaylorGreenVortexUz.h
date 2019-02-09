#ifndef SIMULATION_INFO_TAYLORGREENVORTEX_UZ_H
#define SIMULATION_INFO_TAYLORGREENVORTEX_UZ_H

#include "Utilities/SimulationInfo/SimulationInfoImp.h"

#include <memory>

class TimeInfo;

struct TaylorGreenVortexUzParameterStruct;
struct GridInformationStruct;

class SimulationInfoTaylorGreenUz : public SimulationInfoImp
{
public:
	static std::shared_ptr<SimulationInfoTaylorGreenUz> getNewInstance(int simID, std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations);

private:
	SimulationInfoTaylorGreenUz() {};
	SimulationInfoTaylorGreenUz(int simID, std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations);
	
};
#endif 
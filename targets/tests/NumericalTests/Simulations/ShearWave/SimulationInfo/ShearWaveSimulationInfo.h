#ifndef SHEARWAVE_SIMULATION_INFO_H
#define SHEARWAVE_SIMULATION_INFO_H

#include "Utilities\SimulationInfo\SimulationInfoImp.h"

#include <memory>

class TimeInfo;

struct ShearWaveParameterStruct;
struct GridInformationStruct;

class ShearWaveSimulationInfo : public SimulationInfoImp
{
public:
	static std::shared_ptr<ShearWaveSimulationInfo> getNewInstance(int simID, std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations);

private:
	ShearWaveSimulationInfo() {};
	ShearWaveSimulationInfo(int simID, std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations);
};
#endif 
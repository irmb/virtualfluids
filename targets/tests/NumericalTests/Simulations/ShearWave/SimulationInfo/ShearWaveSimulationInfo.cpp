#include "ShearWaveSimulationInfo.h"

#include "Simulations\ShearWave\ShearWaveParameterStruct.h"
#include "Utilities\Structs\GridInformationStruct.h"

#include <sstream>

std::shared_ptr<SimulationInfo> ShearWaveSimulationInfo::getNewInstance(int simID, std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations)
{
	return std::shared_ptr<SimulationInfo>(new ShearWaveSimulationInfo(simID, kernelName,viscosity, simParaStruct, gridInfoStruct, numberOfSimulations));
}

ShearWaveSimulationInfo::ShearWaveSimulationInfo(int simID, std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations)
	: SimulationInfoImp(simID, kernelName, viscosity, gridInfoStruct->lx, numberOfSimulations, "ShearWave")
{
	std::ostringstream oss;
	oss << " ux: " << simParaStruct->ux / (gridInfoStruct->lx / simParaStruct->l0) << " uz: " << simParaStruct->uz / (gridInfoStruct->lx / simParaStruct->l0);
	this->simulationParameterString = oss.str();

}
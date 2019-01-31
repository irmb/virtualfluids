#include "SimulationInfoTaylorGreenVortexUx.h"

#include "Simulations\TaylorGreenVortexUx\TaylorGreenVortexUxParameterStruct.h"
#include "Utilities\Structs\GridInformationStruct.h"

#include <sstream>

std::shared_ptr<SimulationInfo> SimulationInfoTaylorGreenUx::getNewInstance(int simID, std::string kernelName, double viscosity, std::shared_ptr< TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr< GridInformationStruct> gridInfoStruct, int numberOfSimulations)
{
	return std::shared_ptr<SimulationInfo>(new SimulationInfoTaylorGreenUx(simID, kernelName, viscosity, simParaStruct, gridInfoStruct, numberOfSimulations));
}

SimulationInfoTaylorGreenUx::SimulationInfoTaylorGreenUx(int simID, std::string kernelName, double viscosity, std::shared_ptr< TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr< GridInformationStruct> gridInfoStruct, int numberOfSimulations)
	: SimulationInfoImp(simID, kernelName, viscosity, gridInfoStruct->lx, numberOfSimulations, "TaylorGreenVortex Ux")
{
	std::ostringstream oss;
	oss << " ux: " << simParaStruct->ux / (gridInfoStruct->lx / simParaStruct->l0) << " Amplitude: " << simParaStruct->amplitude / (gridInfoStruct->lx / simParaStruct->l0);
	this->simulationParameterString = oss.str();
}
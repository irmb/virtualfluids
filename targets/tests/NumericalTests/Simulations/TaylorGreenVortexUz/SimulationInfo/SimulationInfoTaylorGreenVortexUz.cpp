#include "SimulationInfoTaylorGreenVortexUz.h"

#include "Simulations/TaylorGreenVortexUz/TaylorGreenVortexUzParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#include <sstream>

std::shared_ptr<SimulationInfoTaylorGreenUz> SimulationInfoTaylorGreenUz::getNewInstance(int simID, std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations)
{
	return std::shared_ptr<SimulationInfoTaylorGreenUz>(new SimulationInfoTaylorGreenUz(simID, kernelName, viscosity, simParaStruct, gridInfoStruct, numberOfSimulations));
}

SimulationInfoTaylorGreenUz::SimulationInfoTaylorGreenUz(int simID, std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations)
	: SimulationInfoImp(simID, kernelName, viscosity, gridInfoStruct->lx, numberOfSimulations, "TaylorGreenVortex Uz", simParaStruct->dataToCalcTests)
{
	std::ostringstream oss;
	oss << " uz: " << simParaStruct->uz / (gridInfoStruct->lz / simParaStruct->l0) << " Amplitude: " << simParaStruct->amplitude / (gridInfoStruct->lz / simParaStruct->l0);
	this->simulationParameterString = oss.str();
}
#include "SimulationInfoTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUx/TaylorGreenVortexUxParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#include <sstream>

std::shared_ptr<SimulationInfoTaylorGreenUx> SimulationInfoTaylorGreenUx::getNewInstance(int simID, KernelType kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations)
{
	return std::shared_ptr<SimulationInfoTaylorGreenUx>(new SimulationInfoTaylorGreenUx(simID, kernel, viscosity, simParaStruct, gridInfoStruct, numberOfSimulations));
}

SimulationInfoTaylorGreenUx::SimulationInfoTaylorGreenUx(int simID, KernelType kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations)
	: SimulationInfoImp(simID, kernel, viscosity, gridInfoStruct->lx, numberOfSimulations, "TaylorGreenVortex Ux", simParaStruct->dataToCalcTests)
{
	std::ostringstream oss;
	oss << " ux: " << simParaStruct->ux / (gridInfoStruct->lx / simParaStruct->l0) << " Amplitude: " << simParaStruct->amplitude / (gridInfoStruct->lx / simParaStruct->l0);
	this->simulationParameterString = oss.str();
}
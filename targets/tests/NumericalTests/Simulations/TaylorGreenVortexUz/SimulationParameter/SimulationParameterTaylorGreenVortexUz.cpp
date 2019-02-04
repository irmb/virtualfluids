#include "SimulationParameterTaylorGreenVortexUz.h"

#include "Simulations/TaylorGreenVortexUz/InitialConditions/InitialConditionTaylorGreenVortexUz.h"
#include "Simulations\TaylorGreenVortexUz\TaylorGreenVortexUzParameterStruct.h"

#include <sstream>

std::shared_ptr<SimulationParameterTaylorGreenUz> SimulationParameterTaylorGreenUz::getNewInstance(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
{
	return std::shared_ptr<SimulationParameterTaylorGreenUz>(new SimulationParameterTaylorGreenUz(kernelName, viscosity, tgvParameterStruct, gridInfo));
}

double SimulationParameterTaylorGreenUz::getMaxVelocity()
{
	return uz / (lz / l0);
}

SimulationParameterTaylorGreenUz::SimulationParameterTaylorGreenUz(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
:SimulationParameterImp(kernelName, viscosity, tgvParameterStruct->basicSimulationParameter, gridInfo)
{
	this->uz = tgvParameterStruct->uz;
	this->amplitude = tgvParameterStruct->amplitude;
	this->l0 = tgvParameterStruct->l0;
	this->timeStepLength = tgvParameterStruct->basicTimeStepLength * (gridInfo->lz / l0)*(gridInfo->lz / l0);
	this->rho0 = tgvParameterStruct->rho0;

	std::ostringstream oss;
	oss << tgvParameterStruct->vtkFilePath << "\\TaylorGreenVortex Uz\\viscosity_" << viscosity << "\\uz_" << uz << "_amplitude_" << amplitude << "\\" << kernelName << "\\grid" << lx;
	generateFileDirectionInMyStystem(oss.str());
	this->filePath = oss.str();
}
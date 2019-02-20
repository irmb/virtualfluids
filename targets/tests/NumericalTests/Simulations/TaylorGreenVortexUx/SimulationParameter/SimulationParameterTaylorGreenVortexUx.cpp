#include "SimulationParameterTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUx/TaylorGreenVortexUxParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#include <sstream>

std::shared_ptr<SimulationParameter> SimulationParameterTaylorGreenUx::getNewInstance(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
{
	return std::shared_ptr<SimulationParameter>(new SimulationParameterTaylorGreenUx(kernelName, viscosity, tgvParameterStruct, gridInfo));
}

double SimulationParameterTaylorGreenUx::getMaxVelocity()
{
	return ux / (lx / l0);
}

SimulationParameterTaylorGreenUx::SimulationParameterTaylorGreenUx(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
:SimulationParameterImp(kernelName, viscosity, tgvParameterStruct->basicSimulationParameter, gridInfo)
{
	this->ux = tgvParameterStruct->ux;
	this->amplitude = tgvParameterStruct->amplitude;
	this->l0 = tgvParameterStruct->l0;
	this->timeStepLength = tgvParameterStruct->basicTimeStepLength * (gridInfo->lx / l0)*(gridInfo->lx / l0);
	this->rho0 = tgvParameterStruct->rho0;

	std::ostringstream oss;
	oss << tgvParameterStruct->vtkFilePath << "/TaylorGreenVortex Ux/Viscosity_" << viscosity << "/ux_" << ux << "_amplitude_" << amplitude << "/" << kernelName << "/grid" << lx;
	generateFileDirectionInMyStystem(oss.str());
	this->filePath = oss.str();
}

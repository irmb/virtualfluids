#include "ShearWaveSimulationParameter.h"

#include "Simulations/ShearWave/InitialConditions/InitialConditionShearWave.h"
#include "Simulations/ShearWave/ShearWaveParameterStruct.h"

#include "Utilities/Structs/GridInformationStruct.h"

#include <sstream>

std::shared_ptr<SimulationParameter> ShearWaveSimulationParameter::getNewInstance(std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
{
	return std::shared_ptr<SimulationParameter>(new ShearWaveSimulationParameter(kernelName, viscosity, parameterStruct, gridInfo));
}

double ShearWaveSimulationParameter::getMaxVelocity()
{
	if(ux > uz)
		return ux / (lx / l0);
	return uz / (lx / l0);
}

ShearWaveSimulationParameter::ShearWaveSimulationParameter(std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
:SimulationParameterImp(kernelName, viscosity, parameterStruct->basicSimulationParameter, gridInfo)
{
	this->ux = parameterStruct->ux;
	this->uz = parameterStruct->uz;
	this->l0 = parameterStruct->l0;
	this->timeStepLength = parameterStruct->basicTimeStepLength * (gridInfo->lx / l0)*(gridInfo->lx / l0);
	this->rho0 = parameterStruct->rho0;

	std::ostringstream oss;
	oss << parameterStruct->vtkFilePath << "/ShearWave/Viscosity_" << viscosity << "/ux_" << ux << "_uz_" << uz << "/" << kernelName << "/grid" << lx;
	generateFileDirectionInMyStystem(oss.str());
	this->filePath = oss.str();
}
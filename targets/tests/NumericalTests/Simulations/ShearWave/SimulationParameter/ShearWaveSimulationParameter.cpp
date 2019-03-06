#include "ShearWaveSimulationParameter.h"

#include "Simulations/ShearWave/InitialConditions/InitialConditionShearWave.h"
#include "Simulations/ShearWave/ShearWaveParameterStruct.h"

#include "Utilities/Structs/GridInformationStruct.h"

#include <sstream>

std::shared_ptr<SimulationParameter> ShearWaveSimulationParameter::getNewInstance(std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
{
	return std::shared_ptr<SimulationParameter>(new ShearWaveSimulationParameter(kernelName, viscosity, parameterStruct, gridInfo));
}

ShearWaveSimulationParameter::ShearWaveSimulationParameter(std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
:SimulationParameterImp(kernelName, viscosity, parameterStruct->basicSimulationParameter, gridInfo)
{
	this->timeStepLength = parameterStruct->basicTimeStepLength * (gridInfo->lx / l0)*(gridInfo->lx / l0);

	if (parameterStruct->ux > parameterStruct->uz)
		this->maxVelocity = parameterStruct->ux / (lx / l0);
	else
		this->maxVelocity = parameterStruct->uz / (lx / l0);
 
	std::ostringstream oss;
	oss << parameterStruct->vtkFilePath << "/ShearWave/Viscosity_" << viscosity << "/ux_" << parameterStruct->ux << "_uz_" << parameterStruct->uz << "/" << kernelName << "/grid" << lx;
	generateFileDirectionInMyStystem(oss.str());
	this->filePath = oss.str();
}
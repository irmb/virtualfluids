#include "SimulationParameterTaylorGreenVortexUz.h"

#include "Simulations/TaylorGreenVortexUz/InitialConditions/InitialConditionTaylorGreenVortexUz.h"
#include "Simulations/TaylorGreenVortexUz/TaylorGreenVortexUzParameterStruct.h"

#include "VirtualFluids_GPU/Kernel/Utilities/Mapper/KernelMapper/KernelMapper.h"

#include <sstream>

std::shared_ptr<SimulationParameterTaylorGreenUz> SimulationParameterTaylorGreenUz::getNewInstance(KernelType kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
{
	return std::shared_ptr<SimulationParameterTaylorGreenUz>(new SimulationParameterTaylorGreenUz(kernel, viscosity, tgvParameterStruct, gridInfo));
}

SimulationParameterTaylorGreenUz::SimulationParameterTaylorGreenUz(KernelType kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
:SimulationParameterImp(kernel, viscosity, tgvParameterStruct->basicSimulationParameter, gridInfo)
{
	this->timeStepLength = tgvParameterStruct->basicTimeStepLength * (gridInfo->lz / l0)*(gridInfo->lz / l0);
	this->maxVelocity = tgvParameterStruct->uz / (lz / l0);

	std::shared_ptr<KernelMapper> myKernelMapper = KernelMapper::getInstance();
	std::string kernelName = myKernelMapper->getString(kernel);

	std::ostringstream oss;
	oss << tgvParameterStruct->vtkFilePath << "/TaylorGreenVortex Uz/Viscosity_" << viscosity << "/uz_" << tgvParameterStruct->uz << "_amplitude_" << tgvParameterStruct->amplitude << "/" << kernelName << "/grid" << lx;
	generateFileDirectionInMyStystem(oss.str());
	this->filePath = oss.str();
}
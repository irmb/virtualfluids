#include "SimulationParameterTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUx/TaylorGreenVortexUxParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#include <sstream>

std::shared_ptr<SimulationParameter> SimulationParameterTaylorGreenUx::getNewInstance(std::string kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
{
    return std::shared_ptr<SimulationParameter>(new SimulationParameterTaylorGreenUx(kernel, viscosity, tgvParameterStruct, gridInfo));
}

SimulationParameterTaylorGreenUx::SimulationParameterTaylorGreenUx(std::string kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
:SimulationParameterImp(kernel, viscosity, tgvParameterStruct->basicSimulationParameter, gridInfo)
{
    this->maxVelocity = tgvParameterStruct->ux / (lx / l0);
    this->timeStepLength = tgvParameterStruct->basicTimeStepLength * (gridInfo->lx / l0)*(gridInfo->lx / l0);

    std::string kernelName = kernel;

    std::ostringstream oss;
    oss << tgvParameterStruct->vtkFilePath << "/TaylorGreenVortex Ux/Viscosity_" << viscosity << "/ux_" << tgvParameterStruct->ux << "_amplitude_" << tgvParameterStruct->amplitude << "/" << kernelName << "/grid" << lx;
    generateFileDirectionInMyStystem(oss.str());
    this->filePath = oss.str();
}

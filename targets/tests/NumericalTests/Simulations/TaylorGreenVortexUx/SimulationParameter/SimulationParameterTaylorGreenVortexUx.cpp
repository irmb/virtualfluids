#include "SimulationParameterTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUx/InitialConditions/InitialConditionTaylorGreenVortexUx.h"
#include "Utilities\KernelConfiguration\KernelConfigurationImp.h"

#include <sstream>

std::shared_ptr<SimulationParameter> SimulationParameterTaylorGreenUx::getNewInstance(std::string kernelName, real ux, real amplitude, real viscosity, real rho0, real lx, real lz, real l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles, unsigned int startStepFileWriter, std::string filePath, std::vector<int> devices)
{
	return std::shared_ptr<SimulationParameter>(new SimulationParameterTaylorGreenUx(kernelName, ux, amplitude, viscosity, rho0, lx, lz, l0, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, filePath, devices));
}

double SimulationParameterTaylorGreenUx::getMaxVelocity()
{
	return ux / (lx / l0);
}

SimulationParameterTaylorGreenUx::SimulationParameterTaylorGreenUx(std::string kernelName, real ux, real amplitude, real viscosity, real rho0, real lx, real lz, real l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles, unsigned int startStepFileWriter, std::string filePath, std::vector<int> devices)
:SimulationParameterImp("TaylorGreenVortex Ux", viscosity, lx, lz, l0, lx, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, devices), ux(ux), amplitude(amplitude), rho0(rho0)
{
	std::ostringstream oss;
	oss << filePath << "\\TaylorGreenVortex Ux\\" << viscosity << "\\ux_" << ux << "_amplitude_" << amplitude << "\\" << kernelName << "\\grid" << lx;
	generateFilePath(oss.str());
	this->filePath = oss.str();

	initialCondition = InitialConditionTaylorGreenUx::getNewInstance(lx, lz, l0, ux, amplitude, rho0);
	kernelConfig = KernelConfigurationImp::getNewInstance(kernelName);
}

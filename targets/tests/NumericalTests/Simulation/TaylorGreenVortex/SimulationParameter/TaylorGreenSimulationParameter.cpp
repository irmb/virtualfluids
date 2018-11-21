#include "TaylorGreenSimulationParameter.h"

#include "Simulation/TaylorGreenVortex/InitialConditions/InitialConditionTaylorGreenVortex.h"
#include "Utilities\KernelConfiguration\KernelConfigurationImp.h"

#include <sstream>

std::shared_ptr<SimulationParameter> TaylorGreenSimulationParameter::getNewInstance(std::string kernelName, real u0, real amplitude, real viscosity, real rho0, real lx, real lz, real l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles, unsigned int startStepFileWriter, std::string filePath, std::vector<int> devices)
{
	return std::shared_ptr<SimulationParameter>(new TaylorGreenSimulationParameter(kernelName, u0, amplitude, viscosity, rho0, lx, lz, l0, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, filePath, devices));
}

double TaylorGreenSimulationParameter::getMaxVelocity()
{
	return u0 / (lx / l0);
}

TaylorGreenSimulationParameter::TaylorGreenSimulationParameter(std::string kernelName, real u0, real amplitude, real viscosity, real rho0, real lx, real lz, real l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles, unsigned int startStepFileWriter, std::string filePath, std::vector<int> devices)
:SimulationParameterImp("TaylorGreenVortex", viscosity, lx, lz, l0, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, devices), u0(u0), amplitude(amplitude), rho0(rho0)
{
	std::ostringstream oss;
	oss << filePath << "\\" << kernelName << "\\TaylorGreenVortex\\" << viscosity << "\\u0_" << u0 << "_amplitude_" << amplitude << "\\grid" << lx;
	generateFilePath(oss.str());
	this->filePath = oss.str();

	initialCondition = InitialConditionTaylorGreen::getNewInstance(lx, lz, l0, u0, amplitude, rho0);
	kernelConfig = KernelConfigurationImp::getNewInstance(kernelName);
}

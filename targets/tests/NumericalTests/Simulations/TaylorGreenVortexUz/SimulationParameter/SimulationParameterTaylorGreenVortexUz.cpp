#include "SimulationParameterTaylorGreenVortexUz.h"

#include "Simulations/TaylorGreenVortexUz/InitialConditions/InitialConditionTaylorGreenVortexUz.h"
#include "Utilities\KernelConfiguration\KernelConfigurationImp.h"

#include <sstream>

std::shared_ptr<SimulationParameterTaylorGreenUz> SimulationParameterTaylorGreenUz::getNewInstance(std::string kernelName, real uz, real amplitude, real viscosity, real rho0, real lx, real lz, real l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles, unsigned int startStepFileWriter, std::string filePath, std::vector<int> devices)
{
	return std::shared_ptr<SimulationParameterTaylorGreenUz>(new SimulationParameterTaylorGreenUz(kernelName, uz, amplitude, viscosity, rho0, lx, lz, l0, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, filePath, devices));
}

double SimulationParameterTaylorGreenUz::getMaxVelocity()
{
	return uz / (lz / l0);
}

SimulationParameterTaylorGreenUz::SimulationParameterTaylorGreenUz(std::string kernelName, real uz, real amplitude, real viscosity, real rho0, real lx, real lz, real l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles, unsigned int startStepFileWriter, std::string filePath, std::vector<int> devices)
:SimulationParameterImp("TaylorGreenVortex Uz", viscosity, lx, lz, l0, lz, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, devices), uz(uz), amplitude(amplitude), rho0(rho0)
{
	std::ostringstream oss;
	oss << filePath << "\\TaylorGreenVortex Uz\\viscosity_" << viscosity << "\\uz_" << uz << "_amplitude_" << amplitude << "\\" << kernelName << "\\grid" << lx;
	generateFilePath(oss.str());
	this->filePath = oss.str();

	initialCondition = InitialConditionTaylorGreenUz::getNewInstance(lx, lz, l0, uz, amplitude, rho0);
	kernelConfig = KernelConfigurationImp::getNewInstance(kernelName);
}

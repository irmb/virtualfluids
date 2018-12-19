#ifndef SIMULATION_PARAMETER_TAYLORGREENVORTEX_UZ_H
#define SIMULATION_PARAMETER_TAYLORGREENVORTEX_UZ_H

#include "Utilities/SimulationParameter/SimulationParameterImp.h"

#include <string>
#include <memory>


class SimulationParameterTaylorGreenUz : public SimulationParameterImp
{
public:
	static std::shared_ptr<SimulationParameterTaylorGreenUz> getNewInstance(std::string kernelName, real uz, real amplitude, real viscosity, real rho0,
														real lx, real lz, real l0,
														unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
														unsigned int startStepCalculation, unsigned int ySliceForCalculation,
														std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
														bool writeFiles, unsigned int startStepFileWriter, std::string filePath,
														std::vector<int> devices);
	double getMaxVelocity();
	
protected:
	SimulationParameterTaylorGreenUz(std::string kernelName, real uz, real amplitude,
							real viscosity, real rho0, real lx, real lz, real l0,
							unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
							unsigned int startStepCalculation, unsigned int ySliceForCalculation,
							std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
							bool writeFiles, unsigned int startStepFileWriter, std::string filePath, 
							std::vector<int> devices);

private:
	real uz, amplitude, rho0;
};
#endif 

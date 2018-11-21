#ifndef SHEAR_WAVE_SIMULATION_PARAMETER_H
#define SHEAR_WAVE_SIMULATION_PARAMETER_H

#include "Utilities/SimulationParameter/SimulationParameterImp.h"

class  PhiAndNuTest;

class ShearWaveSimulationParameter : public SimulationParameterImp
{
public:
	static std::shared_ptr<SimulationParameter> getNewInstance(std::string kernelName, real u0, real v0, real viscosity, real rho0,
														real lx, real lz, real l0,
														unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
														unsigned int startStepCalculation, unsigned int ySliceForCalculation,
														std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
														bool writeFiles, unsigned int startStepFileWriter, std::string filePath,
														std::vector<int> devices);
	double getMaxVelocity();
	std::string getSimulationParameterString();

protected:
	ShearWaveSimulationParameter() {};
	ShearWaveSimulationParameter(std::string kernelName, real u0, real v0, real viscosity, real rho0,
							real lx, real lz, real l0,
							unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
							unsigned int startStepCalculation, unsigned int ySliceForCalculation,
							std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
							bool writeFiles, unsigned int startStepFileWriter, std::string filePath,
							std::vector<int> devices);

private:
	real u0, v0, rho0;
};

#endif

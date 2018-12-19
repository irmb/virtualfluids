#ifndef SIMULATION_PARAMETER_TaylorGreen_Ux_H
#define SIMULATION_PARAMETER_TaylorGreen_Ux_H

#include "Utilities/SimulationParameter/SimulationParameterImp.h"

#include <string>
#include <memory>

class PhiAndNuTest;

class SimulationParameterTaylorGreenUx : public SimulationParameterImp
{
public:
	static std::shared_ptr<SimulationParameter> getNewInstance(std::string kernelName, real ux, real amplitude, real viscosity, real rho0,
														real lx, real lz, real l0,
														unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
														unsigned int startStepCalculation, unsigned int ySliceForCalculation,
														std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
														bool writeFiles, unsigned int startStepFileWriter, std::string filePath,
														std::vector<int> devices);
	double getMaxVelocity();
	
protected:
	SimulationParameterTaylorGreenUx(std::string kernelName, real ux, real amplitude,
							real viscosity, real rho0, real lx, real lz, real l0,
							unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
							unsigned int startStepCalculation, unsigned int ySliceForCalculation,
							std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
							bool writeFiles, unsigned int startStepFileWriter, std::string filePath, 
							std::vector<int> devices);

private:
	real ux, amplitude, rho0;
};
#endif 

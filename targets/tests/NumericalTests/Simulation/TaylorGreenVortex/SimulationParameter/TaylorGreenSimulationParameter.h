#ifndef TGV_SIMULATION_PARAMETER_H
#define TGV_SIMULATION_PARAMETER_H

#include "Utilities/SimulationParameter/SimulationParameterImp.h"

#include <string>
#include <memory>

class PhiAndNuTest;

class TaylorGreenSimulationParameter : public SimulationParameterImp
{
public:
	static std::shared_ptr<SimulationParameter> getNewInstance(real u0, real amplitude, real viscosity, real rho0,
														real lx, real lz, real l0,
														unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
														unsigned int startStepCalculation, unsigned int ySliceForCalculation,
														std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
														bool writeFiles, unsigned int startStepFileWriter, std::string filePath,
														std::shared_ptr<PhiAndNuTest> testResults, std::vector<int> devices);
	double getMaxVelocity();
	
protected:
	TaylorGreenSimulationParameter(real u0, real amplitude,
							real viscosity, real rho0, real lx, real lz, real l0,
							unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
							unsigned int startStepCalculation, unsigned int ySliceForCalculation,
							std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
							bool writeFiles, unsigned int startStepFileWriter, std::string filePath, 
							std::shared_ptr<PhiAndNuTest> testResults, std::vector<int> devices);

private:
	real u0, amplitude, rho0;
};
#endif 

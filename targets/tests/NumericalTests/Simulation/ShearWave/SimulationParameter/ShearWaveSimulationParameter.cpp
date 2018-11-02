#include "ShearWaveSimulationParameter.h"

#include "Simulation/ShearWave/InitialConditions/InitialConditionShearWave.h"
#include "Utilities/Calculator/FFTCalculator/VzFFTCalculator/VzFFTCalculator.h"
#include "Tests/PhiAndNuTest/PhiAndNuTest.h"
#include "Utilities/SimulationResults/SimulationResults.h"

#include <sstream>

std::shared_ptr<SimulationParameter> ShearWaveSimulationParameter::getNewInstance(	real u0, real v0, real viscosity, real rho0, real lx, real lz, real l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
																		unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles,
																		unsigned int startStepFileWriter, std::string filePath, std::shared_ptr<PhiAndNuTest> testResults,
																		std::vector<int> devices)
{
	return std::shared_ptr<SimulationParameter>(new ShearWaveSimulationParameter(	u0, v0, viscosity, rho0, lx, lz, l0, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels,
																		writeFiles, startStepFileWriter, filePath,testResults, devices));
}

double ShearWaveSimulationParameter::getMaxVelocity()
{
	if(u0 > v0)
		return u0 / (lx / l0);
	return v0 / (lx / l0);
}

ShearWaveSimulationParameter::ShearWaveSimulationParameter(	real u0, real v0, real viscosity, real rho0, real lx, real lz, real l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles, unsigned int startStepFileWriter, std::string filePath, std::shared_ptr<PhiAndNuTest> testResults, std::vector<int> devices)
:SimulationParameterImp(viscosity, lx, lz, l0, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, testResults, devices), u0(u0), v0(v0), rho0(rho0)
{
	std::ostringstream oss;
	oss << filePath + "/ShearWave/grid" << lx;
	this->filePath = oss.str();

	initialCondition = InitialConditionShearWave::getNewInstance(lx, lz, l0, u0, v0, rho0);
	calculator = VzFFTCalculator::getNewInstance(viscosity, testResults);
}
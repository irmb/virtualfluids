#include "TaylorGreenTestParameter.h"

#include "Simulation/TaylorGreenVortex/InitialConditions/InitialConditionTaylorGreenVortex.h"
#include "Utilities/Calculator/FFTCalculator/VxFFTCalculator/VxFFTCalculator.h"
#include "Tests/PhiAndNuTest/PhiAndNuTest.h"
#include "Utilities/SimulationResults/SimulationResults.h"

#include <sstream>

std::shared_ptr<TestParameter> TaylorGreenTestParameter::getNewInstance(real u0, real amplitude, real viscosity, real rho0, unsigned int lx, unsigned int lz, unsigned int l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles, unsigned int startStepFileWriter, std::string filePath, std::shared_ptr<PhiAndNuTest> testResults, std::vector<int> devices)
{
	return std::shared_ptr<TestParameter>(new TaylorGreenTestParameter(u0, amplitude, viscosity, rho0, lx, lz, l0, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, filePath, testResults, devices));
}

double TaylorGreenTestParameter::getMaxVelocity()
{
	return u0 / (lx / l0);
}

TaylorGreenTestParameter::TaylorGreenTestParameter(real u0, real amplitude, real viscosity, real rho0, unsigned int lx, unsigned int lz, unsigned int l0, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels, bool writeFiles, unsigned int startStepFileWriter, std::string filePath, std::shared_ptr<PhiAndNuTest> testResults, std::vector<int> devices)
:TestParameterImp(viscosity, lx, lz, l0, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, maxLevel, numberOfGridLevels, writeFiles, startStepFileWriter, testResults, devices), u0(u0), amplitude(amplitude), rho0(rho0)
{
	std::ostringstream oss;
	oss << filePath << "/TaylorGreenVortex/grid" << lx;
	this->filePath = oss.str();

	initialCondition = InitialConditionTaylorGreen::getNewInstance((double)lx, (double)lz, (double)l0, u0, amplitude, rho0);
	simResults = SimulationResults::getNewInstance(lx, lz, timeStepLength);
	calculator = VxFFTCalculator::getNewInstance(viscosity, testResults);
}

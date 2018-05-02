#include "TaylorGreenTestParameter.h"

#include "Utilities\InitialCondition\TaylorGreenVortex\InitialconditionTaylorGreenVortex.h"

#include <sstream>

std::shared_ptr<TestParameter> TaylorGreenTestParameter::getNewInstance(real u0, real amplitude, real viscosity, unsigned int lx, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, bool writeFiles, unsigned int startStepFileWriter, std::string filePath)
{
	return std::shared_ptr<TestParameter>(new TaylorGreenTestParameter(u0, amplitude,
		viscosity, lx,
		numberOfTimeSteps, basisTimeStepLength,
		startStepCalculation, ySliceForCalculation,
		gridPath,
		writeFiles, startStepFileWriter, filePath));
}

TaylorGreenTestParameter::TaylorGreenTestParameter(real u0, real amplitude, real viscosity, unsigned int lx, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, bool writeFiles, unsigned int startStepFileWriter, std::string filePath)
	:TestParameterImp(viscosity, lx, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, writeFiles, startStepFileWriter), u0(u0), amplitude(amplitude)
{
	std::ostringstream oss;
	oss << filePath + "/TaylorGreenVortex/grid" << lx;
	this->filePath = oss.str();

	initialCondition = std::shared_ptr<InitialConditionTaylorGreen>(new InitialConditionTaylorGreen((double)lx, (double)lz, (double)l0, u0, amplitude, rho0));
}

std::shared_ptr<InitialCondition> TaylorGreenTestParameter::getInitialCondition()
{
	return initialCondition;
}

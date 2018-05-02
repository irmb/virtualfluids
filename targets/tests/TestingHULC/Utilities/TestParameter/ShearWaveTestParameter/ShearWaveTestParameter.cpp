#include "ShearWaveTestParameter.h"

#include "Utilities\InitialCondition\ShearWave\InitialConditionShearWave.h"

#include <sstream>

std::shared_ptr<TestParameter> ShearWaveTestParameter::getNewInstance(real u0, real v0, real viscosity, unsigned int lx, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, bool writeFiles, unsigned int startStepFileWriter, std::string filePath)
{
	return std::shared_ptr<TestParameter>(new ShearWaveTestParameter(u0, v0,
		viscosity, lx,
		numberOfTimeSteps, basisTimeStepLength,
		startStepCalculation, ySliceForCalculation,
		gridPath,
		writeFiles, startStepFileWriter, filePath));
}

ShearWaveTestParameter::ShearWaveTestParameter(real u0, real v0, real viscosity, unsigned int lx, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength, unsigned int startStepCalculation, unsigned int ySliceForCalculation, std::string gridPath, bool writeFiles, unsigned int startStepFileWriter, std::string filePath)
	:TestParameterImp(viscosity, lx, numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, gridPath, writeFiles, startStepFileWriter), u0(u0), v0(v0)
{
	std::ostringstream oss;
	oss << filePath + "\\ShearWave\\grid" << lx;
	this->filePath = oss.str();

	initialCondition = std::shared_ptr<InitialConditionShearWave>(new InitialConditionShearWave((double)lx, (double)lz, (double)l0, u0, v0, rho0));
}

std::shared_ptr<InitialCondition> ShearWaveTestParameter::getInitialCondition()
{
	return initialCondition;
}

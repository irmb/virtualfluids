#ifndef SHEARWAVETESTPARAMETER_H
#define SHEARWAVETESTPARAMETER_H

#include "../TestParameterImp.h"

class ShearWaveTestParameter : public TestParameterImp
{
public:
	static std::shared_ptr<TestParameter> getNewInstance(real u0, real v0,
		real viscosity, unsigned int lx,
		unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
		unsigned int startStepCalculation, unsigned int ySliceForCalculation,
		std::string gridPath,
		bool writeFiles, unsigned int startStepFileWriter, std::string filePath);

protected:
	ShearWaveTestParameter() {};
	ShearWaveTestParameter(real u0, real v0,
		real viscosity, unsigned int lx,
		unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
		unsigned int startStepCalculation, unsigned int ySliceForCalculation,
		std::string gridPath,
		bool writeFiles, unsigned int startStepFileWriter, std::string filePath);

	std::shared_ptr<InitialCondition> getInitialCondition();

private:
	std::shared_ptr<InitialCondition> initialCondition;

	real u0, v0;
};

#endif // !SHEARWAVETESTPARAMETER_H

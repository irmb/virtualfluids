#ifndef TGVTESTPARAMETER_H
#define TGVTESTPARAMETER_H

#include "../TestParameterImp.h"

#include <string>
#include <memory>

class TaylorGreenTestParameter : public TestParameterImp
{
public:
	static std::shared_ptr<TestParameter> getNewInstance(real u0, real amplitude,
		real viscosity, unsigned int lx,
		unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
		unsigned int startStepCalculation, unsigned int ySliceForCalculation,
		std::string gridPath,
		bool writeFiles, unsigned int startStepFileWriter, std::string filePath);
	
protected:
	TaylorGreenTestParameter(real u0, real amplitude,
		real viscosity, unsigned int lx,
		unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
		unsigned int startStepCalculation, unsigned int ySliceForCalculation,
		std::string gridPath,
		bool writeFiles, unsigned int startStepFileWriter, std::string filePath);

	std::shared_ptr<InitialCondition> getInitialCondition();

private:
	std::shared_ptr<InitialCondition> initialCondition;

	real u0, amplitude;

};
#endif // !TGVTESTPARAMETER_H

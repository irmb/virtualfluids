#ifndef SHEAR_WAVE_TEST_PARAMETER_H
#define SHEAR_WAVE_TEST_PARAMETER_H

#include "Utilities/TestParameter/TestParameterImp.h"
class  PhiAndNuTest;

class ShearWaveTestParameter : public TestParameterImp
{
public:
	static std::shared_ptr<TestParameter> getNewInstance(real u0, real v0, real viscosity, real rho0, 
														real lx, real lz, real l0,
														unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
														unsigned int startStepCalculation, unsigned int ySliceForCalculation,
														std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
														bool writeFiles, unsigned int startStepFileWriter, std::string filePath,
														std::shared_ptr<PhiAndNuTest> testResults, std::vector<int> devices);
	double getMaxVelocity();

protected:
	ShearWaveTestParameter() {};
	ShearWaveTestParameter(real u0, real v0, real viscosity, real rho0,
							real lx, real lz, real l0,
							unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
							unsigned int startStepCalculation, unsigned int ySliceForCalculation,
							std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
							bool writeFiles, unsigned int startStepFileWriter, std::string filePath,
							std::shared_ptr<PhiAndNuTest> testResults, std::vector<int> devices);

private:
	real u0, v0, rho0;
};

#endif

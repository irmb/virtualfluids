#ifndef TGV_TEST_PARAMETER_H
#define TGV_TEST_PARAMETER_H

#include "Utilities/TestParameter/TestParameterImp.h"

#include <string>
#include <memory>

class PhiAndNuTest;

class TaylorGreenTestParameter : public TestParameterImp
{
public:
	static std::shared_ptr<TestParameter> getNewInstance(real u0, real amplitude, real viscosity, real rho0,
														unsigned int lx, unsigned int lz, unsigned int l0,
														unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
														unsigned int startStepCalculation, unsigned int ySliceForCalculation,
														std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
														bool writeFiles, unsigned int startStepFileWriter, std::string filePath,
														std::shared_ptr<PhiAndNuTest> testResults, std::vector<int> devices);
	double getMaxVelocity();
	
protected:
	TaylorGreenTestParameter(real u0, real amplitude,
							real viscosity, real rho0, unsigned int lx, unsigned int lz, unsigned int l0,
							unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
							unsigned int startStepCalculation, unsigned int ySliceForCalculation,
							std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
							bool writeFiles, unsigned int startStepFileWriter, std::string filePath, 
							std::shared_ptr<PhiAndNuTest> testResults, std::vector<int> devices);

private:
	real u0, amplitude, rho0;

};
#endif 

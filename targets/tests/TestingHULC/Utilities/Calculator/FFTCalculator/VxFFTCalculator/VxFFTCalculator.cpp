#include "VxFFTCalculator.h"

#include "Utilities\Results\Results.h"

std::shared_ptr<VxFFTCalculator> VxFFTCalculator::getNewInstance(double viscosity, std::shared_ptr<PhiAndNuTestResults> testResults)
{
	return std::shared_ptr<VxFFTCalculator>(new VxFFTCalculator(viscosity, testResults));
}

void VxFFTCalculator::setVectorToCalc()
{
	data = simResults->getVx();
}

VxFFTCalculator::VxFFTCalculator(double viscosity, std::shared_ptr<PhiAndNuTestResults> testResults) : FFTCalculator(viscosity, testResults)
{
	
}
#ifndef VZFFTCALCULATOR_H
#define VZFFTCALCULATOR_H

#include "../FFTCalculator.h"

class VzFFTCalculator : public FFTCalculator
{
public:
	static std::shared_ptr<VzFFTCalculator> getNewInstance(double viscosity, std::shared_ptr<PhiAndNuTestResults> testResults);

protected:
	void setVectorToCalc();
private:
	VzFFTCalculator(double viscosity, std::shared_ptr<PhiAndNuTestResults> testResults);
};
#endif 
#ifndef VXFFTCALCULATOR_H
#define VXFFTCALCULATOR_H

#include "../FFTCalculator.h"

class VxFFTCalculator : public FFTCalculator
{
public:
	static std::shared_ptr<VxFFTCalculator> getNewInstance(double viscosity, std::shared_ptr<PhiAndNuTest> testResults);

protected:
	void setVectorToCalc();
private:
	VxFFTCalculator(double viscosity, std::shared_ptr<PhiAndNuTest> testResults);
};
#endif 
#ifndef L2NORM_CALCULATOR_IMP_H
#define L2NORM_CALCULATOR_IMP_H

#include "L2NormCalculator.h"

#include <memory>

class FFTCalculator;

class L2NormCalculatorImp : public L2NormCalculator
{
public:
	virtual double calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz, double timeStepLength) = 0;
	std::string getErrorMessage();

protected:
	L2NormCalculatorImp(std::string errorMessage);

	bool equalDouble(double num1, double num2);
	double calcCounter(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz, double timeStepLength);

	std::string errorMessage;

private:
	L2NormCalculatorImp();
};
#endif
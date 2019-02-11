#include "L2NormCalculator.h"

#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"
#include "Utilities/Calculator/AlmostEquals.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

std::shared_ptr<L2NormCalculator> L2NormCalculator::getInstance()
{
	static std::shared_ptr<L2NormCalculator> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<L2NormCalculator>(new L2NormCalculator());
	return uniqueInstance;
}

double L2NormCalculator::calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz, double timeStepLength)
{
	fftCalc = FFTCalculator::getNewInstance(lx, lz, timeStepLength);
	double amplitude = fftCalc->calcAmplitudeForTimeStep(basicData, false);
	if (equalDouble(amplitude, 0.0))
		return -1.0;
	double zaehler = 0.0;
	//double nenner = 0.0;
	for (int i = 0; i < basicData.size(); i++) {
		double area = (1 / pow(2.0, level.at(i))) * (1 / pow(2.0, level.at(i)));
		zaehler += ((divergentData.at(i) - basicData.at(i))*(divergentData.at(i) - basicData.at(i))) * area;
		//nenner += (basicData.at(i)*basicData.at(i)) * flaeche;
	}
	return sqrt(zaehler / (amplitude * amplitude));
}

bool L2NormCalculator::equalDouble(double num1, double num2)
{
	const FloatingPoint<double> lhs(num1), rhs(num2);

	if (lhs.AlmostEquals(rhs))
		return true;
	return false;
}

L2NormCalculator::L2NormCalculator()
{
	
}
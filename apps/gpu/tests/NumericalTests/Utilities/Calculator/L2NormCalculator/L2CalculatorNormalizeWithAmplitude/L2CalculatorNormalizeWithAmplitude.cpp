#include "L2CalculatorNormalizeWithAmplitude.h"

#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"

#include <cmath>

std::shared_ptr<L2NormCalculator> L2CalculatorNormalizeWithAmplitude::getInstance()
{
	static std::shared_ptr<L2NormCalculator> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<L2NormCalculator>(new L2CalculatorNormalizeWithAmplitude());
	return uniqueInstance;
}

double L2CalculatorNormalizeWithAmplitude::calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz, double l0)
{
	std::shared_ptr<FFTCalculator> fftCalc = FFTCalculator::getInstance();
	double amplitude = fftCalc->calcAmplitudeForTimeStep(basicData, false, lx, lz) * lx / l0;
	if (equalDouble(amplitude, 0.0))
		return -1.0;
	double counter = calcCounter(basicData, divergentData, level, lx, lz);
	return std::sqrt(counter / (amplitude * amplitude));
}

L2CalculatorNormalizeWithAmplitude::L2CalculatorNormalizeWithAmplitude() : L2NormCalculatorImp("Test could not run. Amplitude is zero. Normalization of the data is not possible.")
{
}

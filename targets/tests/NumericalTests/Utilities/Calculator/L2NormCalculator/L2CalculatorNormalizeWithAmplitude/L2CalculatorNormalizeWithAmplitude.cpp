#include "L2CalculatorNormalizeWithAmplitude.h"

#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"

std::shared_ptr<L2NormCalculator> L2CalculatorNormalizeWithAmplitude::getInstance()
{
	static std::shared_ptr<L2NormCalculator> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<L2NormCalculator>(new L2CalculatorNormalizeWithAmplitude());
	return uniqueInstance;
}

double L2CalculatorNormalizeWithAmplitude::calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz)
{
	std::shared_ptr<FFTCalculator> fftCalc = FFTCalculator::getInstance();
	double amplitude = fftCalc->calcAmplitudeForTimeStep(basicData, false, lx, lz);
	if (equalDouble(amplitude, 0.0))
		return -1.0;
	double counter = calcCounter(basicData, divergentData, level, lx, lz);
	return sqrt(counter / (amplitude * amplitude));
}

L2CalculatorNormalizeWithAmplitude::L2CalculatorNormalizeWithAmplitude() : L2NormCalculatorImp("Test could not run. Amplitude is zero. Normalization of the data is not possible.")
{
}

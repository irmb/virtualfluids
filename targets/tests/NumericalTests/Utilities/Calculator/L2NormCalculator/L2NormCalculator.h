#ifndef L2NORM_CALCULATOR_H
#define L2NORM_CALCULATOR_H

#include <vector>
#include <memory>

class FFTCalculator;

class L2NormCalculator
{
public:
	static std::shared_ptr< L2NormCalculator> getInstance();

	double calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz, double timeStepLength);

private:
	L2NormCalculator();

	std::shared_ptr< FFTCalculator> fftCalc;
};
#endif
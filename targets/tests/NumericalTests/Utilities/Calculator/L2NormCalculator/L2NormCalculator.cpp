#include "L2NormCalculator.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<L2NormCalculator> L2NormCalculator::getNewInstance()
{
	return std::shared_ptr<L2NormCalculator>(new L2NormCalculator());
}

double L2NormCalculator::calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level)
{
	double zaehler = 0.0;
	double nenner = 0.0;
	for (int i = 0; i < basicData.size(); i++) {
		double flaeche = (1 / pow(2.0, level.at(i))) * (1 / pow(2.0, level.at(i)));
		zaehler += ((divergentData.at(i) - basicData.at(i))*(divergentData.at(i) - basicData.at(i))) * flaeche;
		nenner += (basicData.at(i)*basicData.at(i)) * flaeche;
	}
	if (nenner == 0.0)
		return sqrt(zaehler);
	return sqrt(zaehler / nenner);
}

L2NormCalculator::L2NormCalculator()
{

}
#include "L2CalculatorNormalizeWithBasicData.h"

std::shared_ptr<L2NormCalculator> L2CalculatorNormalizeWithBasicData::getInstance()
{
	static std::shared_ptr<L2NormCalculator> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<L2NormCalculator>(new L2CalculatorNormalizeWithBasicData());
	return uniqueInstance;
}

double L2CalculatorNormalizeWithBasicData::calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz, double l0)
{
	double counter = calcCounter(basicData, divergentData, level, lx, lz);
	double denominator = 0.0;
	for (int i = 0; i < basicData.size(); i++) {
		double area = (1 / pow(2.0, level.at(i))) * (1 / pow(2.0, level.at(i)));
		denominator += (basicData.at(i)*basicData.at(i)) * area;
	}
	if (equalDouble(denominator, 0.0))
		return -1.0;

	return sqrt(counter / denominator);
}

L2CalculatorNormalizeWithBasicData::L2CalculatorNormalizeWithBasicData() : L2NormCalculatorImp("Test could not run. BasicData is zero. Normalization of the data is not possible.")
{

}

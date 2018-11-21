#include "L2NormCalculator.h"

std::shared_ptr<L2NormCalculator> L2NormCalculator::getNewInstance()
{
	return std::shared_ptr<L2NormCalculator>(new L2NormCalculator());
}

std::vector< double> L2NormCalculator::calc(std::vector<std::vector<double>> basicData, std::vector<std::vector<double>> divergentData)
{
	std::vector< double> result;

	return result;
}

L2NormCalculator::L2NormCalculator()
{

}
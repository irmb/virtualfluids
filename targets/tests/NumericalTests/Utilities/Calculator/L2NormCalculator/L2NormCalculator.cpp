#include "L2NormCalculator.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<L2NormCalculator> L2NormCalculator::getNewInstance()
{
	return std::shared_ptr<L2NormCalculator>(new L2NormCalculator());
}

std::vector< double> L2NormCalculator::calc(std::vector<std::vector<double>> basicData, std::vector<std::vector<double>> divergentData)
{
	std::vector< double> results;

	for (int i = 0; i < basicData.size(); i++) {
		double sum = 0;
		for (int j = 0; j < basicData.at(i).size(); j++) {
			sum += (divergentData.at(i).at(j) - basicData.at(i).at(j))*(divergentData.at(i).at(j) - basicData.at(i).at(j)) / (basicData.at(i).at(j)*basicData.at(i).at(j));
		}
		results.push_back(sqrt(sum));
	}
	
	return results;
}

L2NormCalculator::L2NormCalculator()
{

}
#include "L2NormCalculator.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<L2NormCalculator> L2NormCalculator::getNewInstance()
{
	return std::shared_ptr<L2NormCalculator>(new L2NormCalculator());
}

std::vector< double> L2NormCalculator::calc(std::vector<std::vector<double>> basicData, std::vector<std::vector<double>> divergentData, std::vector<std::vector<unsigned int>> level)
{
	std::vector< double> results;

	for (int i = 0; i < basicData.size(); i++) {
		double zaehler = 0;
		double nenner = 0;
		for (int j = 0; j < basicData.at(i).size(); j++) {
			double flaeche = (1 / pow(2.0, level.at(i).at(j))) * (1 / pow(2.0, level.at(i).at(j)));
			zaehler += ((divergentData.at(i).at(j) - basicData.at(i).at(j))*(divergentData.at(i).at(j) - basicData.at(i).at(j))) * flaeche;
			nenner += (basicData.at(i).at(j)*basicData.at(i).at(j)) * flaeche;
		}
		results.push_back(sqrt(zaehler/nenner));
	}
	
	return results;
}

L2NormCalculator::L2NormCalculator()
{

}
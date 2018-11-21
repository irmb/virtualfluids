#ifndef L2NORM_CALCULATOR_H
#define L2NORM_CALCULATOR_H

#include <vector>
#include <memory>

class L2NormCalculator
{
public:
	static std::shared_ptr< L2NormCalculator> getNewInstance();

	std::vector< double> calc(std::vector<std::vector<double>> basicData, std::vector<std::vector<double>> divergentData);

private:
	L2NormCalculator();
};
#endif
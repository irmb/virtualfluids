#ifndef L2_NORM_CALCULATOR_NORMALIZE_WITH_BASIC_DATA_H
#define L2_NORM_CALCULATOR_NORMALIZE_WITH_BASIC_DATA_H

#include "../L2NormCalculatorImp.h"

class L2CalculatorNormalizeWithBasicData : public L2NormCalculatorImp
{
public:
	static std::shared_ptr<L2NormCalculator> getInstance();

	double calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz, double l0);

private:
	L2CalculatorNormalizeWithBasicData();
};
#endif
#ifndef L2NORM_CALCULATOR_FACTORY_H
#define L2NORM_CALCULATOR_FACTORY_H

#include <memory>
#include <string>

class L2NormCalculator;

class L2NormCalculatorFactory
{
public:
	virtual std::shared_ptr<L2NormCalculator> makeL2NormCalculator(std::string type) = 0;

};
#endif 
#ifndef L2NORM_CALCULATOR_FACTORY_H
#define L2NORM_CALCULATOR_FACTORY_H

#include <memory>
#include <string>

class L2NormCalculator;

class L2NormCalculatorFactory
{
public:
	static std::shared_ptr<L2NormCalculatorFactory> getInstance();

	std::shared_ptr<L2NormCalculator> makeL2NormCalculator(std::string type);

private:
	L2NormCalculatorFactory();
};
#endif 
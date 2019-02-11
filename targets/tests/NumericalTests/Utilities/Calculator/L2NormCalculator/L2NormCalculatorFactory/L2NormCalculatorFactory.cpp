#include "L2NormCalculatorFactory.h"

#include "Utilities/Calculator/L2NormCalculator/L2CalculatorNormalizeWithAmplitude/L2CalculatorNormalizeWithAmplitude.h"
#include "Utilities/Calculator/L2NormCalculator/L2CalculatorNormalizeWithBasicData/L2CalculatorNormalizeWithBasicData.h"

std::shared_ptr<L2NormCalculatorFactory> L2NormCalculatorFactory::getInstance()
{
	static std::shared_ptr<L2NormCalculatorFactory> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<L2NormCalculatorFactory>(new L2NormCalculatorFactory());
	return uniqueInstance;
}

std::shared_ptr<L2NormCalculator> L2NormCalculatorFactory::makeL2NormCalculator(std::string type)
{
	if(type == "amplitude")
		return L2CalculatorNormalizeWithAmplitude::getInstance();
	if (type == "basicData")
		return L2CalculatorNormalizeWithBasicData::getInstance();
}

L2NormCalculatorFactory::L2NormCalculatorFactory()
{
}

#include "L2NormCalculatorFactoryImp.h"

#include "Utilities/Calculator/L2NormCalculator/L2CalculatorNormalizeWithAmplitude/L2CalculatorNormalizeWithAmplitude.h"
#include "Utilities/Calculator/L2NormCalculator/L2CalculatorNormalizeWithBasicData/L2CalculatorNormalizeWithBasicData.h"

std::shared_ptr<L2NormCalculatorFactory> L2NormCalculatorFactoryImp::getInstance()
{
    static std::shared_ptr<L2NormCalculatorFactory> uniqueInstance;
    if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<L2NormCalculatorFactory>(new L2NormCalculatorFactoryImp());
    return uniqueInstance;
}

std::shared_ptr<L2NormCalculator> L2NormCalculatorFactoryImp::makeL2NormCalculator(std::string type)
{
    if(type == "Amplitude")
        return L2CalculatorNormalizeWithAmplitude::getInstance();
    if (type == "BasicData")
        return L2CalculatorNormalizeWithBasicData::getInstance();
}

L2NormCalculatorFactoryImp::L2NormCalculatorFactoryImp()
{
}

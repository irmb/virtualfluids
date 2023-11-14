#ifndef L2NORM_CALCULATOR_FACTORY_IMP_H
#define L2NORM_CALCULATOR_FACTORY_IMP_H

#include "L2NormCalculatorFactory.h"

class L2NormCalculatorFactoryImp : public L2NormCalculatorFactory
{
public:
    static std::shared_ptr<L2NormCalculatorFactory> getInstance();

    std::shared_ptr<L2NormCalculator> makeL2NormCalculator(std::string type);

private:
    L2NormCalculatorFactoryImp();
};
#endif 
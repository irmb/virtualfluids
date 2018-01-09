#ifndef ConcreteCalculatorFactory_H
#define ConcreteCalculatorFactory_H


#include <memory>

#include "Grid3D.h"

#include "Calculator.h"
#include "MPICalculator.h"
#include "PrePostBcCalculator.h"

#include "CalculatorFactory.h"


class ConcreteCalculatorFactory : public CalculatorFactory
{
public:
    ConcreteCalculatorFactory(UbSchedulerPtr visScheduler) : CalculatorFactory(visScheduler){}
    virtual ~ConcreteCalculatorFactory() {}

    std::shared_ptr<Calculator> makeCalculator(Grid3DPtr grid, SynchronizerPtr sync, bool isMainThread, CalculatorType type) override
    {
        std::shared_ptr<Calculator> calculator;
        switch(type)
        {
        case CalculatorType::HYBRID:
            calculator = std::make_shared<Calculator>(grid, sync, isMainThread);
            break;
        case CalculatorType::MPI:
            calculator = std::make_shared<MPICalculator>(grid);
            break;
        case CalculatorType::PREPOSTBC:
            calculator = std::make_shared<PrePostBcCalculator>(grid, sync, isMainThread);
            break;
        #if defined CalculatorType::VF_FETOL
        case FETOL:
            calculator = std::make_shared<FetolCalculator>(grid);
            break;
        #endif
        default: 
            throw std::runtime_error("CalculatorType not valid in ConcreteCalculatorFactory");
        }
 
        calculator->setVisScheduler(visScheduler);
        return calculator;
    }


};



#endif 


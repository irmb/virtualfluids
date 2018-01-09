/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef DEM_CALCULATOR_FACTORY_H
#define DEM_CALCULATOR_FACTORY_H

#include <memory>

#include "CalculatorFactory.h"

#include "Calculator.h"
#include "Grid3D.h"

#include <dem_coupling/DemCalculator.h>

class DemCalculatorFactory : public CalculatorFactory
{
public:
    DemCalculatorFactory(UbSchedulerPtr visScheduler) : CalculatorFactory(visScheduler){}
    virtual ~DemCalculatorFactory() {}

    std::shared_ptr<Calculator> makeCalculator(Grid3DPtr grid, SynchronizerPtr sync, bool isMainThread, CalculatorType type) override
    {
        std::shared_ptr<Calculator> calculator;
        switch(type)
        {
        case CalculatorType::DEM:
            calculator = std::make_shared<DemCalculator>(grid, sync, isMainThread);
            break;
        default: 
            throw std::runtime_error("CalculatorType not valid in DemCalculatorFactory");
        }
 
        calculator->setVisScheduler(visScheduler);
        return calculator;
    }


};



#endif 


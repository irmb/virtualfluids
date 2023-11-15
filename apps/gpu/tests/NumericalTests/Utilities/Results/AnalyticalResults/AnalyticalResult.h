#ifndef ANALYTICAL_RESULTS_H
#define ANALYTICAL_RESULTS_H

#include "../ResultsImp.h"

#include <memory>

class SimulationResults;

class AnalyticalResults : public ResultsImp
{
public:
    virtual ~AnalyticalResults() = default;
    virtual void calc(std::shared_ptr<SimulationResults> simResults) = 0;
    virtual bool isCalculated() = 0;
};
#endif 
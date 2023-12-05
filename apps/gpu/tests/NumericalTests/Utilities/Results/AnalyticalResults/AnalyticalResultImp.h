#ifndef ANALYTICAL_RESULTS_IMP_H
#define ANALYTICAL_RESULTS_IMP_H

#include "AnalyticalResult.h"

#include "gpu/core/Calculation/Calculation.h"

class AnalyticalResultsImp : public AnalyticalResults
{
public:
    virtual void calc(std::shared_ptr<SimulationResults> simResults) = 0;
    bool isCalculated();
    int getL0();

protected:
    AnalyticalResultsImp(int l0);
    void init(std::shared_ptr<SimulationResults> simResults);

    bool calculated;
    int l0;

private:
    AnalyticalResultsImp();
};
#endif 
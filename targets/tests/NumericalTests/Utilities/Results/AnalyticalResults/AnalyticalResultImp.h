#ifndef ANALYTICAL_RESULTS_IMP_H
#define ANALYTICAL_RESULTS_IMP_H

#include "AnalyticalResult.h"

class AnalyticalResultsImp : public AnalyticalResults
{
public:
	virtual void calc(std::shared_ptr< SimulationResults> simResults) = 0;
	bool isCalculated();

protected:
	AnalyticalResultsImp();
	void init(std::shared_ptr< SimulationResults> simResults);

	bool calculated;
private:

};
#endif 
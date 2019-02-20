#ifndef SHEARWAVE_ANALYTICAL_RESULTS_H
#define SHEARWAVE_ANALYTICAL_RESULTS_H

#include "Utilities/Results/AnalyticalResults/AnalyticalResultImp.h"

struct ShearWaveParameterStruct;

class ShearWaveAnalyticalResults : public AnalyticalResultsImp
{
public:
	static std::shared_ptr<AnalyticalResults> getNewInstance(double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct);
	void calc(std::shared_ptr<SimulationResults> simResults);

private:
	ShearWaveAnalyticalResults() {};
	ShearWaveAnalyticalResults(double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct);

	double viscosity, rho0;
	double l0;
	double u0, v0;
};
#endif 
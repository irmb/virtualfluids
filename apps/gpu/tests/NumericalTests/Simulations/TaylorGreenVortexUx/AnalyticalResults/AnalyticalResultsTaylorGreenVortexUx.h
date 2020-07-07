#ifndef ANALYTICAL_RESULTS_TAYLORGREENVORTEX_U0_H
#define ANALYTICAL_RESULTS_TAYLORGREENVORTEX_U0_H

#include "Utilities/Results/AnalyticalResults/AnalyticalResultImp.h"

struct TaylorGreenVortexUxParameterStruct;

class AnalyticalResultsTaylorGreenUx : public AnalyticalResultsImp
{
public:
	static std::shared_ptr<AnalyticalResults> getNewInstance(double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct);
	void calc(std::shared_ptr<SimulationResults> simResults);

private:
	AnalyticalResultsTaylorGreenUx();
	AnalyticalResultsTaylorGreenUx(double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct);

	double viscosity, rho0;
	double l0;
	double ux, amplitude; 
};
#endif 
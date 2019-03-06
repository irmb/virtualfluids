#ifndef ANALYTICAL_RESULTS_TAYLORGREENVORTEX_Uz_H
#define ANALYTICAL_RESULTS_TAYLORGREENVORTEX_Uz_H

#include "Utilities/Results/AnalyticalResults/AnalyticalResultImp.h"

struct TaylorGreenVortexUzParameterStruct;

class AnalyticalResultsTaylorGreenUz : public AnalyticalResultsImp
{
public:
	static std::shared_ptr<AnalyticalResults> getNewInstance(double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct);
	void calc(std::shared_ptr<SimulationResults> simResults);


private:
	AnalyticalResultsTaylorGreenUz() {};
	AnalyticalResultsTaylorGreenUz(double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct);

	double viscosity, rho0;
	double l0;
	double uz, amplitude; 
};
#endif 
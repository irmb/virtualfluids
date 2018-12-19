#ifndef ANALYTICAL_RESULTS_TAYLORGREENVORTEX_Uz_H
#define ANALYTICAL_RESULTS_TAYLORGREENVORTEX_Uz_H

#include "Utilities\Results\AnalyticalResults\AnalyticalResultImp.h"

class AnalyticalResultsTaylorGreenUz : public AnalyticalResultsImp
{
public:
	static std::shared_ptr< AnalyticalResults> getNewInstance(double viscosity, double uz, double amplitude, double l0, double rho0);
	void calc(std::shared_ptr< SimulationResults> simResults);


private:
	AnalyticalResultsTaylorGreenUz() {};
	AnalyticalResultsTaylorGreenUz(double viscosity, double uz, double amplitude, double l0, double rho0);

	double viscosity, rho0;
	double l0;
	double uz, amplitude; 
};
#endif 
#ifndef ANALYTICAL_RESULTS_TAYLORGREENVORTEX_U0_H
#define ANALYTICAL_RESULTS_TAYLORGREENVORTEX_U0_H

#include "Utilities\Results\AnalyticalResults\AnalyticalResultImp.h"

class AnalyticalResultsTaylorGreenUx : public AnalyticalResultsImp
{
public:
	static std::shared_ptr< AnalyticalResults> getNewInstance(double viscosity, double ux, double amplitude, double l0, double rho0);
	void calc(std::shared_ptr< SimulationResults> simResults);

private:
	AnalyticalResultsTaylorGreenUx() {};
	AnalyticalResultsTaylorGreenUx(double viscosity, double ux, double amplitude, double l0, double rho0);

	double viscosity, rho0;
	double l0;
	double ux, amplitude; 
};
#endif 
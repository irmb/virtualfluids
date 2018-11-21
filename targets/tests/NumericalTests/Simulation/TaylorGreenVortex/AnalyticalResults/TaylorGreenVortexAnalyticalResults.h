#ifndef TAYLORGREENVORTEX_ANALYTICAL_RESULTS_H
#define TAYLORGREENVORTEX_ANALYTICAL_RESULTS_H

#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"

class TaylorGreenAnalyticalResults : public AnalyticalResults
{
public:
	static std::shared_ptr< AnalyticalResults> getNewInstance(double viscosity, double u0, double amplitude, double l0, double rho0);
	void calc(std::shared_ptr< SimulationResults> simResults);


private:
	TaylorGreenAnalyticalResults() {};
	TaylorGreenAnalyticalResults(double viscosity, double u0, double amplitude, double l0, double rho0);

	double viscosity, rho0;
	double l0;
	double u0, amplitude; 
};
#endif 
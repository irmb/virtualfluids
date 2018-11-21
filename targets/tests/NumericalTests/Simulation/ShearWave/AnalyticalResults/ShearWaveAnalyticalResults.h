#ifndef SHEARWAVE_ANALYTICAL_RESULTS_H
#define SHEARWAVE_ANALYTICAL_RESULTS_H

#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"

class ShearWaveAnalyticalResults : public AnalyticalResults
{
public:
	static std::shared_ptr< AnalyticalResults> getNewInstance(double viscosity, double u0, double v0, double l0, double rho0);
	void calc(std::shared_ptr< SimulationResults> simResults);


private:
	ShearWaveAnalyticalResults() {};
	ShearWaveAnalyticalResults(double viscosity, double u0, double v0, double l0, double rho0);

	double viscosity, rho0;
	double l0;
	double u0, v0;
};
#endif 
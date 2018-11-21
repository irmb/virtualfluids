#include "ShearWaveAnalyticalResults.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<AnalyticalResults> ShearWaveAnalyticalResults::getNewInstance(double viscosity, double u0, double v0, double l0, double rho0)
{
	return std::shared_ptr<AnalyticalResults>(new ShearWaveAnalyticalResults(viscosity, u0, v0, l0, rho0));
}

void ShearWaveAnalyticalResults::calc(std::shared_ptr<SimulationResults> simResults)
{
	AnalyticalResults::init(simResults);

	for (int i = 0; i < numberOfTimeSteps; i++) {
		for (int j = 0; j < numberOfNodes; j++) {
			vx.at(i).at(j) = (l0*u0) / xNodes;
			vy.at(i).at(j) = (double)0.0;
			vz.at(i).at(j) = (l0*v0*cos(((double)2.0 * M_PI*z.at(i).at(j)) / zNodes) * sin(((double)2.0 * M_PI*(x.at(i).at(j) + (l0*timeStep.at(i)*u0) / xNodes)) / xNodes)) / (exp(timeStep.at(i)*viscosity*(((double)4.0 * M_PI*M_PI) / xNodes*xNodes + ((double)4.0 * M_PI*M_PI) / zNodes*zNodes))*xNodes);
			press.at(i).at(j) = (double)0.0;
			rho.at(i).at(j) = (l0*l0 * rho0*v0*sin(((double)2.0 * M_PI*z.at(i).at(j)) / zNodes) * ((double)-4.0 * exp(((double)4.0 * M_PI*M_PI * timeStep.at(i)*viscosity*(xNodes*xNodes + zNodes*zNodes)) / (xNodes*xNodes * zNodes*zNodes))*u0*zNodes*cos(((double)2.0 * M_PI*(l0*timeStep.at(i)*u0 + x.at(i).at(j)*xNodes)) / (xNodes*xNodes)) + v0*xNodes*sin(((double)2.0 * M_PI*(l0*timeStep.at(i)*u0 + x.at(i).at(j)*xNodes)) / (xNodes*xNodes))*sin(((double)2.0 * M_PI*(l0*timeStep.at(i)*u0 + x.at(i).at(j)*xNodes)) / (xNodes*xNodes)) * sin((2 * M_PI*z.at(i).at(j)) / zNodes))) / ((double)2.0 * exp(((double)8.0 * M_PI*M_PI * timeStep.at(i)*viscosity*(xNodes*xNodes + zNodes*zNodes)) / (xNodes*xNodes * zNodes*zNodes))*xNodes*xNodes*xNodes);
		}
	}
}

ShearWaveAnalyticalResults::ShearWaveAnalyticalResults(double viscosity, double u0, double v0, double l0, double rho0) : viscosity(viscosity), u0(u0), v0(v0), l0(l0), rho0(rho0)
{

}
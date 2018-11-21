#include "TaylorGreenVortexAnalyticalResults.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<AnalyticalResults> TaylorGreenAnalyticalResults::getNewInstance(double viscosity, double u0, double amplitude, double l0, double rho0)
{
	return std::shared_ptr<AnalyticalResults>(new TaylorGreenAnalyticalResults(viscosity, u0, amplitude, l0, rho0));
}

void TaylorGreenAnalyticalResults::calc(std::shared_ptr< SimulationResults> simResults)
{
	AnalyticalResults::init(simResults);

	for (int i = 0; i < numberOfTimeSteps; i++) {
		for (int j = 0; j < numberOfNodes; j++) {
			vx.at(i).at(j) = (l0*u0) / xNodes + (amplitude * exp(timeStep.at(i)*viscosity*(((double)-4.0 * M_PI*M_PI) / (xNodes*xNodes) - ((double)4.0 * M_PI *M_PI) / (zNodes*zNodes)))*l0*cos(((double)2.0 * M_PI * z.at(i).at(j)) / zNodes) * sin(((double)2.0 * M_PI*x.at(i).at(j)) / xNodes)) / xNodes;
			vy.at(i).at(j) = (double)0.0;
			vz.at(i).at(j) = -((amplitude*exp(timeStep.at(i)*viscosity*(((double)-4.0 * M_PI*M_PI) / (xNodes*xNodes) - ((double)4.0 * M_PI*M_PI) / (zNodes*zNodes)))*l0*zNodes*cos(((double)2 * M_PI*x.at(i).at(j)) / xNodes) * sin((2 * M_PI*z.at(i).at(j)) / zNodes)) / (xNodes*xNodes));
			press.at(i).at(j) = (double)0.0;
			rho.at(i).at(j) = (amplitude*l0*l0 * rho0*(amplitude * zNodes*zNodes * cos(((double)2.0 * M_PI*z.at(i).at(j)) / zNodes)*cos(((double)2.0 * M_PI*z.at(i).at(j)) / zNodes) - (double)2.0 * exp(((double)4.0 * M_PI*M_PI * timeStep.at(i)*(xNodes*xNodes + zNodes*zNodes)*viscosity) / (xNodes*xNodes * zNodes*zNodes))*u0*(xNodes*xNodes - zNodes*zNodes)*cos(((double)2.0 * M_PI*z.at(i).at(j)) / zNodes) * sin(((double)2.0 * M_PI*x.at(i).at(j)) / xNodes) - amplitude*xNodes*xNodes * sin(((double)2.0 * M_PI*x.at(i).at(j)) / xNodes)*sin(((double)2.0 * M_PI*x.at(i).at(j)) / xNodes))) / ((double)2.0 * exp(((double)8.0 * M_PI*M_PI * timeStep.at(i)*(xNodes*xNodes + zNodes*zNodes)*viscosity) / (xNodes*xNodes * zNodes*zNodes))*xNodes*xNodes*xNodes*xNodes);
		}
	}
}

TaylorGreenAnalyticalResults::TaylorGreenAnalyticalResults(double viscosity, double u0, double amplitude, double l0, double rho0) : viscosity(viscosity), u0(u0), amplitude(amplitude), l0(l0), rho0(rho0)
{
	
}
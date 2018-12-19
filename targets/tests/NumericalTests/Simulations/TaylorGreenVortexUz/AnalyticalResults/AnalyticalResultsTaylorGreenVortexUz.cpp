#include "AnalyticalResultsTaylorGreenVortexUz.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<AnalyticalResults> AnalyticalResultsTaylorGreenUz::getNewInstance(double viscosity, double uz, double amplitude, double l0, double rho0)
{
	return std::shared_ptr<AnalyticalResults>(new AnalyticalResultsTaylorGreenUz(viscosity, uz, amplitude, l0, rho0));
}

void AnalyticalResultsTaylorGreenUz::calc(std::shared_ptr< SimulationResults> simResults)
{
	AnalyticalResultsImp::init(simResults);

	for (int i = 0; i < numberOfTimeSteps; i++) {
		for (int j = 0; j < numberOfNodes; j++) {
			vx.at(i).at(j) = (amplitude*exp( time.at(i)*viscosity*((-(double)4.0 * pow(M_PI, (double)2.0)) / pow(xNodes, (double)2.0) - ((double)4.0 * pow(M_PI, (double)2.0)) / pow(zNodes, (double)2.0)))*l0*cos(((double)2.0 * M_PI*((l0*time.at(i)*uz) / xNodes + z.at(i).at(j))) / zNodes)*sin(((double)2.0 * M_PI*x.at(i).at(j)) / xNodes)) / xNodes;
			vy.at(i).at(j) = (double)0.0;
			vz.at(i).at(j) = (l0*uz) / xNodes - (amplitude*exp( time.at(i)*viscosity*((-(double)4.0 * pow(M_PI, (double)2.0)) / pow(xNodes, (double)2.0) - ((double)4.0 * pow(M_PI, (double)2.0)) / pow(zNodes, (double)2.0)))*l0*zNodes*cos(((double)2.0 * M_PI*x.at(i).at(j)) / xNodes)*sin(((double)2.0 * M_PI*((l0*time.at(i)*uz) / xNodes + z.at(i).at(j))) / zNodes)) / pow(xNodes, (double)2.0);
			press.at(i).at(j) = (double)0.0;
			rho.at(i).at(j) = (amplitude*pow(l0, (double)2.0)*rho0*(amplitude*pow(xNodes, (double)2.0)*zNodes*pow(cos(((double)2.0 * M_PI*x.at(i).at(j)) / xNodes), (double)2.0) - (double)4.0 * exp( ((double)4.0 * pow(M_PI, (double)2.0)*time.at(i)*viscosity*(pow(xNodes, (double)2.0) + pow(zNodes, (double)2.0))) / (pow(xNodes, (double)2.0)*pow(zNodes, (double)2.0)))*uz*xNodes * (pow(xNodes, (double)2.0) - pow(zNodes, (double)2.0))*cos(((double)2.0 * M_PI*x.at(i).at(j)) / xNodes)*sin(((double)2.0 * M_PI*((l0*time.at(i)*uz) / xNodes + z.at(i).at(j))) / zNodes) - amplitude*pow(zNodes, 3)*pow(sin(((double)2.0 * M_PI*((l0*time.at(i)*uz) / xNodes + z.at(i).at(j))) / zNodes), (double)2.0))) / ((double)2.0 *exp( (8 * pow(M_PI, (double)2.0)*time.at(i)*viscosity*(pow(xNodes, (double)2.0) + pow(zNodes, (double)2.0))) / (pow(xNodes, (double)2.0)*pow(zNodes, (double)2.0)))*pow(xNodes, (double)4.0)*zNodes);
		}
	}
	calculated = true;
}

AnalyticalResultsTaylorGreenUz::AnalyticalResultsTaylorGreenUz(double viscosity, double uz, double amplitude, double l0, double rho0) : AnalyticalResultsImp(), viscosity(viscosity), uz(uz), amplitude(amplitude), l0(l0), rho0(rho0)
{
	
}
#include "AnalyticalResultsTaylorGreenVortexUz.h"

#include "Simulations/TaylorGreenVortexUz/TaylorGreenVortexUzParameterStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<AnalyticalResults> AnalyticalResultsTaylorGreenUz::getNewInstance(double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct)
{
	return std::shared_ptr<AnalyticalResults>(new AnalyticalResultsTaylorGreenUz(viscosity, simParaStruct));
}

void AnalyticalResultsTaylorGreenUz::calc(std::shared_ptr<SimulationResults> simResults)
{
	AnalyticalResultsImp::init(simResults);

	for (int i = 0; i < numberOfTimeSteps; i++) {
		for (int j = 0; j < numberOfNodes; j++) {
			vx.at(i).at(j) = (amplitude*exp( time.at(i)*viscosity*((-(real)4.0 * pow(M_PI, (real)2.0)) / pow(xNodes, (real)2.0) - ((real)4.0 * pow(M_PI, (real)2.0)) / pow(zNodes, (real)2.0)))*l0*cos(((real)2.0 * M_PI*((l0*time.at(i)*uz) / xNodes + z.at(i).at(j))) / zNodes)*sin(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes)) / xNodes;
			vy.at(i).at(j) = (real)0.0;
			vz.at(i).at(j) = (l0*uz) / zNodes - (amplitude*exp(time.at(i)*viscosity*((-(real)4.0 * pow(M_PI, (real)2.0)) / pow(xNodes, (real)2.0) - ((real)4.0 * pow(M_PI, (real)2.0)) / pow(zNodes, (real)2.0)))*l0*zNodes*cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes)*sin(((real)2.0 * M_PI*(z.at(i).at(j) + (l0*time.at(i)*uz) / zNodes)) / zNodes)) / pow(xNodes, (real)2.0);
			press.at(i).at(j) = (amplitude*pow(l0, (real)2.0)*rho0*(amplitude*pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)*pow(cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes), (real)2.0) - (real)4.0 * exp(((real)4.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*uz*pow(xNodes, (real)2.0)*(pow(xNodes, (real)2.0) - pow(zNodes, (real)2.0))*cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes)*sin(((real)2.0 * M_PI*(l0*time.at(i)*uz + z.at(i).at(j)*zNodes)) / pow(zNodes, (real)2.0)) - amplitude*pow(zNodes, (real)4.0)*pow(sin(((real)2.0 * M_PI*(l0*time.at(i)*uz + z.at(i).at(j)*zNodes)) / pow(zNodes, (real)2.0)), (real)2.0))) / ((real)2.0*exp(((real)8.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*pow(xNodes, (real)4.0)*pow(zNodes, (real)2.0));
			rho.at(i).at(j) = (amplitude*pow(l0, (real)2.0)*rho0*(amplitude*pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)*pow(cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes), (real)2.0) - (real)4.0 * exp(((real)4.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*uz*pow(xNodes, (real)2.0)*(pow(xNodes, (real)2.0) - pow(zNodes, (real)2.0))*cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes)*sin(((real)2.0 * M_PI*(l0*time.at(i)*uz + z.at(i).at(j)*zNodes)) / pow(zNodes, (real)2.0)) - amplitude*pow(zNodes, (real)4.0)*pow(sin(((real)2.0 * M_PI*(l0*time.at(i)*uz + z.at(i).at(j)*zNodes)) / pow(zNodes, (real)2.0)), (real)2.0))) / ((real)2.0*exp(((real)8.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*pow(xNodes, (real)4.0)*pow(zNodes, (real)2.0));
		}
	}
	calculated = true;
}

AnalyticalResultsTaylorGreenUz::AnalyticalResultsTaylorGreenUz(double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct) : AnalyticalResultsImp()
{
	this->viscosity = viscosity;
	this->uz = simParaStruct->uz;
	this->amplitude = simParaStruct->amplitude;
	this->l0 = simParaStruct->l0;
	this->rho0 = simParaStruct->rho0;
}
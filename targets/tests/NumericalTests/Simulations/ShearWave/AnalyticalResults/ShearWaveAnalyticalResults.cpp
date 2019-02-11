#include "ShearWaveAnalyticalResults.h"

#include "Simulations/ShearWave/ShearWaveParameterStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<AnalyticalResults> ShearWaveAnalyticalResults::getNewInstance(double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct)
{
	return std::shared_ptr<AnalyticalResults>(new ShearWaveAnalyticalResults(viscosity, simParaStruct));
}

void ShearWaveAnalyticalResults::calc(std::shared_ptr<SimulationResults> simResults)
{
	AnalyticalResultsImp::init(simResults);

	for (int i = 0; i < numberOfTimeSteps; i++) {
		for (int j = 0; j < numberOfNodes; j++) {
			vx.at(i).at(j) = (l0*u0) / xNodes;
			vy.at(i).at(j) = (double)0.0;
			vz.at(i).at(j) = (l0*v0*cos(((double)2.0 * M_PI*z.at(i).at(j)) / zNodes)*sin(((double)2.0 * M_PI*(x.at(i).at(j) + (l0*time.at(i)*u0) / xNodes)) / xNodes)) / (exp(time.at(i)*viscosity*(((double)4.0 * pow(M_PI, (double)2.0)) / pow(xNodes, (double)2.0) + ((double)4.0 * pow(M_PI, (double)2.0)) / pow(zNodes, (double)2.0)))*xNodes);
			press.at(i).at(j) = (double)0.0;
			rho.at(i).at(j) = (pow(l0, (double)2.0)*rho0*v0*sin(((double)2.0 * M_PI*z.at(i).at(j)) / zNodes)*(-(double)4.0 * exp(((double)4.0 * pow(M_PI, (double)2.0)*time.at(i)*viscosity*(pow(xNodes, (double)2.0) + pow(zNodes, (double)2.0))) / (pow(xNodes, (double)2.0)*pow(zNodes, (double)2.0)))*u0*zNodes*cos(((double)2.0 * M_PI*(l0*time.at(i)*u0 + x.at(i).at(j)*xNodes)) / pow(xNodes, (double)2.0)) + v0*xNodes*pow(sin(((double)2.0 * M_PI*(l0*time.at(i)*u0 + x.at(i).at(j)*xNodes)) / pow(xNodes, (double)2.0)), (double)2.0)*sin(((double)2.0 * M_PI*z.at(i).at(j)) / zNodes))) / ((double)2.0*exp(((double)8.0 * pow(M_PI, (double)2.0)*time.at(i)*viscosity*(pow(xNodes, (double)2.0) + pow(zNodes, (double)2.0))) / (pow(xNodes, (double)2.0)*pow(zNodes, (double)2.0)))*pow(xNodes, (double)3.0));
		}
	}
	calculated = true;
}

ShearWaveAnalyticalResults::ShearWaveAnalyticalResults(double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct)
	: AnalyticalResultsImp()
{
	this->viscosity = viscosity;
	this->u0 = simParaStruct->ux;
	this->v0 = simParaStruct->uz;
	this->l0 = simParaStruct->l0;
	this->rho0 = simParaStruct->rho0;
}
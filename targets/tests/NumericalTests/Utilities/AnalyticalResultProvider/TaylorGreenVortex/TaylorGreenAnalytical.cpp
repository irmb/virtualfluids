/*#include "TaylorGreenAnalytical.h"

#include "../../Results/Results.h"
#include "utilities/TestCondition/TaylorGreenVortex/TaylorGreenTestCondition.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>


TaylorGreenAnalytical::TaylorGreenAnalytical(std::vector< std::shared_ptr<Results> > simulationResults, std::shared_ptr<TaylorGreenTestCondition> testCondition)
{
	/*Amp = testCondition->getAmp();
	L = testCondition->getL();
	L0 = testCondition->getL0();
	rho0 = testCondition->getRho0();
	vis = testCondition->getViscosity();
	Lx = testCondition->getLx();
	Lz = testCondition->getLz();
	numberNodes = Lx * Lz;
	this->init(simulationResults);
	this->calculate();
}

std::vector<std::shared_ptr<Results>> TaylorGreenAnalytical::getAnalyticalResults()
{
	return analyticalResults;
}

void TaylorGreenAnalytical::init(std::vector< std::shared_ptr<Results> > simulationResults)
{
	std::cout << "Initialize analytical Solution...\n";
	analyticalResults.resize(simulationResults.size());
	for (int i = 0; i < analyticalResults.size(); i++)
	{
		std::cout << "Initialize analytical solution TaylorGreenVortex 2D	Timestep: " << i << std::endl;
		analyticalResults[i] = std::shared_ptr<Results>(new Results);
		analyticalResults[i]->setNumberofNodes(numberNodes);
		
		analyticalResults[i]->setTime(simulationResults[i]->getTime());
		analyticalResults[i]->setTimeStep(simulationResults[i]->getTimeStep());
		for (int j = 0; j < numberNodes; j++)
		{
			analyticalResults[i]->setX(j, simulationResults[i]->getX()[j]);
			analyticalResults[i]->setY(j, simulationResults[i]->getY()[j]);
			analyticalResults[i]->setZ(j, simulationResults[i]->getZ()[j]);
		}
	}

}

void TaylorGreenAnalytical::calculate()
{
	for (int i = 0; i < analyticalResults.size(); i++)
	{
		std::cout << "Calculate analytical solution TaylorGreenVortex 2D	Timestep: " << i << std::endl;
		t = analyticalResults[i]->getTime();
		for (int j = 0; j < numberNodes; j++)
		{
			x = analyticalResults[i]->getX()[j];
			z = analyticalResults[i]->getZ()[j];
			
			vx = (Amp * L0 * cos((real)2.0 * M_PI * z * Lx / (Lz * L)) * sin((real)2.0 * M_PI * x / L) / L) * exp((real) -4.0 * (Lx*Lx + Lz*Lz) * M_PI*M_PI * t * vis / (L*L * Lz*Lz));
			vz = (-Amp * L0 * Lz * cos((real)2.0 * M_PI * x / L) * sin((real)2.0 * M_PI * z * Lx / (L * Lz)) / (L*Lx)) * exp((real)-4.0 * (Lx*Lx + Lz*Lz) * M_PI*M_PI * t * vis / (L*L * Lz*Lz));
			press = (Amp*Amp * L0*L0 * Lz*Lz * rho0 / ((real)16.0 * Lx*Lx * M_PI*M_PI)) * (((real)4.0 * Lx*Lx / (L*L * Lz*Lz)) * M_PI*M_PI * cos(((real)4.0 / L) * M_PI * x) + ((real)4.0 / (L*L)) * M_PI*M_PI * cos(((real)4.0 * Lx / (L * Lz)) * M_PI * z)) * exp((real)-8.0 * (Lx*Lx * Lz*Lz) * M_PI*M_PI * t * vis /(L*L * Lz*Lz));

			analyticalResults[i]->setVx(j, vx);
			analyticalResults[i]->setVy(j, (real) 0.0);
			analyticalResults[i]->setVz(j, vz);
			analyticalResults[i]->setPress(j, press);

		}
	}

}*/

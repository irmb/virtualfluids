#include "Results.h"

#define _USE_MATH_DEFINES
#include <math.h>

Results::Results(unsigned int lx, unsigned int lz, unsigned int timeStepLength)
{
	this->xNodes = lx;
	this->zNodes = lz;
	this->numberOfNodes = lx*lz;
	this->timeStepLength = timeStepLength;
	this->numberOfTimeSteps = 0;
}

void Results::addTimeStep(unsigned int timeStep, unsigned int time, std::vector<double> x, std::vector<double> z, std::vector<double> vx, std::vector<double> vz, std::vector<double> press, std::vector<double> rho)
{
	this->timeStep.push_back(timeStep);
	this->time.push_back(time);
	this->x.push_back(x);
	this->z.push_back(z);
	this->vx.push_back(vx);
	this->vz.push_back(vz);
	this->press.push_back(press);
	this->rho.push_back(rho);
	numberOfTimeSteps++;
}

int Results::getNumberOfTimeSteps()
{
	return numberOfTimeSteps;
}

std::vector< std::vector<double> > Results::getVx()
{
	return vx;
}

std::vector<std::vector<double>> Results::getVz()
{
	return vz;
}

int Results::getXNodes()
{
	return xNodes;
}

int Results::getZNodes()
{
	return zNodes;
}

int Results::getTimeStepLength()
{
	return timeStepLength;
}

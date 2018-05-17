#include "SimulationResults.h"

#define _USE_MATH_DEFINES
#include <math.h>

SimulationResults::SimulationResults(unsigned int lx, unsigned int lz, unsigned int timeStepLength)
{
	this->xNodes = lx;
	this->zNodes = lz;
	this->numberOfNodes = lx*lz;
	this->timeStepLength = timeStepLength;
	this->numberOfTimeSteps = 0;
}

std::shared_ptr<SimulationResults> SimulationResults::getNewInstance(unsigned int lx, unsigned int lz, unsigned int timeStepLength)
{
	return std::shared_ptr<SimulationResults>(new SimulationResults(lx, lz, timeStepLength));
}

void SimulationResults::addTimeStep(unsigned int timeStep, unsigned int time, std::vector<double> x, std::vector<double> z, std::vector<double> vx, std::vector<double> vz, std::vector<double> press, std::vector<double> rho)
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

int SimulationResults::getNumberOfTimeSteps()
{
	return numberOfTimeSteps;
}

std::vector< std::vector<double> > SimulationResults::getVx()
{
	return vx;
}

std::vector<std::vector<double>> SimulationResults::getVz()
{
	return vz;
}

int SimulationResults::getXNodes()
{
	return xNodes;
}

int SimulationResults::getZNodes()
{
	return zNodes;
}

int SimulationResults::getTimeStepLength()
{
	return timeStepLength;
}

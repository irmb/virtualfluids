#include "SimulationResults.h"

#define _USE_MATH_DEFINES
#include <math.h>

SimulationResults::SimulationResults(unsigned int lx, unsigned int ly, unsigned int lz, unsigned int timeStepLength)
{
	this->xNodes = lx;
	this->yNodes = ly;
	this->zNodes = lz;
	this->numberOfNodes = lx*ly*lz;
	this->timeStepLength = timeStepLength;
	this->numberOfTimeSteps = 0;
}

std::shared_ptr<SimulationResults> SimulationResults::getNewInstance(unsigned int lx, unsigned int ly, unsigned int lz, unsigned int timeStepLength)
{
	return std::shared_ptr<SimulationResults>(new SimulationResults(lx, ly, lz, timeStepLength));
}

void SimulationResults::addTimeStep(unsigned int timeStep, unsigned int time, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> vx, std::vector<double> vy, std::vector<double> vz, std::vector<double> press, std::vector<double> rho)
{
	this->timeStep.push_back(timeStep);
	this->time.push_back(time);
	this->x.push_back(x);
	this->y.push_back(y);
	this->z.push_back(z);
	this->vx.push_back(vx);
	this->vy.push_back(vy);
	this->vz.push_back(vz);
	this->press.push_back(press);
	this->rho.push_back(rho);
	numberOfTimeSteps++;
}
#include "SimulationResults.h"

#include "Utilities/SimulationParameter/SimulationParameter.h"

#define _USE_MATH_DEFINES
#include <math.h>

SimulationResults::SimulationResults(std::shared_ptr<SimulationParameter> simPara)
{
	this->xNodes = simPara->getLx();
	this->yNodes = 1;
	this->zNodes = simPara->getLz();
	this->numberOfNodes = xNodes*yNodes*zNodes;
	this->timeStepLength = simPara->getTimeStepLength();
	this->numberOfTimeSteps = 0;
}

std::shared_ptr<SimulationResults> SimulationResults::getNewInstance(std::shared_ptr<SimulationParameter> simPara)
{
	return std::shared_ptr<SimulationResults>(new SimulationResults(simPara));
}

void SimulationResults::addTimeStep(unsigned int timeStep, unsigned int time, std::vector<unsigned int> level, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> vx, std::vector<double> vy, std::vector<double> vz, std::vector<double> press, std::vector<double> rho)
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
	this->level.push_back(level);
	numberOfTimeSteps++;
}
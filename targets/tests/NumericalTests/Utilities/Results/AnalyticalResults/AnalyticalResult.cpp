#include "AnalyticalResult.h"

#include "../SimulationResults/SimulationResults.h"

void AnalyticalResults::init(std::shared_ptr<SimulationResults> simResults)
{
	this->xNodes = simResults->getNumberOfXNodes();
	this->yNodes = simResults->getNumberOfYNodes();
	this->zNodes = simResults->getNumberOfZNodes();
	this->numberOfNodes = xNodes*yNodes*zNodes;
	this->timeStepLength = simResults->getTimeStepLength();
	this->numberOfTimeSteps = simResults->getNumberOfTimeSteps();
	this->timeStep = simResults->getTimeSteps();
	this->x = simResults->getXNodes();
	this->y = simResults->getYNodes();
	this->z = simResults->getZNodes();

	this->vx.resize(numberOfTimeSteps);
	this->vy.resize(numberOfTimeSteps);
	this->vz.resize(numberOfTimeSteps);
	this->press.resize(numberOfTimeSteps);
	this->rho.resize(numberOfTimeSteps);

	for (int i = 0; i < numberOfTimeSteps; i++) {
		this->vx.at(i).resize(numberOfNodes);
		this->vy.at(i).resize(numberOfNodes);
		this->vz.at(i).resize(numberOfNodes);
		this->press.at(i).resize(numberOfNodes);
		this->rho.at(i).resize(numberOfNodes);
	}
}
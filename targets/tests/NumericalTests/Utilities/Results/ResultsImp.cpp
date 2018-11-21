#include "ResultsImp.h"

int ResultsImp::getNumberOfTimeSteps()
{
	return numberOfTimeSteps;
}

std::vector<std::vector<double>> ResultsImp::getVx()
{
	return vx;
}

std::vector<std::vector<double>> ResultsImp::getVy()
{
	return vy;
}

std::vector<std::vector<double>> ResultsImp::getVz()
{
	return vz;
}

int ResultsImp::getNumberOfXNodes()
{
	return xNodes;
}

int ResultsImp::getNumberOfYNodes()
{
	return yNodes;
}

int ResultsImp::getNumberOfZNodes()
{
	return zNodes;
}

std::vector<std::vector<double>> ResultsImp::getXNodes()
{
	return x;
}

std::vector<std::vector<double>> ResultsImp::getYNodes()
{
	return y;
}

std::vector<std::vector<double>> ResultsImp::getZNodes()
{
	return z;
}

int ResultsImp::getTimeStepLength()
{
	return timeStepLength;
}

std::vector<unsigned int> ResultsImp::getTimeSteps()
{
	return timeStep;
}
#include "ResultsImp.h"

#include <iostream>

int ResultsImp::getNumberOfTimeSteps()
{
	return numberOfTimeSteps;
}

std::vector<std::vector<double> > ResultsImp::getVx()
{
	return vx;
}

std::vector<std::vector<double> > ResultsImp::getVy()
{
	return vy;
}

std::vector<std::vector<double> > ResultsImp::getVz()
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

std::vector<std::vector<double> > ResultsImp::getXNodes()
{
	return x;
}

std::vector<std::vector<double> > ResultsImp::getYNodes()
{
	return y;
}

std::vector<std::vector<double> > ResultsImp::getZNodes()
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

std::vector<int> ResultsImp::getTime()
{
	return time;
}

std::vector<std::vector<unsigned int> > ResultsImp::getLevels()
{
	return level;
}

std::vector<std::vector<double> > ResultsImp::getPress()
{
	return press;
}

std::vector<std::vector<double> > ResultsImp::getRho()
{
	return rho;
}

int ResultsImp::getL0()
{
	return l0;
}

bool ResultsImp::checkYourData()
{
	std::cout << "checking Simulation Results Data...";
	for (int i = 0; i < vx.size(); i++) {
		for (int j = 0; j < vx.at(i).size(); j++) {
			if (vx.at(i).at(j) != vx.at(i).at(j)) {
				std::cout << "done." << std::endl;
				std::cout << "Simulation Result Data contains failure data." << std::endl;
				std::cout << "Testing not possible." << std::endl;
				return false;
			}
			if (vy.at(i).at(j) != vy.at(i).at(j)) {
				std::cout << "done." << std::endl;
				std::cout << "Simulation Result Data contains failure data." << std::endl;
				std::cout << "Testing not possible." << std::endl;
				return false;
			}
			if (vz.at(i).at(j) != vz.at(i).at(j)) {
				std::cout << "done." << std::endl;
				std::cout << "Simulation Result Data contains failure data." << std::endl;
				std::cout << "Testing not possible." << std::endl;
				return false;
			}
			if (rho.at(i).at(j) != rho.at(i).at(j)) {
				std::cout << "done." << std::endl;
				std::cout << "Simulation Result Data contains failure data." << std::endl;
				std::cout << "Testing not possible." << std::endl;
				return false;
			}
			if (press.at(i).at(j) != press.at(i).at(j)) {
				std::cout << "done." << std::endl;
				std::cout << "Simulation Result Data contains failure data." << std::endl;
				std::cout << "Testing not possible." << std::endl;
				return false;
			}
		}
	}
	std::cout << "done." << std::endl;
	std::cout << "Simulation Result Data contains no failure data." << std::endl;
	return true;
}

ResultsImp::ResultsImp(int l0)
{
	this->l0 = l0;
}

ResultsImp::ResultsImp()
{

}

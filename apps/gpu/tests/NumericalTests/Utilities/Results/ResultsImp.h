#ifndef RESULTS_IMP_H
#define RESULTS_IMP_H

#include "Results.h"

class ResultsImp : public Results
{
public:
	int getNumberOfTimeSteps();
	int getTimeStepLength();
	std::vector<unsigned int> getTimeSteps();
	std::vector<int> getTime();
	std::vector<std::vector<double> > getVx();
	std::vector<std::vector<double> > getVy();
	std::vector<std::vector<double> > getVz();
	int getNumberOfXNodes();
	int getNumberOfYNodes();
	int getNumberOfZNodes();
	std::vector<std::vector<double> > getXNodes();
	std::vector<std::vector<double> > getYNodes();
	std::vector<std::vector<double> > getZNodes();
	std::vector<std::vector<unsigned int> > getLevels();
	std::vector<std::vector<double> > getPress();
	std::vector<std::vector<double> > getRho();
	int getL0();

	bool checkYourData();

protected:
	ResultsImp(int l0);
	ResultsImp() = default;

	unsigned int numberOfTimeSteps;
	unsigned int timeStepLength;
	unsigned int xNodes, yNodes, zNodes;
	unsigned int numberOfNodes;

	std::vector<unsigned int> timeStep;
	std::vector<int> time;
	std::vector<std::vector<double> > x, y, z;
	std::vector<std::vector<double> > vx, vy, vz;
	std::vector<std::vector<double> > press;
	std::vector<std::vector<double> > rho;
	std::vector<std::vector<unsigned int> > level;

	int l0;

private:
};
#endif
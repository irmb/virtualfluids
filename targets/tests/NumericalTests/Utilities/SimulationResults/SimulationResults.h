#ifndef SIMULATION_RESULTS_H
#define SIMULATION_RESULTS_H

#include <vector>
#include <memory>

class SimulationResults
{
public:
	static std::shared_ptr<SimulationResults> getNewInstance(unsigned int lx, unsigned int lz, unsigned int timeStepLength);
	void addTimeStep(unsigned int timeStep, unsigned int time, std::vector<double> x, std::vector<double> z, std::vector<double> vx, std::vector<double> vz, std::vector<double> press, std::vector<double> rho);
	int getNumberOfTimeSteps();
	std::vector<std::vector<double>> getVx();
	std::vector<std::vector<double>> getVz();
	int getXNodes();
	int getZNodes();
	int getTimeStepLength();

private:
	SimulationResults(unsigned int lx, unsigned int lz, unsigned int timeStepLength);

	unsigned int numberOfTimeSteps;
	unsigned int timeStepLength;
	unsigned int xNodes;
	unsigned int zNodes;
	unsigned int numberOfNodes;

	std::vector<unsigned int> timeStep;
	std::vector<unsigned int> time;
	std::vector<std::vector<double>> x, z;
	std::vector<std::vector<double>> vx, vz;
	std::vector<std::vector<double>> press;
	std::vector<std::vector<double>> rho;
};
#endif
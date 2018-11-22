#ifndef SIMULATION_RESULTS_H
#define SIMULATION_RESULTS_H

#include "../ResultsImp.h"
#include <memory>

class SimulationResults : public ResultsImp
{
public:
	static std::shared_ptr<SimulationResults> getNewInstance(unsigned int lx, unsigned int ly, unsigned int lz, unsigned int timeStepLength);
	void addTimeStep(unsigned int timeStep, unsigned int time, std::vector<unsigned int> level, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> vx, std::vector<double> vy, std::vector<double> vz, std::vector<double> press, std::vector<double> rho);

private:
	SimulationResults(unsigned int lx, unsigned int ly, unsigned int lz, unsigned int timeStepLength);

};
#endif
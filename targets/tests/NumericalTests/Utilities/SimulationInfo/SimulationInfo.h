#ifndef SIMULATION_INFO_H
#define SIMULATION_INFO_H

#include <string>

class SimulationInfo
{
public:
	virtual std::string getKernelName() = 0;
	virtual double getViscosity() = 0;
	virtual std::string getSimulationName() = 0;
	virtual std::string getSimulationParameterString() = 0;
	virtual int getLx() = 0;

private:

};
#endif 
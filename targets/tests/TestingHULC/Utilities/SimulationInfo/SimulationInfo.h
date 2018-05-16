#ifndef SIMULATION_INFO_H
#define SIMULATION_INFO_H

#include <iostream>

class SimulationInfo
{
public:
	virtual void makeSimulationHeadOutput() = 0;
	virtual void setStartTime() = 0;
	virtual void setEndTime() = 0;
	virtual std::string getSimulationRunTimeOutput() = 0;
private:

};
#endif
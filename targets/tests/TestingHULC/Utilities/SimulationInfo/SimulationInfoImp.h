#ifndef SIMULATION_INFO_IMP_H
#define SIMULATION_INFO_IMP_H

#include "SimulationInfo.h"

#include <memory>
#include <iostream>
#include <time.h>

class TestCout;

class SimulationInfoImp : public SimulationInfo
{
public:
	static std::shared_ptr<SimulationInfo> getNewInstance(std::shared_ptr<TestCout> output, std::string simName, int l);

	void makeSimulationHeadOutput();
	void setStartTime();
	void setEndTime();
	std::string getSimulationRunTimeOutput();

private:
	SimulationInfoImp(std::shared_ptr<TestCout> output, std::string simName, int l);
	SimulationInfoImp() {};
	double getSimTime();

	std::string simName;
	int l;
	std::shared_ptr<TestCout> output;
	time_t startTime, endTime;
};
#endif
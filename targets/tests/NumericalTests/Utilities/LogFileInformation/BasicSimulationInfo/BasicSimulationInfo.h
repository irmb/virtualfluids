#ifndef BASIC_SIMULATION_INFO_H
#define BASIC_SIMULATION_INFO_H

#include "../LogFileInformationImp.h"

#include <memory>

class BasicSimulationInfo : public LogFileInformationImp
{
public:
	static std::shared_ptr<BasicSimulationInfo> getNewInstance(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, std::string kernelName);
	std::string getOutput();

private:
	BasicSimulationInfo() {};
	BasicSimulationInfo(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, std::string kernelName);

	int numberOfTimeSteps;
	int basicTimeStepLength;
	double viscosity;
	std::string kernelName;
};
#endif
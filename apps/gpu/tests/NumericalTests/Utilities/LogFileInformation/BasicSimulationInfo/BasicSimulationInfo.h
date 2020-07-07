#ifndef BASIC_SIMULATION_INFO_H
#define BASIC_SIMULATION_INFO_H

#include "../LogFileInformationImp.h"

#include "VirtualFluids_GPU/Kernel//Utilities/KernelType.h"

#include <memory>

class BasicSimulationInfo : public LogFileInformationImp
{
public:
	static std::shared_ptr<BasicSimulationInfo> getNewInstance(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, KernelType kernel);
	std::string getOutput();

private:
	BasicSimulationInfo() {};
	BasicSimulationInfo(int numberOfTimeSteps, double viscosity, int basicTimeStepLength, KernelType kernel);

	int numberOfTimeSteps;
	int basicTimeStepLength;
	double viscosity;
	std::string kernelName;
};
#endif
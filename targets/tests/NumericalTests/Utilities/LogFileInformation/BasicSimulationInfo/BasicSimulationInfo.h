#ifndef BASIC_SIMULATION_INFO_H
#define BASIC_SIMULATION_INFO_H

#include "../LogFileInformationImp.h"

#include <memory>

class BasicSimulationInfo : public LogFileInformationImp
{
public:
	static std::shared_ptr<BasicSimulationInfo> getNewInstance(int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation, double viscosity);
	std::string getOutput();

private:
	BasicSimulationInfo() {};
	BasicSimulationInfo(int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation, double viscosity);

	int numberOfTimeSteps;
	int basisTimeStepLength;
	int startStepCalculation;
	double viscosity;
};
#endif
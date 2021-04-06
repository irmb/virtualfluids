#ifndef SIMULATION_INFO_IMP_H
#define SIMULATION_INFO_IMP_H

#include "SimulationInfo.h"

#include <memory>

class TimeInfo;

class SimulationInfoImp : public SimulationInfo
{
public:
	void setTimeInfo(std::shared_ptr<TimeInfo> timeInfo);

	std::string getKernelName();
	double getViscosity();
	std::string getSimulationName();
	std::string getSimulationParameterString();
	int getLx();
	int getNumberOfSimulations();
	int getSimulationID();
	std::string getRunTimeOutput();
	std::vector<std::string> getDataToCalcTests();

protected:
	SimulationInfoImp() {};
	SimulationInfoImp(int simID, std::string kernel, double viscosity, int lx, int numberOfSimulations, std::string simulationName, std::vector<std::string> dataToCalcTests);

	double viscosity;
	std::string kernelName;
	std::string simulationName;
	std::string simulationParameterString;
	int lx;
	int numberOfSimulations, simID;
	std::shared_ptr<TimeInfo> timeInfo;
	std::vector<std::string> dataToCalcTests;

private:

};
#endif 
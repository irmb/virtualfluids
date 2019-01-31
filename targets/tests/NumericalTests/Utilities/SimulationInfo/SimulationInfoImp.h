#ifndef SIMULATION_INFO_IMP_H
#define SIMULATION_INFO_IMP_H

#include "SimulationInfo.h"

class SimulationInfoImp : public SimulationInfo
{
public:
	std::string getKernelName();
	double getViscosity();
	std::string getSimulationName();
	std::string getSimulationParameterString();
	int getLx();
	int getNumberOfSimulations();
	int getSimulationID();

protected:
	SimulationInfoImp() {};
	SimulationInfoImp(int simID, std::string kernelName, double viscosity, int lx, int numberOfSimulations, std::string simulationName);

	double viscosity;
	std::string kernelName;
	std::string simulationName;
	std::string simulationParameterString;
	int lx;
	int numberOfSimulations, simID;

private:

};
#endif 
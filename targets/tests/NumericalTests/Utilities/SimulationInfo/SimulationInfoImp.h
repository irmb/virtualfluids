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

protected:
	SimulationInfoImp() {};
	SimulationInfoImp(int lx, double viscosity, std::string kernelName);

	double viscosity;
	std::string kernelName;
	std::string simulationName;
	std::string simulationParameterString;
	int lx;

private:

};
#endif 
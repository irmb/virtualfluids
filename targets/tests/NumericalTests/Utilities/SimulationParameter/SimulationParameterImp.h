#ifndef SIMULATION_PARAMETER_IMP_H
#define SIMULATION_PARAMETER_IMP_H

#include "SimulationParameter.h"

#include "LBM\LB.h"

struct GridInformationStruct;
struct BasicSimulationParameterStruct;

class SimulationParameterImp : public SimulationParameter
{
public:
	double getViscosity();
	std::string getGridPath();
	std::string getFilePath();
	unsigned int getNumberOfGridLevels();
	unsigned int getEndTime();
	unsigned int getTimeStepLength();
	unsigned int getLx();
	unsigned int getLz();
	std::vector< int> getDevices();
	std::shared_ptr< InitialCondition> getInitialCondition();
	std::shared_ptr< KernelConfiguration> getKernelConfiguration();

protected:
	SimulationParameterImp() {};
	SimulationParameterImp(std::string kernelName, double viscosity, std::shared_ptr<BasicSimulationParameterStruct> basicSimPara, std::shared_ptr<GridInformationStruct> gridInfo);

	void generateFileDirectionInMyStystem(std::string filePath);

	real viscosity;
	real lx, l0, lz;
	unsigned int numberOfTimeSteps, basisTimeStepLength;
	std::string gridPath;
	std::string filePath;
	std::vector<int> devices;

	unsigned int maxLevel, numberOfGridLevels;
	unsigned int timeStepLength;

	std::shared_ptr< InitialCondition> initialCondition;
	std::shared_ptr< KernelConfiguration> kernelConfig;
};

#endif

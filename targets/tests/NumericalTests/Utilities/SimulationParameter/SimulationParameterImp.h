#ifndef SIMULATION_PARAMETER_IMP_H
#define SIMULATION_PARAMETER_IMP_H

#include "SimulationParameter.h"

#include "LBM\LB.h"

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
	unsigned int getYSliceForCalculation();
	unsigned int getStartTimeCalculation();
	bool getWriteFiles();
	unsigned int getStartTimeDataWriter();
	std::vector< int> getDevices();
	std::shared_ptr< InitialCondition> getInitialCondition();
	std::shared_ptr< KernelConfiguration> getKernelConfiguration();

protected:
	SimulationParameterImp() {};
	SimulationParameterImp(std::string simName, real viscosity, real lx, real lz, real l0,
		real lForTimeStepLength, unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
		unsigned int startStepCalculation, unsigned int ySliceForCalculation,
		std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
		bool writeFiles, unsigned int startStepFileWriter,
		std::vector<int> devices);

	void generateFilePath(std::string filePath);

	real viscosity;
	real lx, l0, lz;
	unsigned int numberOfTimeSteps, basisTimeStepLength;
	unsigned int startStepCalculation, startStepFileWriter, ySliceForCalculation;
	std::string gridPath;
	std::string filePath;
	std::string simName;
	bool writeFiles;
	std::vector<int> devices;

	unsigned int maxLevel, numberOfGridLevels;
	unsigned int timeStepLength;
	unsigned int startTimeCalculation, startTimeDataWriter;
	unsigned int endTime;

	std::shared_ptr< InitialCondition> initialCondition;
	std::shared_ptr< KernelConfiguration> kernelConfig;
};

#endif

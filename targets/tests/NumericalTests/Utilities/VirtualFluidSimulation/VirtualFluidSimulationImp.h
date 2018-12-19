#ifndef VIRTUAL_FLUID_SIMULATION_IMP_H
#define VIRTUAL_FLUID_SIMULATION_IMP_H

#include "VirtualFluidSimulation.h"

#include <string>
#include <memory>

class InitialCondition;
class DataWriter;
class Parameter;
class GridProvider;
class KernelConfiguration;
class TestSimulation;

class VirtualFluidSimulationImp : public VirtualFluidSimulation
{
public:
	void run();

	static std::shared_ptr< VirtualFluidSimulationImp> getNewInstance();
	void initParameter(std::shared_ptr<Parameter> para, std::shared_ptr< KernelConfiguration> kernelConfig, real viscosity, std::string gridPath, std::string filePath, int numberOfGridLevels, unsigned int endTime, unsigned int timeStepLength, std::vector<int> devices, real velocity);
	void initInitialConditions(std::shared_ptr< InitialCondition> initialCondition);
	void initGridProvider();

	void setDataWriter(std::shared_ptr< DataWriter> dataWriter);
	void setTestSimulation(std::shared_ptr< TestSimulation> testSim);

protected:
	VirtualFluidSimulationImp() {};
		
private:
	std::shared_ptr< Parameter> para;
	std::shared_ptr< InitialCondition> initialCondition;
	std::shared_ptr< GridProvider> grid;
	std::shared_ptr<TestSimulation> testSim;
	std::shared_ptr<DataWriter> dataWriter;
};
#endif
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
class TimeTracking;
class NumericalTestSuite;

class VirtualFluidSimulationImp : public VirtualFluidSimulation
{
public:
	void run();

	static std::shared_ptr<VirtualFluidSimulationImp> getNewInstance();

	void setParameter(std::shared_ptr<Parameter> para);
	void setGridProvider(std::shared_ptr<GridProvider> grid);
	void setDataWriter(std::shared_ptr<DataWriter> dataWriter);
	void setNumericalTestSuite(std::shared_ptr<NumericalTestSuite> numericalTestSuite);
	void setTimeTracking(std::shared_ptr<TimeTracking> timeTracking);

protected:
	VirtualFluidSimulationImp() {};
		
private:
	std::shared_ptr<Parameter> para;
	std::shared_ptr<InitialCondition> initialCondition;
	std::shared_ptr<GridProvider> grid;
	std::shared_ptr<DataWriter> dataWriter;
	std::shared_ptr<NumericalTestSuite> numericalTestSuite;
	std::shared_ptr<TimeTracking> timeTracking;
};
#endif
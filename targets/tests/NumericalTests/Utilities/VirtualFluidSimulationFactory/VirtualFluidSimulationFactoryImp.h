#ifndef VIRTUAL_FLUID_SIMULATION_FACTORY_IMP_H
#define VIRTUAL_FLUID_SIMULATION_FACTORY_IMP_H

#include "VirtualFluidSimulationFactory.h"

class GridReaderforTesting;
class InitialCondition;
class Parameter;
class SimulationParameter;

class VirtualFluidSimulationFactoryImp: public VirtualFluidSimulationFactory
{
public:
	static std::shared_ptr< VirtualFluidSimulationFactory> getNewInstance();
	std::vector<std::shared_ptr<VirtualFluidSimulation> > makeVirtualFluidSimulations(std::vector<std::shared_ptr<TestSimulation> > testSim);

protected:
	VirtualFluidSimulationFactoryImp();
	
	std::shared_ptr<Parameter> makeParameter(std::shared_ptr<SimulationParameter> simPara);
	std::shared_ptr<GridReaderforTesting> makeGridReaderForTesting(std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<Parameter> para);
	void initInitialConditions(std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<Parameter> para);

private:

};
#endif
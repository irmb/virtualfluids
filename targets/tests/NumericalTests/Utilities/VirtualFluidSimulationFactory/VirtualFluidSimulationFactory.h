#ifndef VIRTUAL_FLUID_SIMULATION_FACTORY_H
#define VIRTUAL_FLUID_SIMULATION_FACTORY_H

#include <memory>
#include <vector>

class VirtualFluidSimulation;
class TestSimulation;

class VirtualFluidSimulationFactory
{
public:
	virtual std::vector< std::shared_ptr< VirtualFluidSimulation> > makeVirtualFluidSimulations(std::vector< std::shared_ptr< TestSimulation> > testSim) = 0;
};
#endif

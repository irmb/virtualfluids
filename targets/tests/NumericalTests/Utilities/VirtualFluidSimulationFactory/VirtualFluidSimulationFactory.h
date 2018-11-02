#ifndef VIRTUAL_FLUID_SIMULATION_FACTORY_H
#define VIRTUAL_FLUID_SIMULATION_FACTORY_H

#include <memory>
#include <vector>

class VirtualFluidSimulation;
class SimulationParameter;

class VirtualFluidSimulationFactory
{
public:
	virtual std::vector< std::shared_ptr< VirtualFluidSimulation> > makeVirtualFluidSimulations(std::vector< std::shared_ptr< SimulationParameter> > simPara) = 0;
};
#endif

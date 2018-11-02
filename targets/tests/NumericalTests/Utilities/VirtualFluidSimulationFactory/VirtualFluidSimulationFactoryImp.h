#ifndef VIRTUAL_FLUID_SIMULATION_FACTORY_IMP_H
#define VIRTUAL_FLUID_SIMULATION_FACTORY_IMP_H

#include "VirtualFluidSimulationFactory.h"


class VirtualFluidSimulationFactoryImp: public VirtualFluidSimulationFactory
{
public:
	static std::shared_ptr< VirtualFluidSimulationFactory> getNewInstance();
	std::vector< std::shared_ptr< VirtualFluidSimulation> > makeVirtualFluidSimulations(std::vector< std::shared_ptr< SimulationParameter> > simPara);

protected:
	VirtualFluidSimulationFactoryImp();

private:
	std::vector< std::shared_ptr< SimulationParameter> > simPara;
};
#endif
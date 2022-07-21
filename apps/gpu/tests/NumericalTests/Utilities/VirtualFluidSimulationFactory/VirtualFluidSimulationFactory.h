#ifndef VIRTUAL_FLUID_SIMULATION_FACTORY_IMP_H
#define VIRTUAL_FLUID_SIMULATION_FACTORY_IMP_H

#include <vector>
#include <memory>
#include "Utilities/VirtualFluidSimulation/VirtualFluidSimulation.h"
#include "Utilities/TestSimulation/TestSimulation.h"

namespace vf::gpu::tests
{
std::vector<std::shared_ptr<VirtualFluidSimulation>> makeVirtualFluidSimulations(std::vector<std::shared_ptr<TestSimulation>> testSim);
}

#endif
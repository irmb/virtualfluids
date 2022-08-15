#ifndef VIRTUAL_FLUID_SIMULATION_FACTORY_IMP_H
#define VIRTUAL_FLUID_SIMULATION_FACTORY_IMP_H

#include "Utilities/TestSimulation/TestSimulation.h"
#include <functional>
#include <memory>
#include <vector>

namespace vf::gpu::tests
{
std::shared_ptr<Parameter> makeParameter(std::shared_ptr<SimulationParameter> simPara);
const std::function<void()> makeVirtualFluidSimulation(std::shared_ptr<Parameter> para,
                                                       std::shared_ptr<InitialCondition> condition,
                                                       std::shared_ptr<DataWriter> dataWriter);
} // namespace vf::gpu::tests

#endif
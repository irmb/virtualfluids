#include "VirtualFluidSimulationImp.h"

#include "Utilities\TestSimulation\TestSimulation.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"

#include <sstream>

void VirtualFluidSimulationImp::run()
{
	testSim->makeSimulationHeadOutput();
	testSim->setSimulationStartTime();

	Simulation sim;
	sim.init(para, grid, dataWriter);
	sim.run();

	testSim->setSimulationEndTimeAndNotifyObserver();

	sim.free();
}

void VirtualFluidSimulationImp::setParameter(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

void VirtualFluidSimulationImp::setGridProvider(std::shared_ptr<GridProvider> grid)
{
	this->grid = grid;
}

std::shared_ptr<VirtualFluidSimulationImp> VirtualFluidSimulationImp::getNewInstance()
{
	return std::shared_ptr<VirtualFluidSimulationImp>(new VirtualFluidSimulationImp());
}

void VirtualFluidSimulationImp::setDataWriter(std::shared_ptr<DataWriter> dataWriter)
{
	this->dataWriter = dataWriter;
}

void VirtualFluidSimulationImp::setTestSimulation(std::shared_ptr<TestSimulation> testSim)
{
	this->testSim = testSim;
}

#include "VirtualFluidSimulation.h"

#include "Utilities/NumericalTestSuite/NumericalTestSuite.h"
#include "Utilities/Time/TimeTracking.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactory.h"

#include <sstream>

void VirtualFluidSimulation::run()
{
	numericalTestSuite->makeSimulationHeadOutput();
	BoundaryConditionFactory bc_factory;
	auto sim = Simulation(para, cudaManager, vf::gpu::Communicator::getInstance(), *grid.get(), &bc_factory);
	sim.setDataWriter(dataWriter);

	timeTracking->setSimulationStartTime();
	sim.run();
	timeTracking->setSimulationEndTime();

	numericalTestSuite->startPostProcessing();
}

void VirtualFluidSimulation::setParameter(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

void VirtualFluidSimulation::setCudaMemoryManager(std::shared_ptr<CudaMemoryManager> cudaManager)
{
	this->cudaManager = cudaManager;
}

void VirtualFluidSimulation::setGridProvider(std::shared_ptr<GridProvider> grid)
{
	this->grid = grid;
}

void VirtualFluidSimulation::setDataWriter(std::shared_ptr<DataWriter> dataWriter)
{
	this->dataWriter = dataWriter;
}

void VirtualFluidSimulation::setNumericalTestSuite(std::shared_ptr<NumericalTestSuite> numericalTestSuite)
{
	this->numericalTestSuite = numericalTestSuite;
}

void VirtualFluidSimulation::setTimeTracking(std::shared_ptr<TimeTracking> timeTracking)
{
	this->timeTracking = timeTracking;
}

void VirtualFluidSimulation::setKernelFactory(std::shared_ptr<KernelFactory> kernelFactory)
{
	this->kernelFactory = kernelFactory;
}

void VirtualFluidSimulation::setPreProcessorFactory(std::shared_ptr<PreProcessorFactory> preProcessorFactory)
{
	this->preProcessorFactory = preProcessorFactory;
}

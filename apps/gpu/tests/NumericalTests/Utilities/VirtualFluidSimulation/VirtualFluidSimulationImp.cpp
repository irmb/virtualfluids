#include "VirtualFluidSimulationImp.h"

#include "Utilities/NumericalTestSuite/NumericalTestSuite.h"
#include "Utilities/Time/TimeTracking.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactory.h"

#include <sstream>

void VirtualFluidSimulationImp::run()
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

void VirtualFluidSimulationImp::setParameter(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

void VirtualFluidSimulationImp::setCudaMemoryManager(std::shared_ptr<CudaMemoryManager> cudaManager)
{
	this->cudaManager = cudaManager;
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

void VirtualFluidSimulationImp::setNumericalTestSuite(std::shared_ptr<NumericalTestSuite> numericalTestSuite)
{
	this->numericalTestSuite = numericalTestSuite;
}

void VirtualFluidSimulationImp::setTimeTracking(std::shared_ptr<TimeTracking> timeTracking)
{
	this->timeTracking = timeTracking;
}

void VirtualFluidSimulationImp::setKernelFactory(std::shared_ptr<KernelFactory> kernelFactory)
{
	this->kernelFactory = kernelFactory;
}

void VirtualFluidSimulationImp::setPreProcessorFactory(std::shared_ptr<PreProcessorFactory> preProcessorFactory)
{
	this->preProcessorFactory = preProcessorFactory;
}

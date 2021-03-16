#include <memory>
#include <string>
#include <set>
#include <utility>
#include <cmath>
#include <omp.h>
#include <mpi.h>

#include <basics/utilities/UbScheduler.h>
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbSystem3D.h>

#include <BoundaryConditions/BCProcessor.h>
#include <CoProcessors/CoProcessor.h>
#include <CoProcessors/NUPSCounterCoProcessor.h>
#include <CoProcessors/WriteBlocksCoProcessor.h>
#include <CoProcessors/WriteBoundaryConditionsCoProcessor.h>
#include <CoProcessors/WriteMacroscopicQuantitiesCoProcessor.h>
#include <Grid/BasicCalculator.h>
#include <Grid/Calculator.h>
#include <Grid/Grid3D.h>
#include <Interactors/InteractorsHelper.h>
#include <LBM/CompressibleOffsetMomentsInterpolationProcessor.h>
#include <LBM/LBMKernel.h>
#include <LBM/LBMUnitConverter.h>
#include <Parallel/MPICommunicator.h>
#include <Visitors/GenBlocksGridVisitor.h>
#include <Visitors/InitDistributionsBlockVisitor.h>
#include <Visitors/MetisPartitioningGridVisitor.h>
#include <Visitors/SetConnectorsBlockVisitor.h>
#include <Visitors/SetKernelBlockVisitor.h>

#include <simulationconfig/SimulationParameters.h>
#include <simulationconfig/Simulation.h>


Simulation::Simulation()
{
    this->communicator = MPICommunicator::getInstance();
    this->grid = std::shared_ptr<Grid3D>(new Grid3D(communicator));
    this->interactors = std::vector<std::shared_ptr<Interactor3D>>();
    this->bcVisitor = BoundaryConditionsBlockVisitor();
    this->registeredAdapters = std::set<std::shared_ptr<BCAdapter>>();
}

void Simulation::setGridParameters(std::shared_ptr<GridParameters> parameters)
{
    this->gridParameters = std::move(parameters);
}

void Simulation::setPhysicalParameters(std::shared_ptr<PhysicalParameters> parameters)
{
    this->physicalParameters = std::move(parameters);
}

void Simulation::setRuntimeParameters(std::shared_ptr<RuntimeParameters> parameters)
{
    this->simulationParameters = std::move(parameters);
}

void
Simulation::addObject(const std::shared_ptr<GbObject3D> &object, const std::shared_ptr<BCAdapter> &bcAdapter, int state,
                      const std::string &folderPath)
{
    const bool is_in = registeredAdapters.find(bcAdapter) != registeredAdapters.end();
    if (!is_in) addBCAdapter(bcAdapter);
    this->interactors.push_back(lbmSystem->makeInteractor(object, this->grid, bcAdapter, state));
    GbSystem3D::writeGeoObject(object, writerConfig.outputPath + folderPath, writerConfig.getWriter());
}

void Simulation::addBCAdapter(const std::shared_ptr<BCAdapter> &bcAdapter)
{
    registeredAdapters.insert(bcAdapter);
    this->bcVisitor.addBC(bcAdapter);
}

void Simulation::setKernelConfiguration(const std::shared_ptr<LBMKernelConfiguration> &kernel)
{
    this->kernelConfig = kernel;
    this->lbmKernel = kernelFactory.makeKernel(kernel->kernelType);
    this->lbmSystem = kernelFactory.makeLBMSystem(kernel->kernelType);
}

void Simulation::setWriterConfiguration(const WriterConfiguration &config)
{
    this->writerConfig = config;
}

WriterConfiguration &Simulation::getWriterConfig()
{
    return writerConfig;
}

void Simulation::run()
{
    UBLOG(logINFO, "Beginning simulation setup for MPI rank " << communicator->getProcessID())
    grid->setDeltaX(gridParameters->nodeDistance);
    grid->setPeriodicX1(gridParameters->periodicBoundaryInX1);
    grid->setPeriodicX2(gridParameters->periodicBoundaryInX2);
    grid->setPeriodicX3(gridParameters->periodicBoundaryInX3);

    //int &numberOfNodesInReferenceDirection = gridParameters->numberOfNodesPerDirection[gridParameters->referenceDirectionIndex];
    std::shared_ptr<LBMUnitConverter> converter = makeLBMUnitConverter();

    int &nodesInX1 = gridParameters->numberOfNodesPerDirection[0];
    int &nodesInX2 = gridParameters->numberOfNodesPerDirection[1];
    int &nodesInX3 = gridParameters->numberOfNodesPerDirection[2];
    logSimulationData(nodesInX1, nodesInX2, nodesInX3);

    setBlockSize(nodesInX1, nodesInX2, nodesInX3);
    auto gridCube = makeSimulationBoundingBox();

    generateBlockGrid(gridCube);

    setKernelForcing(lbmKernel, converter);
    setBoundaryConditionProcessor(lbmKernel);

    auto metisVisitor = std::make_shared<MetisPartitioningGridVisitor>(communicator,
                                                                       MetisPartitioningGridVisitor::LevelBased,
                                                                       D3Q27System::B, MetisPartitioner::RECURSIVE);

    InteractorsHelper intHelper(grid, metisVisitor);
    for (auto const &interactor : interactors)
        intHelper.addInteractor(interactor);

    intHelper.selectBlocks();


    int numberOfProcesses = communicator->getNumberOfProcesses();
    SetKernelBlockVisitor kernelVisitor(lbmKernel, physicalParameters->latticeViscosity,
                                        numberOfProcesses);
    grid->accept(kernelVisitor);
    intHelper.setBC();

    //double bulkViscosity = physicalParameters->latticeViscosity * physicalParameters->bulkViscosityFactor;
    //auto iProcessor = std::make_shared<CompressibleOffsetMomentsInterpolationProcessor>();
    //iProcessor->setBulkViscosity(physicalParameters->latticeViscosity, bulkViscosity);

    //SetConnectorsBlockVisitor setConnsVisitor(communicator, true,
    //                                          lbmSystem->getNumberOfDirections(),
    //                                          physicalParameters->latticeViscosity, iProcessor);

    OneDistributionSetConnectorsBlockVisitor setConnsVisitor(communicator);
    grid->accept(setConnsVisitor);

    InitDistributionsBlockVisitor initVisitor;
    grid->accept(initVisitor);
    grid->accept(setConnsVisitor);
    grid->accept(bcVisitor);

    writeBoundaryConditions();
    // important: run this after metis & intHelper.selectBlocks()
    writeBlocksToFile();

    auto visualizationScheduler = std::make_shared<UbScheduler>(simulationParameters->timeStepLogInterval);
    auto mqCoProcessor = makeMacroscopicQuantitiesCoProcessor(converter,
                                                              visualizationScheduler);

    std::shared_ptr<UbScheduler> nupsScheduler(new UbScheduler(100, 100));
    std::shared_ptr<CoProcessor> nupsCoProcessor(
            new NUPSCounterCoProcessor(grid, nupsScheduler, simulationParameters->numberOfThreads, communicator));


#ifdef _OPENMP
    omp_set_num_threads(simulationParameters->numberOfThreads);
    UBLOG(logINFO, "OpenMP is set to run with " << omp_get_num_threads() << " threads")
#endif

    auto calculator = std::make_shared<BasicCalculator>(grid, visualizationScheduler,
                                                        simulationParameters->numberOfTimeSteps);
    calculator->addCoProcessor(nupsCoProcessor);
    calculator->addCoProcessor(mqCoProcessor);

    UBLOG(logINFO, "Simulation-start")
    calculator->calculate();
    UBLOG(logINFO, "Simulation-end")
}

void
Simulation::setKernelForcing(const std::shared_ptr<LBMKernel> &kernel,
                             std::shared_ptr<LBMUnitConverter> &converter) const
{
    kernel->setWithForcing(kernelConfig->useForcing);
    kernel->setForcingX1(kernelConfig->forcingX1 * converter->getFactorForceWToLb());
    kernel->setForcingX2(kernelConfig->forcingX2 * converter->getFactorForceWToLb());
    kernel->setForcingX3(kernelConfig->forcingX3 * converter->getFactorForceWToLb());
}

void Simulation::logSimulationData(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const
{
    UBLOG(logINFO, "Domain size = " << nodesInX1 << " x " << nodesInX2 << " x " << nodesInX3)
    UBLOG(logINFO, "dx          = " << gridParameters->nodeDistance << " m")
    UBLOG(logINFO, "latticeViscosity    = " << physicalParameters->latticeViscosity)
}

void Simulation::generateBlockGrid(const std::shared_ptr<GbObject3D> &gridCube) const
{
    UBLOG(logINFO, "Generate block grid")
    GenBlocksGridVisitor genBlocks(gridCube);
    grid->accept(genBlocks);
}

void Simulation::setBoundaryConditionProcessor(const std::shared_ptr<LBMKernel> &kernel)
{
    UBLOG(logINFO, "Create boundary conditions processor")
    auto bcProc = std::make_shared<BCProcessor>();
    kernel->setBCProcessor(bcProc);
}

void Simulation::setBlockSize(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const
{
    int blockSizeX1 = nodesInX1 / gridParameters->blocksPerDirection[0];
    int blockSizeX2 = nodesInX2 / gridParameters->blocksPerDirection[1];
    int blockSizeX3 = nodesInX3 / gridParameters->blocksPerDirection[2];
    UBLOG(logINFO, "Block size  = " << blockSizeX1 << " x " << blockSizeX2 << " x " << blockSizeX3)
    grid->setBlockNX(blockSizeX1,
                     blockSizeX2,
                     blockSizeX3);
}

std::shared_ptr<LBMUnitConverter>
Simulation::makeLBMUnitConverter()
{
    return std::make_shared<LBMUnitConverter>();
}

std::shared_ptr<CoProcessor>
Simulation::makeMacroscopicQuantitiesCoProcessor(const std::shared_ptr<LBMUnitConverter> &converter,
                                                 const std::shared_ptr<UbScheduler> &visualizationScheduler) const
{
    auto mqCoProcessor = std::make_shared<WriteMacroscopicQuantitiesCoProcessor>(grid, visualizationScheduler,
                                                                                 writerConfig.outputPath,
                                                                                 writerConfig.getWriter(),
                                                                                 converter,
                                                                                 communicator);
    mqCoProcessor->process(0);
    return mqCoProcessor;
}

void Simulation::writeBoundaryConditions() const
{
    auto geoSch = std::make_shared<UbScheduler>(1);
    WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, writerConfig.outputPath, writerConfig.getWriter(),
                                             communicator);
    ppgeo.process(0);
}

void Simulation::writeBlocksToFile() const
{
    UBLOG(logINFO, "Write block grid to VTK-file")
    auto ppblocks = std::make_shared<WriteBlocksCoProcessor>(grid,
                                                             std::make_shared<UbScheduler>(1),
                                                             writerConfig.outputPath,
                                                             writerConfig.getWriter(),
                                                             communicator);
    ppblocks->process(0);
    ppblocks.reset();
}

std::shared_ptr<GbObject3D>
Simulation::makeSimulationBoundingBox() const
{
    auto box = gridParameters->boundingBox();

    UBLOG(logINFO, "Bounding box dimensions = [("
            << box->minX1 << ", " << box->minX2 << ", " << box->minX3 << "); ("
            << box->maxX1 << ", " << box->maxX2 << ", " << box->maxX3 << ")]")

    auto gridCube = std::make_shared<GbCuboid3D>(box->minX1, box->minX2, box->minX3, box->maxX1, box->maxX2, box->maxX3);
    GbSystem3D::writeGeoObject(gridCube.get(), writerConfig.outputPath + "/geo/gridCube", writerConfig.getWriter());
    return gridCube;
}

Simulation::~Simulation() = default;

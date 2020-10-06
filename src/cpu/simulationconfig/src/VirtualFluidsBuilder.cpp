#include <string>
#include <set>
#include <utility>
#include <cmath>
#include <omp.h>

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

#include <VirtualFluidsBuilder/SimulationParameters.h>
#include <VirtualFluidsBuilder/VirtualFluidsBuilder.h>


VirtualFluidsBuilder::VirtualFluidsBuilder()
{
    this->communicator = MPICommunicator::getInstance();
    this->grid = SPtr<Grid3D>(new Grid3D(communicator));
    this->interactors = std::vector<SPtr<Interactor3D>>();
    this->bcVisitor = BoundaryConditionsBlockVisitor();
    this->registeredAdapters = std::set<SPtr<BCAdapter>>();
}

void VirtualFluidsBuilder::setGridParameters(SPtr<GridParameters> parameters)
{
    this->gridParameters = std::move(parameters);
}

void VirtualFluidsBuilder::setPhysicalParameters(SPtr<PhysicalParameters> parameters)
{
    this->physicalParameters = std::move(parameters);
}

void VirtualFluidsBuilder::setSimulationParameters(SPtr<SimulationParameters> parameters)
{
    this->simulationParameters = std::move(parameters);
}

void
VirtualFluidsBuilder::addObject(const SPtr<GbObject3D> &object, const SPtr<BCAdapter> &bcAdapter, int state,
                                const std::string &folderPath)
{
    const bool is_in = registeredAdapters.find(bcAdapter) != registeredAdapters.end();
    if (!is_in) addBCAdapter(bcAdapter);
    this->interactors.push_back(lbmSystem->makeInteractor(object, this->grid, bcAdapter, state));
    GbSystem3D::writeGeoObject(object, writerConfig.outputPath + folderPath, writerConfig.getWriter());
}

void VirtualFluidsBuilder::addBCAdapter(const SPtr<BCAdapter> &bcAdapter)
{
    registeredAdapters.insert(bcAdapter);
    this->bcVisitor.addBC(bcAdapter);
}

void VirtualFluidsBuilder::setKernelConfig(const SPtr<LBMKernelConfig> &kernel)
{
    this->kernelConfig = kernel;
    this->lbmKernel = kernelFactory.makeKernel(kernel->kernelType);
    this->lbmSystem = kernelFactory.makeLBMSystem(kernel->kernelType);
}

void VirtualFluidsBuilder::setWriterConfig(const WriterConfig &config)
{
    this->writerConfig = config;
}

WriterConfig &VirtualFluidsBuilder::getWriterConfig()
{
    return writerConfig;
}

void VirtualFluidsBuilder::run()
{
    UBLOG(logINFO, "Beginning simulation setup")
    grid->setDeltaX(gridParameters->deltaX);
    grid->setPeriodicX1(gridParameters->periodicBoundaryInX1);
    grid->setPeriodicX2(gridParameters->periodicBoundaryInX2);
    grid->setPeriodicX3(gridParameters->periodicBoundaryInX3);

    int &numberOfNodesInReferenceDirection = gridParameters->numberOfNodesPerDirection[gridParameters->referenceDirectionIndex];
    std::shared_ptr<LBMUnitConverter> converter = makeLBMUnitConverter();

    int &nodesInX1 = gridParameters->numberOfNodesPerDirection[0];
    int &nodesInX2 = gridParameters->numberOfNodesPerDirection[1];
    int &nodesInX3 = gridParameters->numberOfNodesPerDirection[2];
    logSimulationData(nodesInX1, nodesInX2, nodesInX3);

    setBlockSize(nodesInX1, nodesInX2, nodesInX3);
    SPtr<GbObject3D> gridCube = makeSimulationBoundingBox(nodesInX1, nodesInX2, nodesInX3);

    generateBlockGrid(gridCube);

    setKernelForcing(lbmKernel, converter);
    setBoundaryConditionProcessor(lbmKernel);

    SPtr<Grid3DVisitor> metisVisitor(
            new MetisPartitioningGridVisitor(communicator,
                                             MetisPartitioningGridVisitor::LevelBased,
                                             D3Q27System::B));

    InteractorsHelper intHelper(grid, metisVisitor);
    for (auto const &interactor : interactors)
        intHelper.addInteractor(interactor);

    intHelper.selectBlocks();

    int numberOfProcesses = communicator->getNumberOfProcesses();
    SetKernelBlockVisitor kernelVisitor(lbmKernel, physicalParameters->latticeViscosity,
                                        numberOfProcesses);
    grid->accept(kernelVisitor);
    intHelper.setBC();

    SPtr<InterpolationProcessor> iProcessor(new CompressibleOffsetMomentsInterpolationProcessor());
    dynamicPointerCast<CompressibleOffsetMomentsInterpolationProcessor>(iProcessor)->setBulkViscosity(
            physicalParameters->latticeViscosity,
            physicalParameters->latticeViscosity * physicalParameters->bulkViscosityFactor);

    SetConnectorsBlockVisitor setConnsVisitor(communicator, true,
                                              lbmSystem->getNumberOfDirections(),
                                              physicalParameters->latticeViscosity, iProcessor);

    InitDistributionsBlockVisitor initVisitor;
    grid->accept(initVisitor);
    grid->accept(setConnsVisitor);
    grid->accept(bcVisitor);

    writeBoundaryConditions();

    SPtr<UbScheduler> visualizationScheduler(new UbScheduler(simulationParameters->timeStepLogInterval));
    SPtr<CoProcessor> mqCoProcessor = makeMacroscopicQuantitiesCoProcessor(converter, visualizationScheduler);

    SPtr<UbScheduler> nupsScheduler(new UbScheduler(100, 100));
    SPtr<CoProcessor> nupsCoProcessor(
            new NUPSCounterCoProcessor(grid, nupsScheduler, simulationParameters->numberOfThreads, communicator));


#ifdef _OPENMP
    omp_set_num_threads(simulationParameters->numberOfThreads);
#endif

    SPtr<Calculator> calculator(
            new BasicCalculator(grid, visualizationScheduler, simulationParameters->numberOfTimeSteps));
    calculator->addCoProcessor(nupsCoProcessor);
    calculator->addCoProcessor(mqCoProcessor);

    UBLOG(logINFO, "Simulation-start")
    calculator->calculate();
    UBLOG(logINFO, "Simulation-end")
}

void
VirtualFluidsBuilder::setKernelForcing(const SPtr<LBMKernel> &kernel,
                                       std::shared_ptr<LBMUnitConverter> &converter) const
{
    kernel->setWithForcing(kernelConfig->useForcing);
    kernel->setForcingX1(kernelConfig->forcingX1 * converter->getFactorForceWToLb());
    kernel->setForcingX2(kernelConfig->forcingX2 * converter->getFactorForceWToLb());
    kernel->setForcingX3(kernelConfig->forcingX3 * converter->getFactorForceWToLb());
}

void VirtualFluidsBuilder::logSimulationData(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const
{
    UBLOG(logINFO, "Domain size = " << nodesInX1 << " x " << nodesInX2 << " x " << nodesInX3)
    UBLOG(logINFO, "dx          = " << gridParameters->deltaX << " m")
    UBLOG(logINFO, "latticeViscosity    = " << physicalParameters->latticeViscosity)
}

void VirtualFluidsBuilder::generateBlockGrid(const SPtr<GbObject3D> &gridCube) const
{
    UBLOG(logINFO, "Generate block grid")
    GenBlocksGridVisitor genBlocks(gridCube);
    grid->accept(genBlocks);
    writeBlocks();
}

void VirtualFluidsBuilder::setBoundaryConditionProcessor(const SPtr<LBMKernel> &kernel)
{
    UBLOG(logINFO, "Create boundary conditions processor")
    SPtr<BCProcessor> bcProc;
    bcProc = SPtr<BCProcessor>(new BCProcessor());
    kernel->setBCProcessor(bcProc);
}

void VirtualFluidsBuilder::setBlockSize(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const
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
VirtualFluidsBuilder::makeLBMUnitConverter()
{
    return SPtr<LBMUnitConverter>(new LBMUnitConverter());
}

SPtr<CoProcessor>
VirtualFluidsBuilder::makeMacroscopicQuantitiesCoProcessor(const std::shared_ptr<LBMUnitConverter> &converter,
                                                           const SPtr<UbScheduler> &visSch) const
{
    SPtr<CoProcessor> mqCoProcessor(
            new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, writerConfig.outputPath, writerConfig.getWriter(),
                                                      converter,
                                                      communicator));
    mqCoProcessor->process(0);
    return mqCoProcessor;
}

void VirtualFluidsBuilder::writeBoundaryConditions() const
{
    SPtr<UbScheduler> geoSch(new UbScheduler(1));
    WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, writerConfig.outputPath, writerConfig.getWriter(),
                                             communicator);
    ppgeo.process(0);
}

void VirtualFluidsBuilder::writeBlocks() const
{
    UBLOG(logINFO, "Write block grid to VTK-file")
    SPtr<CoProcessor> ppblocks(
            new WriteBlocksCoProcessor(grid,
                                       SPtr<UbScheduler>(new UbScheduler(1)),
                                       writerConfig.outputPath,
                                       writerConfig.getWriter(),
                                       communicator));

    ppblocks->process(0);
    ppblocks.reset();
}

SPtr<GbObject3D>
VirtualFluidsBuilder::makeSimulationBoundingBox(const int &nodesInX1, const int &nodesInX2,
                                                const int &nodesInX3) const
{
    double minX1 = 0, minX2 = 0, minX3 = 0;
    const double maxX1 = minX1 + gridParameters->deltaX * nodesInX1;
    const double maxX2 = minX2 + gridParameters->deltaX * nodesInX2;
    const double maxX3 = minX3 + gridParameters->deltaX * nodesInX3;
    UBLOG(logINFO, "Bounding box dimensions = [("
            << minX1 << ", " << minX2 << ", " << minX3 << "); ("
            << maxX1 << ", " << maxX2 << ", " << maxX3 << ")]")


    SPtr<GbObject3D> gridCube(new GbCuboid3D(minX1, minX2, minX3, maxX1, maxX2, maxX3));
    GbSystem3D::writeGeoObject(gridCube.get(), writerConfig.outputPath + "/geo/gridCube", writerConfig.getWriter());
    return gridCube;
}

VirtualFluidsBuilder::~VirtualFluidsBuilder() = default;

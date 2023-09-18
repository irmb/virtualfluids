#include <memory.h>
#include <string.h>
#include <basics/geometry3d/GbCuboid3D.h>
#include <cpu/VirtualFluids.h>

int main()
{
    const std::string outputPath = "./output";
    const int numberOfThreads = 1;
    const int timeStepLogInterval = 1000;
    const int numberOfTimeSteps = 10000;

    const double latticeViscosity = 5e-3;
    const double bulkViscosityFactor = 1.;
    const int deltaX = 1;
    const int nodesInX1 = 2;
    const int nodesInX2 = 2;
    const int nodesInX3 = 10;

    const double minX1 = 0, minX2 = 0, minX3 = 0;
    const double maxX1 = minX1 + nodesInX1 * deltaX;
    const double maxX2 = minX2 + nodesInX2 * deltaX;
    const double maxX3 = minX3 + nodesInX3 * deltaX;

    const auto lbmUnitConverter = std::make_shared<LBMUnitConverter>();
    const auto writer = WbWriterVtkXmlBinary::getInstance();

    const auto communicator = vf::parallel::MPICommunicator::getInstance();
    const auto kernel = std::make_shared<CompressibleCumulant4thOrderViscosityLBMKernel>();
    kernel->setBCProcessor(std::make_shared<BCProcessor>());
    kernel->setForcingX1(1e-6 * lbmUnitConverter->getFactorForceWToLb());

    auto grid = std::make_shared<Grid3D>();
    grid->setDeltaX(deltaX);
    grid->setPeriodicX1(true);
    grid->setPeriodicX2(true);
    grid->setPeriodicX3(false);
    grid->setBlockNX(nodesInX1, nodesInX2, nodesInX3);

    auto gridCube = std::make_shared<GbCuboid3D>(minX1, minX2, minX3, maxX1, maxX2, maxX3);
    GbSystem3D::writeGeoObject(gridCube.get(), outputPath + "/geo/gridCube", writer);

    GenBlocksGridVisitor genBlocks(gridCube);
    grid->accept(genBlocks);

    const auto scheduler(std::make_shared<UbScheduler>(1));
    auto ppblocks = std::make_shared<WriteBlocksCoProcessor>(grid, scheduler, outputPath, writer, communicator);
    ppblocks->process(0);
    ppblocks.reset();

    const auto metisVisitor(std::make_shared<MetisPartitioningGridVisitor>(communicator,
                                                                           MetisPartitioningGridVisitor::LevelBased,
                                                                           D3Q27System::B));

    const double blockLength = 3 * deltaX;
    const auto topWall = std::make_shared<GbCuboid3D>(minX1 - blockLength, minX2 - blockLength, minX3 - blockLength,
                                                      maxX1 + blockLength, maxX2 + blockLength, minX3);

    const auto bottomWall = std::make_shared<GbCuboid3D>(minX1 - blockLength, minX2 - blockLength, minX3,
                                                         maxX1 + blockLength, maxX2 + blockLength, minX3 + blockLength);

    const auto adapter = std::make_shared<NoSlipBCAdapter>();
    adapter->setBcAlgorithm(std::make_shared<NoSlipBCAlgorithm>());

    const auto topWallInteractor = std::make_shared<D3Q27Interactor>(topWall, grid, adapter, Interactor3D::SOLID);
    const auto bottomWallInteractor = std::make_shared<D3Q27Interactor>(bottomWall, grid, adapter, Interactor3D::SOLID);

    GbSystem3D::writeGeoObject(topWall, outputPath + "/geo/topWall", writer);
    GbSystem3D::writeGeoObject(bottomWall, outputPath + "/geo/bottomWall", writer);


    InteractorsHelper interactorsHelper(grid, metisVisitor);
    interactorsHelper.addInteractor(topWallInteractor);
    interactorsHelper.addInteractor(bottomWallInteractor);
    interactorsHelper.selectBlocks();

    int numberOfProcesses = communicator->getNumberOfProcesses();
    SetKernelBlockVisitor kernelVisitor(kernel, latticeViscosity, numberOfProcesses);
    grid->accept(kernelVisitor);
    interactorsHelper.setBC();

    const auto interpolationProcessor(std::make_shared<CompressibleOffsetMomentsInterpolator>());
    interpolationProcessor->setBulkViscosity(latticeViscosity, latticeViscosity * bulkViscosityFactor);

    SetConnectorsBlockVisitor setConnsVisitor(communicator,
                                              true,
                                              D3Q27System::ENDDIR,
                                              latticeViscosity,
                                              interpolationProcessor);


    InitDistributionsBlockVisitor initVisitor{};
    BoundaryConditionsBlockVisitor boundaryConditionsBlockVisitor{};
    grid->accept(initVisitor);
    grid->accept(setConnsVisitor);
    grid->accept(boundaryConditionsBlockVisitor);

    const auto geoSch(std::make_shared<UbScheduler>(1));
    WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, outputPath, writer, communicator);
    ppgeo.process(0);

    const auto visualizationScheduler(std::make_shared<UbScheduler>(timeStepLogInterval));
    const auto mqCoProcessor(
            std::make_shared<WriteMacroscopicQuantitiesCoProcessor>(grid, visualizationScheduler, outputPath, writer,
                                                                    lbmUnitConverter,
                                                                    communicator));
    mqCoProcessor->process(0);

    const auto nupsScheduler(std::make_shared<UbScheduler>(100, 100));
    const auto nupsCoProcessor(
            std::make_shared<NUPSCounterCoProcessor>(grid, nupsScheduler, numberOfThreads,
                                                     communicator));

#ifdef _OPENMP
    omp_set_num_threads(numberOfThreads);
#endif

    const auto calculator(std::make_shared<BasicCalculator>(grid, visualizationScheduler, numberOfTimeSteps));
    calculator->addCoProcessor(nupsCoProcessor);
    calculator->addCoProcessor(mqCoProcessor);
    calculator->calculate();

    return 0;
}
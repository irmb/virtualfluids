#define _USE_MATH_DEFINES
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

//////////////////////////////////////////////////////////////////////////

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <basics/StringUtilities/StringUtil.h>
#include <basics/config/ConfigurationFile.h>
#include <logger/Logger.h>
#include <parallel/MPICommunicator.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/geometries/Cuboid/Cuboid.h"
#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/MultipleGridBuilderFacade.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"
#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"
#include "GridGenerator/utilities/communication.h"

//////////////////////////////////////////////////////////////////////////

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/Kernel/KernelFactory/KernelFactoryImp.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/LBM/Simulation.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"

void runVirtualFluids(const vf::basics::ConfigurationFile& config)
{
    vf::parallel::Communicator& communicator = *vf::parallel::MPICommunicator::getInstance();
    const auto numberOfProcesses = communicator.getNumberOfProcesses();
    const auto processID = communicator.getProcessID();
    SPtr<Parameter> para = std::make_shared<Parameter>(numberOfProcesses, processID, &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
    GridScalingFactory scalingFactory = GridScalingFactory();

    // configure simulation parameters

    const bool useLevels = true;

    const std::string outPath("output/" + std::to_string(para->getNumprocs()) + "GPU/");
    const std::string gridPath = "output/";
    std::string simulationName("SphereMultiGPU");

    para->useReducedCommunicationAfterFtoC = para->getNumprocs() != 1;

    const real length = 1.0;
    const real reynoldsNumber = 1000.0;
    const real velocity = 1.0;
    const real velocityLB = 0.05;
    const uint numberOfNodesX = 80;

    // compute  parameters in lattcie units

    const real dxGrid = length / real(numberOfNodesX);
    const real deltaT = velocityLB / velocity * dxGrid;
    const real viscosityLB = numberOfNodesX * velocityLB / reynoldsNumber;

    // set parameters

    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(velocity / velocityLB);
    para->setDensityRatio((real)1.0);

    para->setOutputPath(outPath);
    para->setOutputPrefix(simulationName);
    para->setPrintFiles(true);

    para->configureMainKernel(vf::collisionKernel::compressible::K17CompressibleNavierStokes);
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

    vf::logging::Logger::changeLogPath(outPath + "vflog_process" + std::to_string(processID) );
    vf::logging::Logger::initializeLogger();

    // log simulation parameters

    VF_LOG_INFO("LB parameters:");
    VF_LOG_INFO("velocity LB [dx/dt]              = {}", velocityLB);
    VF_LOG_INFO("viscosity LB [dx/dt]             = {}", viscosityLB);
    VF_LOG_INFO("dxGrid [-]                       = {}\n", dxGrid);
    VF_LOG_INFO("deltaT [s]                       = {}", deltaT);
    VF_LOG_INFO("simulation parameters:");
    VF_LOG_INFO("mainKernel                       = {}\n", para->getMainKernel());

    // configure simulation grid

    UPtr<MultipleGridBuilderFacade> gridBuilderFacade;

    auto domainDimensions = std::make_shared<GridDimensions>();
    domainDimensions->minX = -0.5 * length;
    domainDimensions->maxX = 0.5 * length;
    domainDimensions->minY = -0.5 * length;
    domainDimensions->maxY = 0.5 * length;
    domainDimensions->minZ = -0.5 * length;
    domainDimensions->maxZ = 0.5 * length;
    domainDimensions->delta = dxGrid;

    gridBuilderFacade = std::make_unique<MultipleGridBuilderFacade>(std::move(domainDimensions), 8. * dxGrid);

    gridBuilderFacade->addGeometry(std::make_shared<Sphere>(0.0, 0.0, 0.0, 0.1 * length));

    std::shared_ptr<Object> level1 = nullptr;
    if (useLevels) {
        gridBuilderFacade->setNumberOfLayersForRefinement(10, 8);
        level1 = std::make_shared<Sphere>(0.0, 0.0, 0.0, 0.25 * length);
        gridBuilderFacade->addFineGrid(level1, 1);
    }

    // configure subdomains for simulation on multiple gpus

    const real xSplit = 0.0;
    const real ySplit = 0.0;
    const real zSplit = 0.0;

    if (numberOfProcesses == 2) {
        gridBuilderFacade->addDomainSplit(zSplit, Axis::z);
    } else if (numberOfProcesses == 4) {
        gridBuilderFacade->addDomainSplit(xSplit, Axis::y);
        gridBuilderFacade->addDomainSplit(zSplit, Axis::z);
    } else if (numberOfProcesses == 8) {
        gridBuilderFacade->addDomainSplit(xSplit, Axis::x);
        gridBuilderFacade->addDomainSplit(ySplit, Axis::y);
        gridBuilderFacade->addDomainSplit(zSplit, Axis::z);
    }

    // create grids
    gridBuilderFacade->createGrids(processID);

    // configure boundary conditions

    // call after createGrids()
    gridBuilderFacade->setPeriodicBoundaryCondition(false, false, false);
    gridBuilderFacade->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);
    gridBuilderFacade->setVelocityBoundaryCondition(SideType::MY, velocityLB, 0.0, 0.0);
    gridBuilderFacade->setVelocityBoundaryCondition(SideType::PY, velocityLB, 0.0, 0.0);
    gridBuilderFacade->setVelocityBoundaryCondition(SideType::MZ, velocityLB, 0.0, 0.0);
    gridBuilderFacade->setVelocityBoundaryCondition(SideType::PZ, velocityLB, 0.0, 0.0);
    gridBuilderFacade->setPressureBoundaryCondition(SideType::PX, 0.0);

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

    // move grid from grid generator to simulation

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilderFacade->getGridBuilder(), para, cudaMemoryManager, communicator);
    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, &scalingFactory);

    // run simulation
    sim.run();
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    std::string str, str2, configFile;

    if (argv != NULL) {

        try {
            VF_LOG_TRACE("For the default config path to work, execute the app from the project root.");
            vf::basics::ConfigurationFile config = vf::basics::loadConfig(argc, argv, "./apps/gpu/SphereMultiGPU/sphere_1gpu.cfg");
            runVirtualFluids(config);

            //////////////////////////////////////////////////////////////////////////
        } catch (const spdlog::spdlog_ex &ex) {
            std::cout << "Log initialization failed: " << ex.what() << std::endl;
        } catch (const std::bad_alloc &e) {
            VF_LOG_CRITICAL("Bad Alloc: {}", e.what());
        } catch (const std::exception &e) {
            VF_LOG_CRITICAL("exception: {}", e.what());
        } catch (...) {
            VF_LOG_CRITICAL("Unknown exception!");
        }
    }

    MPI_Finalize();
    return 0;
}

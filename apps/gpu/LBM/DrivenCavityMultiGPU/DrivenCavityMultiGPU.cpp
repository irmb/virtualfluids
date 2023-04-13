#define _USE_MATH_DEFINES
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <filesystem>

#include "mpi.h"

//////////////////////////////////////////////////////////////////////////

#include "basics/Core/DataTypes.h"
#include "basics/Core/VectorTypes.h"
#include "basics/PointerDefinitions.h"

#include "basics/Core/Logger/Logger.h"
#include "basics/Core/StringUtilities/StringUtil.h"
#include "basics/config/ConfigurationFile.h"
#include "logger/Logger.h"

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "geometries/Cuboid/Cuboid.h"
#include "geometries/TriangularMesh/TriangularMesh.h"

#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"
#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/Factories/GridScalingFactory.h"
#include "VirtualFluids_GPU/Kernel/Utilities/KernelTypes.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

//////////////////////////////////////////////////////////////////////////

#include "utilities/communication.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(std::filesystem::path& configPath)
{
    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    vf::gpu::Communicator& communicator = vf::gpu::Communicator::getInstance();

    auto gridFactory = GridFactory::make();
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    vf::basics::ConfigurationFile config;
    std::cout << configPath << std::endl;
    config.load(configPath.string());
    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNumberOfProcess(), communicator.getPID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
    GridScalingFactory scalingFactory = GridScalingFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool useGridGenerator = true;
    bool useLevels        = true;

    if (para->getNumprocs() == 1) {
        para->useReducedCommunicationAfterFtoC = false;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const std::string outPath("output/" + std::to_string(para->getNumprocs()) + "GPU/");
    const std::string gridPath = "output/";
    std::string simulationName("DrivenCavityMultiGPU");

    const real L = 1.0;
    const real Re = 1000.0;
    const real velocity = 1.0;
    const real velocityLB = 0.05; // LB units
    const uint nx = 64;

    // para->setTimestepOut(10000);   // set in config
    // para->setTimestepEnd(10000);   // set in config

    const real dxGrid      = L / real(nx);
    const real dt  = velocityLB / velocity * dxGrid;
    const real vxLB        = velocityLB / (real)sqrt(2.0); // LB units
    const real vyLB        = velocityLB / (real)sqrt(2.0); // LB units
    const real viscosityLB = nx * velocityLB / Re;         // LB units

    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(velocity / velocityLB);
    para->setDensityRatio((real)1.0);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (para->getOutputPath() == "output/") {para->setOutputPath(outPath);}
    para->setOutputPrefix(simulationName);

    para->setPrintFiles(true);
    std::cout << "Write result files to " << para->getFName() << std::endl;

    para->setMainKernel(vf::CollisionKernel::Compressible::CumulantK17);
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vf::logging::Logger::changeLogPath(para->getOutputPath());
    VF_LOG_INFO("LB parameters:");
    VF_LOG_INFO("velocity LB [dx/dt]              = {}", vxLB);
    VF_LOG_INFO("viscosity LB [dx/dt]             = {}", viscosityLB);
    VF_LOG_INFO("dxGrid [-]                       = {}\n", dxGrid);
    VF_LOG_INFO("dt [s]                           = {}", dt);
    VF_LOG_INFO("simulation parameters:");
    VF_LOG_INFO("mainKernel                       = {}\n", para->getMainKernel());

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (useGridGenerator) {
        const real xGridMin = -0.5 * L;
        const real xGridMax = 0.5 * L;
        const real yGridMin = -0.5 * L;
        const real yGridMax = 0.5 * L;
        const real zGridMin = -0.5 * L;
        const real zGridMax = 0.5 * L;

        SPtr<Cuboid> level1 = nullptr;
        if (useLevels)
            level1 = std::make_shared<Cuboid>(-0.25 * L, -0.25 * L, -0.25 * L, 0.25 * L, 0.25 * L, 0.25 * L);

        if (para->getNumprocs() > 1) {

            const uint generatePart = vf::gpu::Communicator::getInstance().getPID();
            real overlap            = (real)8.0 * dxGrid;
            gridBuilder->setNumberOfLayers(10, 8);

            const real xSplit = 0.0;
            const real ySplit = 0.0;
            const real zSplit = 0.0;

            if (communicator.getNumberOfProcess() == 2) {

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, yGridMax, zSplit + overlap,
                                               dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zSplit - overlap, xGridMax, yGridMax, zGridMax,
                                               dxGrid);
                }

                if (useLevels) {
                    gridBuilder->addGrid(level1, 1);
                }

                if (generatePart == 0) {
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, yGridMin, yGridMax, zGridMin, zSplit));
                }
                if (generatePart == 1) {
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, yGridMin, yGridMax, zSplit, zGridMax));
                }

                gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

                if (generatePart == 0) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 1);
                }

                if (generatePart == 1) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 0);
                }

                gridBuilder->setPeriodicBoundaryCondition(false, false, false);
                //////////////////////////////////////////////////////////////////////////
                if (generatePart == 0)
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
                if (generatePart == 1)
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                //////////////////////////////////////////////////////////////////////////
            } else if (communicator.getNumberOfProcess() == 4) {

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xSplit + overlap, yGridMax,
                                               zSplit + overlap, dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, yGridMin, zGridMin, xGridMax, yGridMax,
                                               zSplit + overlap, dxGrid);
                }
                if (generatePart == 2) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zSplit - overlap, xSplit + overlap, yGridMax,
                                               zGridMax, dxGrid);
                }
                if (generatePart == 3) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, yGridMin, zSplit - overlap, xGridMax, yGridMax,
                                               zGridMax, dxGrid);
                }

                if (useLevels) {
                    gridBuilder->addGrid(level1, 1);
                }

                if (generatePart == 0)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xSplit, yGridMin, yGridMax, zGridMin, zSplit));
                if (generatePart == 1)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xSplit, xGridMax, yGridMin, yGridMax, zGridMin, zSplit));
                if (generatePart == 2)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xSplit, yGridMin, yGridMax, zSplit, zGridMax));
                if (generatePart == 3)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xSplit, xGridMax, yGridMin, yGridMax, zSplit, zGridMax));

                gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

                if (generatePart == 0) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 1);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 2);
                }
                if (generatePart == 1) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 0);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 3);
                }
                if (generatePart == 2) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 3);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 0);
                }
                if (generatePart == 3) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 2);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 1);
                }

                gridBuilder->setPeriodicBoundaryCondition(false, false, false);
                //////////////////////////////////////////////////////////////////////////
                if (generatePart == 0) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
                }
                if (generatePart == 2) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);
                }
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                if (generatePart == 3) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);
                }
                if (generatePart == 1) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                }
                //////////////////////////////////////////////////////////////////////////
            } else if (communicator.getNumberOfProcess() == 8) {

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xSplit + overlap, ySplit + overlap,
                                               zSplit + overlap, dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid(xGridMin, ySplit - overlap, zGridMin, xSplit + overlap, yGridMax,
                                               zSplit + overlap, dxGrid);
                }
                if (generatePart == 2) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, yGridMin, zGridMin, xGridMax, ySplit + overlap,
                                               zSplit + overlap, dxGrid);
                }
                if (generatePart == 3) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, ySplit - overlap, zGridMin, xGridMax, yGridMax,
                                               zSplit + overlap, dxGrid);
                }
                if (generatePart == 4) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zSplit - overlap, xSplit + overlap, ySplit + overlap,
                                               zGridMax, dxGrid);
                }
                if (generatePart == 5) {
                    gridBuilder->addCoarseGrid(xGridMin, ySplit - overlap, zSplit - overlap, xSplit + overlap, yGridMax,
                                               zGridMax, dxGrid);
                }
                if (generatePart == 6) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, yGridMin, zSplit - overlap, xGridMax, ySplit + overlap,
                                               zGridMax, dxGrid);
                }
                if (generatePart == 7) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, ySplit - overlap, zSplit - overlap, xGridMax, yGridMax,
                                               zGridMax, dxGrid);
                }

                if (useLevels) {
                    gridBuilder->addGrid(level1, 1);
                }

                if (generatePart == 0)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xSplit, yGridMin, ySplit, zGridMin, zSplit));
                if (generatePart == 1)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xSplit, ySplit, yGridMax, zGridMin, zSplit));
                if (generatePart == 2)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xSplit, xGridMax, yGridMin, ySplit, zGridMin, zSplit));
                if (generatePart == 3)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xSplit, xGridMax, ySplit, yGridMax, zGridMin, zSplit));
                if (generatePart == 4)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xSplit, yGridMin, ySplit, zSplit, zGridMax));
                if (generatePart == 5)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xSplit, ySplit, yGridMax, zSplit, zGridMax));
                if (generatePart == 6)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xSplit, xGridMax, yGridMin, ySplit, zSplit, zGridMax));
                if (generatePart == 7)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xSplit, xGridMax, ySplit, yGridMax, zSplit, zGridMax));

                gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!
                gridBuilder->setPeriodicBoundaryCondition(false, false, false);

                if (generatePart == 0) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 1);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 2);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 4);
                }
                if (generatePart == 1) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 0);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 3);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 5);
                }
                if (generatePart == 2) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 3);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 0);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 6);
                }
                if (generatePart == 3) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 2);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 1);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 7);
                }
                if (generatePart == 4) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 5);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 6);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 0);
                }
                if (generatePart == 5) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 4);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 7);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 1);
                }
                if (generatePart == 6) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 7);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 4);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 2);
                }
                if (generatePart == 7) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 6);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 5);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 3);
                }

                //////////////////////////////////////////////////////////////////////////
                if (generatePart == 0) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
                }
                if (generatePart == 1) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
                }
                if (generatePart == 2) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                }
                if (generatePart == 3) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                }
                if (generatePart == 4) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);
                }
                if (generatePart == 5) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);
                }
                if (generatePart == 6) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);
                }
                if (generatePart == 7) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);
                }
                //////////////////////////////////////////////////////////////////////////
            }

            // gridBuilder->writeGridsToVtk(outPath +  "/grid/part" + std::to_string(generatePart) + "_");
            // gridBuilder->writeArrows(outPath + "/" + std::to_string(generatePart) + " /arrow");

            SimulationFileWriter::write(gridPath + std::to_string(generatePart) + "/", gridBuilder, FILEFORMAT::BINARY);
        } else {

            gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, yGridMax, zGridMax, dxGrid);

            if (useLevels) {
                gridBuilder->setNumberOfLayers(10, 8);
                gridBuilder->addGrid(level1, 1);
            }

            gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!
            gridBuilder->setPeriodicBoundaryCondition(false, false, false);
            //////////////////////////////////////////////////////////////////////////
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);

            //////////////////////////////////////////////////////////////////////////
            gridBuilder->writeGridsToVtk(outPath + "/grid/");
            // gridBuilder->writeArrows(outPath + "/arrow");

            SimulationFileWriter::write(gridPath, gridBuilder, FILEFORMAT::BINARY);
        }
        bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    SPtr<GridProvider> gridGenerator;
    if (useGridGenerator)
        gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);
    else {
        gridGenerator = GridProvider::makeGridReader(FILEFORMAT::BINARY, para, cudaMemoryManager);
    }

    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, &scalingFactory);
    sim.run();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    std::string str, str2, configFile;

    if (argv != NULL) {

        try {
            //////////////////////////////////////////////////////////////////////////
            // assuming that a config files is stored parallel to this file.
            std::filesystem::path configPath = __FILE__;

            // the config file's default name can be replaced by passing a command line argument
            std::string configName("configDrivenCavityMultiGPU.txt");
            if (argc == 2) {
                configName = argv[1];
                std::cout << "Using configFile command line argument: " << configName << std::endl;
            }

            configPath.replace_filename(configName);

            multipleLevel(configPath);

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

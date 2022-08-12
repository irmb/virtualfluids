#define _USE_MATH_DEFINES
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "mpi.h"

//////////////////////////////////////////////////////////////////////////

#include "basics/Core/DataTypes.h"
#include "basics/Core/VectorTypes.h"
#include "basics/PointerDefinitions.h"

#include "basics/Core/LbmOrGks.h"
#include "basics/Core/Logger/Logger.h"
#include "basics/Core/StringUtilities/StringUtil.h"
#include "basics/config/ConfigurationFile.h"

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
#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

//////////////////////////////////////////////////////////////////////////

#include "utilities/communication.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          U s e r    s e t t i n g s
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//  Tesla 03
// const std::string outPath("E:/temp/DrivenCavityMultiGPUResults/");
// const std::string gridPath = "D:/STLs/DrivenCavity";
// const std::string simulationName("DrivenCavityMultiGPU");

// Phoenix
// const std::string outPath("/work/y0078217/Results/DrivenCavityMultiGPUResults/");
// const std::string gridPath = "/work/y0078217/Grids/GridDrivenCavityMultiGPU/";
// const std::string simulationName("DrivenCavityMultiGPU");

//  Aragorn
const std::string outPath("output/DrivenCavity_Results/");
const std::string gridPath = "output/DrivenCavity_Results/grid/";
const std::string simulationName("DrivenCavity");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string &configPath)
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
    config.load(configPath);
    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNummberOfProcess(), communicator.getPID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool useGridGenerator = true;
    bool useLevels        = true;
    // para->setUseStreams(useStreams);                  // set in config
    // para->useReducedCommunicationAfterFtoC = true;    // set in config
    para->setCalcTurbulenceIntensity(false);

    if (para->getNumprocs() == 1) {
        para->useReducedCommunicationAfterFtoC = false;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real L        = 1.0;
    const real Re       = 1000.0; // 1000
    const real velocity = 1.0;
    const real dt       = (real)1.0e-3; // 0.5e-3;
    const uint nx       = 64;
    std::string simulationName("DrivenCavityChimMultiGPU");

    // para->setTimestepOut(10000);   // set in config
    // para->setTimestepEnd(10000);   // set in config

    const real dxGrid      = L / real(nx);
    const real velocityLB  = velocity * dt / dxGrid;       // LB units
    const real vxLB        = velocityLB / (real)sqrt(2.0); // LB units
    const real vyLB        = velocityLB / (real)sqrt(2.0); // LB units
    const real viscosityLB = nx * velocityLB / Re;         // LB units

    *logging::out << logging::Logger::INFO_HIGH << "velocity  [dx/dt] = " << velocityLB << " \n";
    *logging::out << logging::Logger::INFO_HIGH << "viscosity [dx^2/dt] = " << viscosityLB << "\n";

    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(velocity / velocityLB);
    para->setDensityRatio((real)1.0); // correct value?

    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
        rho = (real)1.0;
        vx  = (real)(coordX * velocityLB);
        vy  = (real)(coordY * velocityLB);
        vz  = (real)(coordZ * velocityLB);
    });

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    para->setCalcDragLift(false);
    para->setUseWale(false);

    para->setOutputPrefix(simulationName);

    para->setPrintFiles(true);
    std::cout << "Write result files to " << para->getFName() << std::endl;

    if (useLevels)
        para->setMaxLevel(2);
    else
        para->setMaxLevel(1);

    // para->setMainKernel("CumulantK17CompChim");
    para->setMainKernel("CumulantK17CompChimStream");
    *logging::out << logging::Logger::INFO_HIGH << "Kernel: " << para->getMainKernel() << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (useGridGenerator) {
        const real xGridMin = -0.5 * L;
        const real xGridMax = 0.5 * L;
        const real yGridMin = -0.5 * L;
        const real yGridMax = 0.5 * L;
        const real zGridMin = -0.5 * L;
        const real zGridMax = 0.5 * L;

        Cuboid *level1 = nullptr;
        if (useLevels)
            level1 = new Cuboid(-0.25 * L, -0.25 * L, -0.25 * L, 0.25 * L, 0.25 * L, 0.25 * L);

        if (para->getNumprocs() > 1) {

            const uint generatePart = vf::gpu::Communicator::getInstance().getPID();
            real overlap            = (real)8.0 * dxGrid;
            gridBuilder->setNumberOfLayers(10, 8);

            const real xSplit = 0.0;
            const real ySplit = 0.0;
            const real zSplit = 0.0;

            if (communicator.getNummberOfProcess() == 2) {

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

                gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

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
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                //////////////////////////////////////////////////////////////////////////
            } else if (communicator.getNummberOfProcess() == 4) {

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

                gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

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
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                if (generatePart == 3) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 1) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                }
                //////////////////////////////////////////////////////////////////////////
            } else if (communicator.getNummberOfProcess() == 8) {

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

                gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!
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
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 5) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 6) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 7) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
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

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!
            gridBuilder->setPeriodicBoundaryCondition(false, false, false);
            //////////////////////////////////////////////////////////////////////////
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);

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

    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory);
    sim.run();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

int main(int argc, char *argv[])
{
    std::string str, str2, configFile;

    if (argv != NULL) {

        try {
            //////////////////////////////////////////////////////////////////////////

            std::string targetPath;

            targetPath = __FILE__;

            if (argc == 2) {
                configFile = argv[1];
                std::cout << "Using configFile command line argument: " << configFile << std::endl;
            }

#ifdef _WIN32
            targetPath = targetPath.substr(0, targetPath.find_last_of('\\') + 1);
#else
            targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);
#endif

            std::cout << targetPath << std::endl;

            if (configFile.size() == 0) {
                configFile = targetPath + "configDrivenCavityMultiGPU.txt";
            }

            multipleLevel(configFile);

            //////////////////////////////////////////////////////////////////////////
        } catch (const std::bad_alloc &e) {
            *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
        } catch (const std::exception &e) {
            *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
        } catch (...) {
            *logging::out << logging::Logger::LOGGER_ERROR << "Unknown exception!\n";
        }
    }
    return 0;
}
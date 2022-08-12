#define _USE_MATH_DEFINES
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <math.h>
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
#include "logger/Logger.h"

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "geometries/Sphere/Sphere.h"
#include "geometries/TriangularMesh/TriangularMesh.h"

#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"
#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"
#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"

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

// Tesla 03
// const std::string outPath("E:/temp/MusselOysterResults/");
// const std::string gridPathParent = "E:/temp/GridMussel/";
// const std::string stlPath("C:/Users/Master/Documents/MasterAnna/STL/");
// const std::string simulationName("MusselOyster");

// Aragorn
const std::string outPath("./output/MusselOysterResults/");
const std::string gridPathParent = "./output/MusselOysterResults/grid/";
const std::string stlPath("./stl/MusselOyster/");
const std::string simulationName("MusselOyster");

// Phoenix
// const std::string outPath("/work/y0078217/Results/MusselOysterResults/");
// const std::string gridPathParent = "/work/y0078217/Grids/GridMusselOyster/";
// const std::string stlPath("/home/y0078217/STL/MusselOyster/");
// const std::string simulationName("MusselOyster");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(std::filesystem::path &configPath)
{
    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    vf::gpu::Communicator &communicator = vf::gpu::Communicator::getInstance();

    auto gridFactory = GridFactory::make();
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    vf::basics::ConfigurationFile config;
    config.load(configPath.string());
    SPtr<Parameter> para =
        std::make_shared<Parameter>(communicator.getNummberOfProcess(), communicator.getPID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool useGridGenerator                  = true;
    bool useStreams                        = true;
    bool useLevels                         = true;
    para->useReducedCommunicationAfterFtoC = true;
    para->setCalcTurbulenceIntensity(true);

    if (para->getNumprocs() == 1) {
        para->useReducedCommunicationAfterFtoC = false;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::string bivalveType = "MUSSEL"; // "MUSSEL" "OYSTER"
    std::string gridPath(
        gridPathParent +
        bivalveType); // only for GridGenerator, for GridReader the gridPath needs to be set in the config file

    real dxGrid = (real)2.0; // 2.0
    // real dxGrid = (real)1.0; // 1.0
    if (para->getNumprocs() == 8)
        dxGrid = 0.5;
    real vxLB = (real)0.051; // LB units
    real Re = (real)300.0;
    real referenceLength = 1.0 / dxGrid; // heightBivalve / dxGrid
    real viscosityLB = (vxLB * referenceLength) / Re;

    para->setVelocityLB(vxLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio((real)58.82352941);
    para->setViscosityRatio((real)0.058823529);
    para->setDensityRatio((real)998.0);

    VF_LOG_INFO("LB parameters:");
    VF_LOG_INFO("velocity LB [dx/dt]              = {}", vxLB);
    VF_LOG_INFO("viscosity LB [dx/dt]             = {}", viscosityLB);
    VF_LOG_INFO("dxGrid [-]                       = {}\n", dxGrid);

    VF_LOG_INFO("world parameters:");
    VF_LOG_INFO("velocity [m/s]                   = {}", vxLB * para->getVelocityRatio());
    VF_LOG_INFO("viscosity [m^2/s]                = {}\n", viscosityLB * para->getViscosityRatio());

    VF_LOG_INFO("simulation parameters:");
    VF_LOG_INFO("useGridGenerator                 = {}", useGridGenerator);
    VF_LOG_INFO("useStreams                       = {}", para->getUseStreams());
    VF_LOG_INFO("number of processes              = {}", para->getNumprocs());
    VF_LOG_INFO("useReducedCommunicationAfterFtoC = {}", para->useReducedCommunicationAfterFtoC);
    VF_LOG_INFO("bivalveType                      = {}", bivalveType);
    VF_LOG_INFO("mainKernel                       = {}\n", para->getMainKernel());

    // para->setTimestepOut(1000);
    // para->setTimestepEnd(10000);

    para->setCalcDragLift(false);
    para->setUseWale(false);

    para->setOutputPrefix(simulationName);

    para->setPrintFiles(true);
    std::cout << "Write result files to " << para->getFName() << std::endl;

    para->setUseStreams(useStreams);
    // para->setMainKernel("CumulantK17CompChim");
    para->setMainKernel("CumulantK17CompChimStream");
    *logging::out << logging::Logger::INFO_HIGH << "Kernel: " << para->getMainKernel() << "\n";

    //////////////////////////////////////////////////////////////////////////

    if (useGridGenerator) {
        const real xGridMin = -100.0; // -100.0;
        const real xGridMax = 470.0;  // alt 540.0 // neu 440 // mit groesserem Level 1 470
        const real yGridMin = 1.0;    // 1.0;
        const real yGridMax = 350.0;  // alt 440.0; // neu 350
        const real zGridMin = -85;    // -85;
        const real zGridMax = 85.0;   // 85;

        // height MUSSEL = 35.0
        // height Oyster = 72.0

        TriangularMesh *bivalveSTL = TriangularMesh::make(stlPath + bivalveType + ".stl");
        TriangularMesh *bivalveRef_1_STL = nullptr;
        if (useLevels)
            bivalveRef_1_STL = TriangularMesh::make(stlPath + bivalveType + "_Level1.stl");

        if (para->getNumprocs() > 1) {
            const uint generatePart = vf::gpu::Communicator::getInstance().getPID();

            real overlap = (real)8.0 * dxGrid;
            gridBuilder->setNumberOfLayers(10, 8);

            if (communicator.getNummberOfProcess() == 2) {
                const real zSplit = 0.0; // round(((double)bbzp + bbzm) * 0.5);

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, yGridMax, zSplit + overlap,
                                               dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zSplit - overlap, xGridMax, yGridMax, zGridMax,
                                               dxGrid);
                }

                if (useLevels) {
                    gridBuilder->addGrid(bivalveRef_1_STL, 1);
                }

                gridBuilder->addGeometry(bivalveSTL);

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
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                if (generatePart == 1)
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
                gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs
                //////////////////////////////////////////////////////////////////////////
            } else if (communicator.getNummberOfProcess() == 4) {

                const real xSplit = 100.0;
                const real zSplit = 0.0;

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
                    gridBuilder->addGrid(bivalveRef_1_STL, 1);
                }

                gridBuilder->addGeometry(bivalveSTL);

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
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                }
                if (generatePart == 2) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                if (generatePart == 3) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs
                }
                if (generatePart == 1) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs
                }
                //////////////////////////////////////////////////////////////////////////
            } else if (communicator.getNummberOfProcess() == 8) {
                real xSplit = 140.0; // 100.0 // mit groesserem Level 1 140.0
                real ySplit = 32.0;  // 32.0
                real zSplit = 0.0;

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
                    gridBuilder->addGrid(bivalveRef_1_STL, 1);
                }

                gridBuilder->addGeometry(bivalveSTL);

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
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 1) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 2) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs
                }
                if (generatePart == 3) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs
                }
                if (generatePart == 4) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 5) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 6) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs
                }
                if (generatePart == 7) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs
                }
                // gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
                //////////////////////////////////////////////////////////////////////////
            }

            // gridBuilder->writeGridsToVtk(outPath +  bivalveType + "/grid/part" + std::to_string(generatePart) + "_");
            // gridBuilder->writeArrows(outPath + bivalveType + "/" + std::to_string(generatePart) + " /arrow");
            // SimulationFileWriter::write(gridPath + std::to_string(generatePart) + "/", gridBuilder,
            //                             FILEFORMAT::BINARY);
        } else {

            gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, yGridMax, zGridMax, dxGrid);

            if (useLevels) {
                gridBuilder->setNumberOfLayers(10, 8);
                gridBuilder->addGrid(bivalveRef_1_STL, 1);
            }

            gridBuilder->addGeometry(bivalveSTL);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);
            //////////////////////////////////////////////////////////////////////////
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs

            //////////////////////////////////////////////////////////////////////////

            // gridBuilder->writeGridsToVtk(outPath +  bivalveType + "/grid/");
            // gridBuilder->writeArrows ((outPath + bivalveType + "/arrow");

            SimulationFileWriter::write(gridPath, gridBuilder, FILEFORMAT::BINARY);
        }

        bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
        bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);
        bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
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
    MPI_Init(&argc, &argv);
    std::string str, str2, configFile;

    if (argv != NULL) {

        try {
            //////////////////////////////////////////////////////////////////////////
            // assuming that a config files is stored parallel to this file.
            std::filesystem::path configPath = __FILE__;

            // the config file's default name can be replaced by passing a command line argument
            std::string configName("config.txt");
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

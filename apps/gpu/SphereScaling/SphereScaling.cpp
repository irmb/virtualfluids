#define _USE_MATH_DEFINES
#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <filesystem>

#include "mpi.h"

//////////////////////////////////////////////////////////////////////////

#include "basics/DataTypes.h"
#include "basics/PointerDefinitions.h"

#include "basics/StringUtilities/StringUtil.h"
#include "basics/config/ConfigurationFile.h"
#include <logger/Logger.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

#include "geometries/Conglomerate/Conglomerate.h"
#include "geometries/Cuboid/Cuboid.h"
#include "geometries/Sphere/Sphere.h"
#include "geometries/TriangularMesh/TriangularMesh.h"

#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"
#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/Communication/MpiCommunicator.h"
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

void runVirtualFluids(const vf::basics::ConfigurationFile& config)
{
    vf::gpu::Communicator& communicator = vf::gpu::MpiCommunicator::getInstance();

    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNumberOfProcess(), communicator.getPID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
    GridScalingFactory scalingFactory = GridScalingFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool useGridGenerator   = true;
    bool useLevels = true;
    std::string scalingType = "strong"; // "strong" // "weak"

    const std::string outPath("output/" + std::to_string(para->getNumprocs()) + "GPU/");
    const std::string simulationName("SphereScaling");
    const std::string gridPath = "./output/grids/";
    const std::string stlPath("./stl/SphereScaling/");

    if (para->getNumprocs() == 1) {
        para->useReducedCommunicationAfterFtoC = false;
    }
    if (scalingType != "weak" && scalingType != "strong")
        std::cerr << "unknown scaling type" << std::endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real dxGrid      = (real)1.0;
    real vxLB        = (real)0.005;  // LB units
    real viscosityLB = 0.001;        //(vxLB * dxGrid) / Re;

    para->setVelocityLB(vxLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio((real)58.82352941);
    para->setViscosityRatio((real)0.058823529);
    para->setDensityRatio((real)998.0);

    para->setCalcDragLift(false);
    para->setUseWale(false);

    para->setOutputPrefix(simulationName);
    if (para->getOutputPath() == "output/") {para->setOutputPath(outPath);}
    para->setPrintFiles(true);

    if (useLevels)
        para->setMaxLevel(2);
    else
        para->setMaxLevel(1);

    para->setMainKernel(vf::CollisionKernel::Compressible::CumulantK17);
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    VF_LOG_INFO("Number of processes: {}", para->getNumprocs());

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
    VF_LOG_INFO("scalingType                      = {}", scalingType);
    VF_LOG_INFO("mainKernel                       = {}\n", para->getMainKernel());

    //////////////////////////////////////////////////////////////////////////
    auto gridBuilder = std::make_shared<MultipleGridBuilder>();

    if (useGridGenerator) {
        real sideLengthCube;
        if (useLevels) {
            if (scalingType == "strong")
                sideLengthCube = 76.0; // Phoenix: strong scaling with two levels = 76.0
            else if (scalingType == "weak")
                sideLengthCube = 70.0; // Phoenix: weak scaling with two levels = 70.0
        } else
            sideLengthCube = 92.0; // Phoenix: 86.0
        real xGridMin          = 0.0;
        real xGridMax          = sideLengthCube;
        real yGridMin          = 0.0;
        real yGridMax          = sideLengthCube;
        real zGridMin          = 0.0;
        real zGridMax          = sideLengthCube;
        const real dSphere     = 10.0;
        const real dSphereLev1 = 22.0; // Phoenix: 22.0
        const real dCubeLev1   = 72.0; // Phoenix: 72.0

        if (para->getNumprocs() > 1) {
            const uint generatePart = vf::gpu::MpiCommunicator::getInstance().getPID();

            real overlap = (real)8.0 * dxGrid;
            gridBuilder->setNumberOfLayers(10, 8);

            if (communicator.getNumberOfProcess() == 2) {
                real zSplit = 0.5 * sideLengthCube;

                if (scalingType == "weak") {
                    zSplit   = zGridMax;
                    zGridMax = zGridMax + sideLengthCube;
                }

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, yGridMax, zSplit + overlap,
                                               dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zSplit - overlap, xGridMax, yGridMax, zGridMax,
                                               dxGrid);
                }

                if (useLevels) {
                    if (scalingType == "strong") {
                        gridBuilder->addGrid(
                            std::make_shared<Sphere>(0.5 * sideLengthCube, 0.5 * sideLengthCube, 0.5 * sideLengthCube, dSphereLev1),
                            1);
                    } else if (scalingType == "weak") {
                        gridBuilder->addGrid(std::make_shared<Cuboid>(-0.5 * dCubeLev1, -0.5 * dCubeLev1,
                                                        sideLengthCube - 0.5 * dCubeLev1, 0.5 * dCubeLev1,
                                                        0.5 * dCubeLev1, sideLengthCube + 0.5 * dCubeLev1),
                                             1);
                    }
                }

                if (scalingType == "weak") {
                    if (useLevels) {
                        gridBuilder->addGeometry(std::make_shared<Sphere>(0.0, 0.0, sideLengthCube, dSphere));
                    } else {
                        auto sphereSTL = std::make_shared<TriangularMesh>(stlPath + "Spheres_2GPU.stl");
                        gridBuilder->addGeometry(sphereSTL);
                    }
                } else if (scalingType == "strong") {
                    gridBuilder->addGeometry(
                        std::make_shared<Sphere>(0.5 * sideLengthCube, 0.5 * sideLengthCube, 0.5 * sideLengthCube, dSphere));
                }

                if (generatePart == 0)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, yGridMin, yGridMax, zGridMin, zSplit));
                if (generatePart == 1)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, yGridMin, yGridMax, zSplit, zGridMax));

                gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

                if (generatePart == 0) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 1);
                }

                if (generatePart == 1) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 0);
                }

                gridBuilder->setPeriodicBoundaryCondition(false, false, false);
                //////////////////////////////////////////////////////////////////////////
                gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                if (generatePart == 0)
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                if (generatePart == 1)
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs
                // gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
                //////////////////////////////////////////////////////////////////////////

            } else if (communicator.getNumberOfProcess() == 4) {
                real ySplit = 0.5 * sideLengthCube;
                real zSplit = 0.5 * sideLengthCube;

                if (scalingType == "weak") {
                    ySplit   = yGridMax;
                    yGridMax = yGridMax + (yGridMax - yGridMin);
                    zSplit   = zGridMax;
                    zGridMax = zGridMax + (zGridMax - zGridMin);
                }

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, ySplit + overlap,
                                               zSplit + overlap, dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid(xGridMin, ySplit - overlap, zGridMin, xGridMax, yGridMax,
                                               zSplit + overlap, dxGrid);
                }
                if (generatePart == 2) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zSplit - overlap, xGridMax, ySplit + overlap,
                                               zGridMax, dxGrid);
                }
                if (generatePart == 3) {
                    gridBuilder->addCoarseGrid(xGridMin, ySplit - overlap, zSplit - overlap, xGridMax, yGridMax,
                                               zGridMax, dxGrid);
                }

                if (useLevels) {
                    if (scalingType == "strong") {
                        gridBuilder->addGrid(
                            std::make_shared<Sphere>(0.5 * sideLengthCube, 0.5 * sideLengthCube, 0.5 * sideLengthCube, dSphereLev1),
                            1);
                    } else if (scalingType == "weak") {
                        gridBuilder->addGrid(std::make_shared<Cuboid>(-0.5 * dCubeLev1, sideLengthCube - 0.5 * dCubeLev1,
                                                        sideLengthCube - 0.5 * dCubeLev1, 0.5 * dCubeLev1,
                                                        sideLengthCube + 0.5 * dCubeLev1,
                                                        sideLengthCube + 0.5 * dCubeLev1),
                                             1);
                    }
                }

                if (scalingType == "weak") {
                    if (useLevels) {
                        gridBuilder->addGeometry(std::make_shared<Sphere>(0.0, sideLengthCube, sideLengthCube, dSphere));
                    } else {
                        auto sphereSTL = std::make_shared<TriangularMesh>(stlPath + "Spheres_4GPU.stl");
                        gridBuilder->addGeometry(sphereSTL);
                    }
                } else if (scalingType == "strong") {
                    gridBuilder->addGeometry(
                        std::make_shared<Sphere>(0.5 * sideLengthCube, 0.5 * sideLengthCube, 0.5 * sideLengthCube, dSphere));
                }

                if (generatePart == 0)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, yGridMin, ySplit, zGridMin, zSplit));
                if (generatePart == 1)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, ySplit, yGridMax, zGridMin, zSplit));
                if (generatePart == 2)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, yGridMin, ySplit, zSplit, zGridMax));
                if (generatePart == 3)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, ySplit, yGridMax, zSplit, zGridMax));

                gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!
                gridBuilder->setPeriodicBoundaryCondition(false, false, false);

                if (generatePart == 0) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 1);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 2);
                }
                if (generatePart == 1) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 0);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 3);
                }
                if (generatePart == 2) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 3);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 0);
                }
                if (generatePart == 3) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 2);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 1);
                }

                //////////////////////////////////////////////////////////////////////////
                if (generatePart == 0) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 1) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 2) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 3) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs
                // gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
                //////////////////////////////////////////////////////////////////////////
            } else if (communicator.getNumberOfProcess() == 8) {
                real xSplit = 0.5 * sideLengthCube;
                real ySplit = 0.5 * sideLengthCube;
                real zSplit = 0.5 * sideLengthCube;

                if (scalingType == "weak") {
                    xSplit   = xGridMax;
                    xGridMax = xGridMax + (xGridMax - xGridMin);
                    ySplit   = yGridMax;
                    yGridMax = yGridMax + (yGridMax - yGridMin);
                    zSplit   = zGridMax;
                    zGridMax = zGridMax + (zGridMax - zGridMin);
                }

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
                    if (scalingType == "strong") {
                        gridBuilder->addGrid(
                            std::make_shared<Sphere>(0.5 * sideLengthCube, 0.5 * sideLengthCube, 0.5 * sideLengthCube, dSphereLev1),
                            1);
                    } else if (scalingType == "weak") {
                        gridBuilder->addGrid(
                            std::make_shared<Cuboid>(sideLengthCube - 0.5 * dCubeLev1, sideLengthCube - 0.5 * dCubeLev1,
                                       sideLengthCube - 0.5 * dCubeLev1, sideLengthCube + 0.5 * dCubeLev1,
                                       sideLengthCube + 0.5 * dCubeLev1, sideLengthCube + 0.5 * dCubeLev1),
                            1);
                    }
                }

                if (scalingType == "weak") {
                    if (useLevels) {
                        gridBuilder->addGeometry(std::make_shared<Sphere>(sideLengthCube, sideLengthCube, sideLengthCube, dSphere));
                    } else {
                        auto sphereSTL = std::make_shared<TriangularMesh>(stlPath + "Spheres_8GPU.stl");
                        gridBuilder->addGeometry(sphereSTL);
                    }
                } else if (scalingType == "strong") {
                    gridBuilder->addGeometry(
                        std::make_shared<Sphere>(0.5 * sideLengthCube, 0.5 * sideLengthCube, 0.5 * sideLengthCube, dSphere));
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
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 1);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 2);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 4);
                }
                if (generatePart == 1) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 0);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 3);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 5);
                }
                if (generatePart == 2) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 3);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 0);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 6);
                }
                if (generatePart == 3) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 2);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 1);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, 7);
                }
                if (generatePart == 4) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 5);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 6);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 0);
                }
                if (generatePart == 5) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 4);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 7);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 1);
                }
                if (generatePart == 6) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 7);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 4);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 2);
                }
                if (generatePart == 7) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 6);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 5);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MZ);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, 3);
                }

                //////////////////////////////////////////////////////////////////////////
                if (generatePart == 0) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 1) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 2) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
                }
                if (generatePart == 3) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
                }
                if (generatePart == 4) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 5) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 6) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
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

            // gridBuilder->writeGridsToVtk(outPath + "grid/part" + std::to_string(generatePart) + "_");
            // gridBuilder->writeGridsToVtk(outPath +std::to_string(generatePart) + "/grid/");
            // gridBuilder->writeArrows(outPath + std::to_string(generatePart) + " /arrow");

            SimulationFileWriter::write(gridPath + std::to_string(generatePart) + "/", gridBuilder, FILEFORMAT::BINARY);
        } else {

            gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, yGridMax, zGridMax, dxGrid);

            if (useLevels) {
                gridBuilder->setNumberOfLayers(10, 8);
                if (scalingType == "strong") {
                    gridBuilder->addGrid(
                        std::make_shared<Sphere>(0.5 * sideLengthCube, 0.5 * sideLengthCube, 0.5 * sideLengthCube, dSphereLev1), 1);
                } else if (scalingType == "weak")
                    gridBuilder->addGrid(std::make_shared<Cuboid>(sideLengthCube - 0.5 * dCubeLev1, sideLengthCube - 0.5 * dCubeLev1,
                                                    sideLengthCube - 0.5 * dCubeLev1, sideLengthCube + 0.5 * dCubeLev1,
                                                    sideLengthCube + 0.5 * dCubeLev1, sideLengthCube + 0.5 * dCubeLev1),
                                         1);
            }

            if (scalingType == "weak") {
                if (useLevels) {
                    gridBuilder->addGeometry(std::make_shared<Sphere>(sideLengthCube, sideLengthCube, sideLengthCube, dSphere));
                } else {
                    auto sphereSTL = std::make_shared<TriangularMesh>(stlPath + "Spheres_1GPU.stl");
                    gridBuilder->addGeometry(sphereSTL);
                }
            } else {
                gridBuilder->addGeometry(
                    std::make_shared<Sphere>(0.5 * sideLengthCube, 0.5 * sideLengthCube, 0.5 * sideLengthCube, dSphere));
            }

            gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);
            //////////////////////////////////////////////////////////////////////////
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure BC after velocity BCs

            // gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
            //////////////////////////////////////////////////////////////////////////

            // gridBuilder->writeGridsToVtk("E:/temp/MusselOyster/" + "/grid/");
            // gridBuilder->writeArrows ("E:/temp/MusselOyster/" + "/arrow");

            SimulationFileWriter::write(gridPath, gridBuilder, FILEFORMAT::BINARY);
        }

        bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
        bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

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
            VF_LOG_INFO("For the default config path to work, execute the app from the project root.");
            vf::basics::ConfigurationFile config = vf::basics::loadConfig(argc, argv, "./apps/gpu/SphereScaling/config.txt");
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


#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <exception>
#include <memory>

#include "mpi.h"

//////////////////////////////////////////////////////////////////////////

#include "basics/Core/DataTypes.h"
#include "basics/PointerDefinitions.h"
#include "basics/Core/VectorTypes.h"

#include "basics/Core/LbmOrGks.h"
#include "basics/Core/StringUtilities/StringUtil.h"
#include "basics/config/ConfigurationFile.h"
#include "basics/Core/Logger/Logger.h"

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridFactory.h"

#include "geometries/Sphere/Sphere.h"
#include "geometries/TriangularMesh/TriangularMesh.h"

#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"

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
// std::string outPath("E:/temp/MusselOysterResults/");
// std::string gridPathParent = "E:/temp/GridMussel/";
// std::string stlPath("C:/Users/Master/Documents/MasterAnna/STL/");
// std::string simulationName("MusselOyster");

// Phoenix
std::string outPath("/work/y0078217/Results/MusselOysterResults/");
std::string gridPathParent = "/work/y0078217/Grids/GridMusselOyster/";
std::string stlPath("/home/y0078217/STL/MusselOyster/");
std::string simulationName("MusselOyster");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath)
{
    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    
	vf::gpu::Communicator* comm = vf::gpu::Communicator::getInstanz();
    vf::basics::ConfigurationFile config;
    std::cout << configPath << std::endl;
    config.load(configPath);
    SPtr<Parameter> para = std::make_shared<Parameter>(config, comm->getNummberOfProcess(), comm->getPID());



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool useGridGenerator = true;
    bool useStreams       = false;
    bool useLevels        = false;
    para->useReducedCommunicationAfterFtoC = false;

    if (para->getNumprocs() == 1) {
       useStreams       = false;
       para->useReducedCommunicationAfterFtoC = false;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::string bivalveType = "OYSTER"; // "MUSSEL" "OYSTER"
    std::string gridPath(gridPathParent + bivalveType); // only for GridGenerator, for GridReader the gridPath needs to be set in the config file

    real dxGrid = (real)0.1; // 0.1
    if (para->getNumprocs() == 4)
        dxGrid =0.66666667;
    else if (para->getNumprocs() == 8)  
        dxGrid =0.44444444;  
    real vxLB = (real)0.051; // LB units
    real Re = (real)300.0;

    real heightBivalve;
    if (bivalveType == "MUSSEL")
        heightBivalve = (real)35.0; 
    else if (bivalveType == "OYSTER")
        heightBivalve = (real)72.0;
    else
        std::cerr << "Error: unknown bivalveType" << std::endl;
    real length = 1.0 / dxGrid; // heightBivalve / dxGrid
    real viscosityLB = (vxLB * length) / Re;

    para->setVelocity(vxLB);
    para->setViscosity(viscosityLB);
    para->setVelocityRatio((real) 58.82352941);
    para->setViscosityRatio((real) 0.058823529);
    para->setDensityRatio((real) 998.0);

    *logging::out << logging::Logger::INFO_HIGH << "velocity LB [dx/dt] = " << vxLB << " \n";
    *logging::out << logging::Logger::INFO_HIGH << "viscosity LB [dx^2/dt] = " << viscosityLB << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "velocity real [m/s] = " << vxLB * para->getVelocityRatio()<< " \n";
    *logging::out << logging::Logger::INFO_HIGH << "viscosity real [m^2/s] = " << viscosityLB * para->getViscosityRatio() << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "dxGrid = " << dxGrid << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "useGridGenerator = " << useGridGenerator << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "useStreams = " << useStreams << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "number of processes = " << para->getNumprocs() << "\n";

    
    // para->setTOut(1000);
    // para->setTEnd(10000);

    para->setCalcDragLift(false);
    para->setUseWale(false);

    if (para->getOutputPath().size() == 0) {
        para->setOutputPath(outPath);
    }
    para->setOutputPrefix(simulationName);
    para->setFName(para->getOutputPath() + para->getOutputPrefix());
    para->setPrintFiles(true);
    std::cout << "Write result files to " << para->getFName() << std::endl;


    if (useLevels)
        para->setMaxLevel(2);
    else
        para->setMaxLevel(1);


    if (useStreams)
        para->setUseStreams();
    //para->setMainKernel("CumulantK17CompChim");
    para->setMainKernel("CumulantK17CompChimStream");
    *logging::out << logging::Logger::INFO_HIGH << "Kernel: " << para->getMainKernel() << "\n";

    // if (para->getNumprocs() > 1) {
    //     para->setDevices(std::vector<uint>{ (uint)0, (uint)1 });
    //     para->setMaxDev(2);
    // } else 
    //     para->setDevices(std::vector<uint>{ (uint)0 });



    //////////////////////////////////////////////////////////////////////////


    if (useGridGenerator) {
        const real bbzm = 0.0;
        real bbzp;
        if (bivalveType == "MUSSEL")
            bbzp = 18.3;
        if (bivalveType == "OYSTER")
            bbzp = 26.0;
        // bounding box mussel:
        // const real bbxm = 0.0;
        // const real bbxp = 76.0;
        // const real bbym = 0.0;
        // const real bbyp = 35.0;
        // const real bbzm = 0.0;
        // const real bbzp = 18.3;
        // bounding box oyster:
        // const real bbxm = 0.0;
        // const real bbxp = 102.0;
        // const real bbym = 0.0;
        // const real bbyp = 72.0;
        // const real bbzm = 0.0;
        // const real bbzp = 26.0;

        const real xGridMin  = -100;       // -100.0;
        const real xGridMax  = 540.0;      // 540.0
        const real yGridMin  = 1.0;        // 1.0;
        const real yGridMax  = 440;        // 440.0;
        const real zGridMin  = - 70.0;     // -70;
        const real zGridMax  = 100.0;      // 100;

        TriangularMesh *bivalveSTL       = TriangularMesh::make(stlPath + bivalveType + ".stl");
        TriangularMesh *bivalveRef_1_STL = nullptr;
        if (useLevels)
            bivalveRef_1_STL = TriangularMesh::make(stlPath + bivalveType + "_Level1.stl");


        if (para->getNumprocs() > 1) {
            const uint generatePart = vf::gpu::Communicator::getInstanz()->getPID();

            real overlap = (real)8.0 * dxGrid;
            gridBuilder->setNumberOfLayers(10, 8);

            if (comm->getNummberOfProcess() == 2) {
                const real zSplit = round(((double)bbzp + bbzm) * 0.5);          

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid( xGridMin,   yGridMin,     zGridMin, 
                                                xGridMax,   yGridMax,     zSplit+overlap,   dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid( xGridMin,    yGridMin,     zSplit-overlap, 
                                                xGridMax,    yGridMax,     zGridMax,        dxGrid);
                }


                if (useLevels) {
                    gridBuilder->addGrid(bivalveRef_1_STL, 1);
                }

                gridBuilder->addGeometry(bivalveSTL);

                if (generatePart == 0)
                    gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(xGridMin,    xGridMax,
                                                                               yGridMin,    yGridMax, 
                                                                               zGridMin,    zSplit));
                if (generatePart == 1)
                    gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(xGridMin,    xGridMax, 
                                                                               yGridMin,    yGridMax, 
                                                                               zSplit,      zGridMax));            

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
                gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
                //////////////////////////////////////////////////////////////////////////
                if (para->getKernelNeedsFluidNodeIndicesToRun())
                    gridBuilder->findFluidNodes(useStreams);

                //gridBuilder->writeGridsToVtk(outPath + "/" + bivalveType + "/grid/part" + std::to_string(generatePart) + "_");
                //gridBuilder->writeGridsToVtk(outPath + "/" + bivalveType + "/" + std::to_string(generatePart) + "/grid/");
                //gridBuilder->writeArrows(outPath + "/" + bivalveType + "/" + std::to_string(generatePart) + " /arrow");

                SimulationFileWriter::write(gridPath + "/" + std::to_string(generatePart) + "/", gridBuilder, FILEFORMAT::BINARY);
           
            } else if (comm->getNummberOfProcess() == 4) {

                const real xSplit = 100.0;
                const real zSplit = round(((double)bbzp + bbzm) * 0.5);

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xSplit + overlap, yGridMax,
                                               zSplit + overlap, dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, yGridMin, zGridMin, xGridMax, yGridMax,
                                               zSplit+overlap, dxGrid);
                }
                if (generatePart == 2) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zSplit-overlap, xSplit + overlap, yGridMax,
                                               zGridMax, dxGrid);
                }
                if (generatePart == 3) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, yGridMin, zSplit-overlap, xGridMax, yGridMax,
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
                if (generatePart == 1) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);                    
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
                }
                if (generatePart == 2) {                    
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                if (generatePart == 3) {
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);
                }
                gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                //////////////////////////////////////////////////////////////////////////
            }
            if (para->getKernelNeedsFluidNodeIndicesToRun())
                gridBuilder->findFluidNodes(useStreams);

            gridBuilder->writeGridsToVtk(outPath +  bivalveType + "/grid/part" + std::to_string(generatePart) + "_"); 
            // gridBuilder->writeGridsToVtk(outPath + bivalveType + "/" + std::to_string(generatePart) + "/grid/"); 
            // gridBuilder->writeArrows(outPath + bivalveType + "/" + std::to_string(generatePart) + " /arrow");

            SimulationFileWriter::write(gridPath + std::to_string(generatePart) + "/", gridBuilder,
                                        FILEFORMAT::BINARY);
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
            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, 0.0, 0.0);

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
            //////////////////////////////////////////////////////////////////////////
            if (para->getKernelNeedsFluidNodeIndicesToRun())
                gridBuilder->findFluidNodes(useStreams);

            gridBuilder->writeGridsToVtk(outPath +  bivalveType + "/grid/");
            // gridBuilder->writeArrows ((outPath + bivalveType + "/arrow");

            SimulationFileWriter::write(gridPath, gridBuilder, FILEFORMAT::BINARY);
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // const real velocityLB = velocity * dt / dx; // LB units

        // const real vx = velocityLB / (real)sqrt(2.0); // LB units
        // const real vy = velocityLB / (real)sqrt(2.0); // LB units

        // const real viscosityLB = nx * velocityLB / Re; // LB units

        //*logging::out << logging::Logger::INFO_HIGH << "velocity  [dx/dt] = " << velocityLB << " \n";
        //*logging::out << logging::Logger::INFO_HIGH << "viscosity [dx^2/dt] = " << viscosityLB << "\n";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // para->setVelocity(velocityLB);
        // para->setViscosity(viscosityLB);

        // para->setVelocityRatio(velocity/ velocityLB);

        // para->setMainKernel("CumulantK17CompChim");

        // para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz)
        // {
        //          rho = (real)0.0;
        //          vx  = (real)0.0; //(6 * velocityLB * coordZ * (L - coordZ) / (L * L));
        //          vy  = (real)0.0;
        //          vz  = (real)0.0;
        //      });



       //return;
    }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<CudaMemoryManager> cudaMemoryManager = CudaMemoryManager::make(para);

    SPtr<GridProvider> gridGenerator;
    if (useGridGenerator)
        gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager);
    else {
        gridGenerator = GridProvider::makeGridReader(FILEFORMAT::BINARY, para, cudaMemoryManager);
    }
           
    Simulation sim;
    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    SPtr<KernelFactoryImp> kernelFactory = KernelFactoryImp::getInstance();
    SPtr<PreProcessorFactoryImp> preProcessorFactory = PreProcessorFactoryImp::getInstance();
    sim.setFactories(kernelFactory, preProcessorFactory);
    sim.init(para, gridGenerator, fileWriter, cudaMemoryManager);
    sim.run();
    sim.free();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
}

int main( int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    std::string str, str2, configFile;

    if ( argv != NULL )
    {
        //str = static_cast<std::string>(argv[0]);
        
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
                configFile = targetPath + "configMusselOyster.txt";
            }

            multipleLevel(configFile);

            //////////////////////////////////////////////////////////////////////////
        }
        catch (const std::bad_alloc& e)
        { 
            *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
        }
        catch (const std::exception& e)
        {   
            *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
        }
        catch (...)
        {
            *logging::out << logging::Logger::LOGGER_ERROR << "Unknown exception!\n";
        }
    }

   MPI_Finalize();
   return 0;
}

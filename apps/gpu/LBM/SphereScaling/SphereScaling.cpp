
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
 std::string outPath("E:/temp/SphereScalingResults");
 std::string gridPathParent = "E:/temp/GridSphereScaling/";
 std::string simulationName("SphereScaling");

//// Phoenix
//std::string outPath("/work/y0078217/Results/SphereScalingResults/");
//std::string gridPathParent = "/work/y0078217/Grids/GridSphereScaling/";
//std::string simulationName("SphereScaling");

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

    bool useGridGenerator                  = true;
    bool useStreams                        = false;
    bool useLevels                         = true;
    para->useReducedCommunicationAfterFtoC = true;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::string gridPath(gridPathParent); // only for GridGenerator, for GridReader the gridPath needs to be set in the config file

    real dxGrid      = (real)0.2;
    real vxLB = (real)0.051; // LB units
    real Re = (real)300.0;
    real viscosityLB = (vxLB * dxGrid) / Re;

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

    
    para->setTOut(10);
    para->setTEnd(10);

    para->setCalcDragLift(false);
    para->setUseWale(false);

    if (para->getOutputPath().size() == 0) {
        para->setOutputPath(outPath);
    }
    para->setOutputPrefix(simulationName);
    para->setFName(para->getOutputPath() + "/" + para->getOutputPrefix());
    para->setPrintFiles(true);

    if (useLevels)
        para->setMaxLevel(2);
    else
        para->setMaxLevel(1);


    if (useStreams)
        para->setUseStreams();
    //para->setMainKernel("CumulantK17CompChim");
    para->setMainKernel("CumulantK17CompChimStream");
    *logging::out << logging::Logger::INFO_HIGH << "Kernel: " << para->getMainKernel() << "\n";

    if (para->getNumprocs() == 4) {
        para->setDevices(std::vector<uint>{ 0u, 1u, 2u, 3u });
        para->setMaxDev(4);
    } else if (para->getNumprocs() == 2) {
        para->setDevices(std::vector<uint>{ 0u, 1u });
        para->setMaxDev(2);
    } else 
        para->setDevices(std::vector<uint>{ 0u });
        para->setMaxDev(1);



    //////////////////////////////////////////////////////////////////////////


    if (useGridGenerator) {

        const real xGridMin    = -32;
        const real xGridMax    = 32;
        const real yGridMin    = xGridMin;
        const real yGridMax    = xGridMax;
        const real zGridMin    = xGridMin;
        const real zGridMax    = xGridMax;
        const real dSphere     = 10.0;
        const real dSphereLev1 = 20.0;

        if (para->getNumprocs() > 1) {
            const uint generatePart = vf::gpu::Communicator::getInstanz()->getPID();

            real overlap = (real)8.0 * dxGrid;
            gridBuilder->setNumberOfLayers(10, 8);

            if (comm->getNummberOfProcess() == 2) {
                const real ySplit = 0.0;

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, ySplit + overlap, zGridMax,
                                               dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid(xGridMin, ySplit - overlap, zGridMin, xGridMax, yGridMax, zGridMax,
                                               dxGrid);
                }

                if (useLevels) {
                    gridBuilder->addGrid(new Sphere(0.0, 0.0, 0.0, dSphere), 1);
                }

                gridBuilder->addGeometry(new Sphere(0.0, 0.0, 0.0, dSphereLev1));

                if (generatePart == 0)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, yGridMin, ySplit, zGridMin, zGridMax));
                if (generatePart == 1)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xGridMax, ySplit, yGridMax, zGridMin, zGridMax));

                gridBuilder->setPeriodicBoundaryCondition(false, false, true);

                gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

                if (generatePart == 0) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 1);
                }

                if (generatePart == 1) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 0);
                }

                //////////////////////////////////////////////////////////////////////////
                gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
                if (generatePart == 0)
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
                if (generatePart == 1)
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
                //////////////////////////////////////////////////////////////////////////
           
            } else if (comm->getNummberOfProcess() == 4) {

                const real ySplit = 0.0;
                const real xSplit = 0.0;

                if (generatePart == 0) {
                    gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xSplit + overlap, ySplit + overlap,
                                               zGridMax, dxGrid);
                }
                if (generatePart == 1) {
                    gridBuilder->addCoarseGrid(xGridMin, ySplit - overlap, zGridMin, xSplit + overlap, yGridMax,
                                               zGridMax, dxGrid);
                }
                if (generatePart == 2) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, yGridMin, zGridMin, xGridMax, ySplit + overlap,
                                               zGridMax, dxGrid);
                }
                if (generatePart == 3) {
                    gridBuilder->addCoarseGrid(xSplit - overlap, ySplit - overlap, zGridMin, xGridMax, yGridMax,
                                               zGridMax, dxGrid);
                }

                if (useLevels) {
                    gridBuilder->addGrid(new Sphere(0.0, 0.0, 0.0, dSphere), 1);
                }

                gridBuilder->addGeometry(new Sphere(0.0, 0.0, 0.0, dSphereLev1));

                if (generatePart == 0)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xSplit, yGridMin, ySplit, zGridMin, zGridMax));
                if (generatePart == 1)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xGridMin, xSplit, ySplit, yGridMax, zGridMin, zGridMax));
                if (generatePart == 2)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xSplit, xGridMax, yGridMin, ySplit, zGridMin, zGridMax));
                if (generatePart == 3)
                    gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xSplit, xGridMax, ySplit, yGridMax, zGridMin, zGridMax));

                gridBuilder->setPeriodicBoundaryCondition(false, false, true);

                gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

                if (generatePart == 0) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 1);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 2);
                }
                if (generatePart == 1) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 0);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 3);
                }
                if (generatePart == 2) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::PY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 3);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 0);
                }
                if (generatePart == 3) {
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MY, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 2);
                    gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                    gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 1);
                }

                //////////////////////////////////////////////////////////////////////////
                if (generatePart == 0) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
                }
                if (generatePart == 1) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                }
                if (generatePart == 2) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
                }
                if (generatePart == 3) {
                    gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
                    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
                }
                gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
                //////////////////////////////////////////////////////////////////////////
            }
            if (para->getKernelNeedsFluidNodeIndicesToRun())
                gridBuilder->findFluidNodes(useStreams);

            // gridBuilder->writeGridsToVtk(outPath + "/" + "/grid/part" +
            // std::to_string(generatePart) + "_"); gridBuilder->writeGridsToVtk(outPath + "/" + "/" +
            // std::to_string(generatePart) + "/grid/"); gridBuilder->writeArrows(outPath + "/" +  "/" +
            // std::to_string(generatePart) + " /arrow");

            SimulationFileWriter::write(gridPath + "/" + std::to_string(generatePart) + "/", gridBuilder,
                                        FILEFORMAT::BINARY);
        } else {

            gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, yGridMax, zGridMax, dxGrid);

            if (useLevels) {
                gridBuilder->setNumberOfLayers(10, 8);
                gridBuilder->addGrid(new Sphere(0.0, 0.0, 0.0, dSphereLev1), 1);
            }

            gridBuilder->addGeometry(new Sphere(0.0, 0.0, 0.0, dSphere));

            gridBuilder->setPeriodicBoundaryCondition(false, false, true);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            //////////////////////////////////////////////////////////////////////////
            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, vxLB, 0.0, 0.0);
            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
            //////////////////////////////////////////////////////////////////////////
            if (para->getKernelNeedsFluidNodeIndicesToRun())
                gridBuilder->findFluidNodes(useStreams);

            // gridBuilder->writeGridsToVtk("E:/temp/MusselOyster/" + "/grid/");
            // gridBuilder->writeArrows ("E:/temp/MusselOyster/" + "/arrow");

            SimulationFileWriter::write(gridPath, gridBuilder, FILEFORMAT::BINARY);
        }        
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
        
        try
        {
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

            if (configFile.size()==0) {
                configFile = targetPath + "config.txt";
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

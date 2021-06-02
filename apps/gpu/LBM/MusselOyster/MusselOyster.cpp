
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

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"

#include "Core/LbmOrGks.h"
#include "Core/Input/Input.h"
#include "Core/StringUtilities/StringUtil.h"
#include "Core/Input/ConfigFileReader/ConfigFileReader.h"

#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"

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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          U s e r    s e t t i n g s
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//const real L  = 1.0;
//
//const real Re = 500.0;// 1000.0;
//
//const real velocity  = 1.0;
//
//const real dt = (real)1.0e-3; //0.5e-3;
//
//const uint nx = 64;

std::string path("E:/temp/MusselOyster");

std::string simulationName("MusselOysterChim");

//const uint timeStepOut = 10000;
//const uint timeStepEnd = 250000;

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
	SPtr<ConfigFileReader> configReader = ConfigFileReader::getNewInstance();

    std::cout << configPath << std::endl;
	SPtr<ConfigData> configData = configReader->readConfigFile(configPath.c_str());

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//real dx = L / real(nx);

	//gridBuilder->addCoarseGrid(-0.5 * L, -0.5 * L, -0.5 * L,
	//							0.5 * L,  0.5 * L,  0.5 * L, dx);

	//gridBuilder->setPeriodicBoundaryCondition(false, false, false);

	//gridBuilder->buildGrids(lbmOrGks, false); // buildGrids() has to be called before setting the BCs!!!!

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        SPtr<Parameter>    para         = Parameter::make(configData, comm);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        real dx = 1.0;
        real vx = (real) 0.005;
        real Re = 100;

        para->setVelocity(vx);
        para->setViscosity((vx * dx) / Re);
        para->setVelocityRatio(1.0);

        para->setTOut(50000);
        para->setTEnd(250000);

        para->setCalcDragLift(false);
        para->setUseWale(false);

        para->setMainKernel("CumulantK15Comp");

        TriangularMesh *musselSTL =
            TriangularMesh::make("C:/Users/Master/Documents/MasterAnna/STL/MUSSEL_Paraview.stl");
        TriangularMesh *musselRef_1_STL =
            TriangularMesh::make("C:/Users/Master/Documents/MasterAnna/STL/MUSSEL_Level1.stl");

        //bounding box mussel:
        //x = -18, 58
        //y = -17, 18    
        //z = -5, 13

        const real f = 3.0;
        gridBuilder->addCoarseGrid(-18.0 * f, -17 * f, -5 * f, 
                                    58 * f * 2, 18 * f, 13 * f, dx); 

        gridBuilder->setNumberOfLayers(6, 8);
        gridBuilder->addGrid(musselRef_1_STL, 1);

        // gridBuilder->addGrid(sphereRef_2_STL, 4);
        // gridBuilder->setNumberOfLayers(10,8);

        gridBuilder->addGeometry(musselSTL);

        gridBuilder->setPeriodicBoundaryCondition(false, false, true);

        gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!
        //////////////////////////////////////////////////////////////////////////
        gridBuilder->setVelocityBoundaryCondition(SideType::PY, vx, 0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::MY, vx, 0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx, 0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx, 0.0, 0.0);

        gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::MX, vx, 0.0, 0.0);

        gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

        //////////////////////////////////////////////////////////////////////////
        SPtr<Grid> grid = gridBuilder->getGrid(gridBuilder->getNumberOfLevels() - 1);
        //////////////////////////////////////////////////////////////////////////

        gridBuilder->writeGridsToVtk("E:/temp/MusselOyster/grid/");
        // gridBuilder->writeArrows    ("F:/Work/Computations/out/Sphere/arrow");

        SimulationFileWriter::write("E:/temp/MusselOyster/grid/", gridBuilder, FILEFORMAT::BINARY);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //const real velocityLB = velocity * dt / dx; // LB units

	    //const real vx = velocityLB / (real)sqrt(2.0); // LB units
	    //const real vy = velocityLB / (real)sqrt(2.0); // LB units

        //const real viscosityLB = nx * velocityLB / Re; // LB units

        //*logging::out << logging::Logger::INFO_HIGH << "velocity  [dx/dt] = " << velocityLB << " \n";
        //*logging::out << logging::Logger::INFO_HIGH << "viscosity [dx^2/dt] = " << viscosityLB << "\n";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		para->setDevices(std::vector<uint>{(uint)0});

        para->setOutputPath( path );
        para->setOutputPrefix( simulationName );

        para->setFName(para->getOutputPath() + "/" + para->getOutputPrefix());

        para->setPrintFiles(true);

        para->setMaxLevel(1);

        //para->setVelocity(velocityLB);
        //para->setViscosity(viscosityLB);

        //para->setVelocityRatio(velocity/ velocityLB);

		//para->setMainKernel("CumulantK17CompChim");

		//para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
  //          rho = (real)0.0;
  //          vx  = (real)0.0; //(6 * velocityLB * coordZ * (L - coordZ) / (L * L));
  //          vy  = (real)0.0;
  //          vz  = (real)0.0;
  //      });


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
		//gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
		//gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
	    //gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
	    //gridBuilder->setVelocityBoundaryCondition(SideType::PZ,  vx,  vx, 0.0);
	    //gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        SPtr<CudaMemoryManager> cudaMemoryManager = CudaMemoryManager::make(para);

        SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager);

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
    std::string str, str2; 
    if ( argv != NULL )
    {
        //str = static_cast<std::string>(argv[0]);
        
        try
        {
            //////////////////////////////////////////////////////////////////////////

			std::string targetPath;

			targetPath = __FILE__;

#ifdef _WIN32
            targetPath = targetPath.substr(0, targetPath.find_last_of('\\') + 1);
#else
            targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);
#endif

			std::cout << targetPath << std::endl;

			multipleLevel(targetPath + "configMusselOyster.txt");

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

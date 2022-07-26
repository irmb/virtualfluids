//#define MPI_LOGGING

//Martin Branch

#include <mpi.h>
#if defined( MPI_LOGGING )
	#include <mpe.h>
#endif

#include <string>
#include <iostream>
#include <stdexcept>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Core/LbmOrGks.h"
#include "Core/StringUtilities/StringUtil.h"
#include "basics/config/ConfigurationFile.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

#include "global.h"

#include "geometries/Sphere/Sphere.h"
#include "geometries/VerticalCylinder/VerticalCylinder.h"
#include "geometries/Cuboid/Cuboid.h"
#include "geometries/TriangularMesh/TriangularMesh.h"
#include "geometries/Conglomerate/Conglomerate.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"

#include "grid/GridBuilder/LevelGridBuilder.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/GridFactory.h"

#include "io/SimulationFileWriter/SimulationFileWriter.h"
#include "io/GridVTKWriter/GridVTKWriter.h"
#include "io/STLReaderWriter/STLReader.h"
#include "io/STLReaderWriter/STLWriter.h"

#include "utilities/math/Math.h"
#include "utilities/communication.h"
#include "utilities/transformator/TransformatorImp.h"


void multipleLevel(const std::string& configPath)
{
    std::ofstream logFile( "F:/Work/Computations/out/Sphere/gridGeneratorLog.txt" );
    //std::ofstream logFile( "grid/gridGeneratorLog.txt" );
    logging::Logger::addStream(&logFile);

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    //UbLog::reportingLevel() = UbLog::logLevelFromString("DEBUG5");

    auto gridFactory = GridFactory::make();
    //gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::RAYCASTING);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    //gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_UNDER_TRIANGLE);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    
    vf::gpu::Communicator& communicator = vf::gpu::Communicator::getInstance();
    vf::basics::ConfigurationFile config;
    config.load(configPath);
    SPtr<Parameter> para = std::make_shared<Parameter>(config, communicator.getNummberOfProcess(), communicator.getPID());
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool useGridGenerator = true;

    if(useGridGenerator){

        enum testCase{
			TGV,
			TGV3D,
			SphereTest,
			DrivAer,
            PaperPlane,
            DLC,
            MultiGPU,
            StlGroupTest
        };

        int testcase = SphereTest;
        
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (testcase == TGV)
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		{
			real dx = 1.0;
			// real vx = 0.049;
			//////////////////////////////////////////////////////////////////////////
			//32
			gridBuilder->addCoarseGrid(-24, -2, -16,
										24,  2,  16, dx);
			//////////////////////////////////////////////////////////////////////////
			gridBuilder->setPeriodicBoundaryCondition(true, true, true);
			//////////////////////////////////////////////////////////////////////////
			gridBuilder->buildGrids(LBM, true); 
			//////////////////////////////////////////////////////////////////////////
			SPtr<Grid> grid = gridBuilder->getGrid(gridBuilder->getNumberOfLevels() - 1);
			//////////////////////////////////////////////////////////////////////////
			gridBuilder->writeGridsToVtk("E:/temp/TaylorGreenVortex/results/32/TGV32turned_Grid");
			gridBuilder->writeArrows("E:/temp/TaylorGreenVortex/results/32/TGV32turned_Grid_arrow");
			//////////////////////////////////////////////////////////////////////////
			SimulationFileWriter::write("E:/temp/TaylorGreenVortex/grids/turned/gridUni48x4x32/", gridBuilder, FILEFORMAT::BINARY);
			//////////////////////////////////////////////////////////////////////////
			return;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (testcase == TGV3D)
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		{
			const real PI = 3.141592653589793238462643383279;

			real dx = 2.0 * PI / 32.0; // 32^3 nodes
			//real dx = 2.0 * PI / 64.0; // 64^3 nodes
			//real dx = 2.0 * PI / 128.0; // 128^3 nodes
			//real dx = 2.0 * PI / 256.0; // 128^3 nodes
			// real vx = 0.049;

			gridBuilder->addCoarseGrid(-PI, -PI, -PI,
										PI,  PI,  PI, dx);

			gridBuilder->setPeriodicBoundaryCondition(true, true, true);

			gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

			//////////////////////////////////////////////////////////////////////////
			SPtr<Grid> grid = gridBuilder->getGrid(gridBuilder->getNumberOfLevels() - 1);
			//////////////////////////////////////////////////////////////////////////
			//32
			gridBuilder->writeGridsToVtk("E:/temp/TaylorGreenVortex/results3D/32/TGV3D_Grid");
			gridBuilder->writeArrows("E:/temp/TaylorGreenVortex/results3D/32/TGV3D_Grid_arrow");
			SimulationFileWriter::write("E:/temp/TaylorGreenVortex/grids3D/gridTGV3D/32/", gridBuilder, FILEFORMAT::BINARY); //FILEFORMAT::ASCII
			//256
		    //gridBuilder->writeGridsToVtk("E:/temp/TaylorGreenVortex/results3D/256/TGV3D_Grid");
			//gridBuilder->writeArrows("E:/temp/TaylorGreenVortex/results3D/256/TGV3D_Grid_arrow");
			//SimulationFileWriter::write("E:/temp/TaylorGreenVortex/grids3D/gridTGV3D/256/", gridBuilder, FILEFORMAT::BINARY); //FILEFORMAT::ASCII

			return;

		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( testcase == SphereTest)
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        {
            real dx = 1.0 / 16.0;
            real vx = 0.05;

            real D = 1.0;
            real Re = 100;

            para->setOutputPath( "F:/Work/Computations/out/Sphere/" );
            para->setOutputPrefix( "Sphere" );

            para->setPathAndFilename(para->getOutputPath() + "/" + para->getOutputPrefix());

            para->setPrintFiles(true);
    
            para->setVelocityLB( vx );
            para->setViscosityLB( ( vx * D / dx ) / Re );

            para->setVelocityRatio(1.0);

            para->setTimestepOut( 1000 );
            para->setTimestepEnd( 100000 );

            para->setCalcDragLift(true);

            para->setUseWale(false);

            para->setMainKernel("CumulantK15Comp");

            //////////////////////////////////////////////////////////////////////////

            // TriangularMesh* sphereSTL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/Sphere/SphereNotOptimal.stl");

            TriangularMesh* sphereRef_1_STL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/Sphere/SphereRef_1.stl");

            // TriangularMesh* sphereRef_2_STL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/Sphere/SphereRef_2.stl");

            Object* sphere = new Sphere( 0, 0, 0, 0.5*D );

            gridBuilder->addCoarseGrid(-2.0*D, -2.5*D, -2.5*D,
                                        9.0*D,  2.5*D,  2.5*D, dx);  // DrivAer

            //gridBuilder->setNumberOfLayers(10,8);
            //gridBuilder->addGrid(SphereSTL, 2);
            
            gridBuilder->setNumberOfLayers(4,8);
            gridBuilder->addGrid(sphereRef_1_STL, 1);
            //gridBuilder->addGrid(sphereRef_2_STL, 4);

            //gridBuilder->setNumberOfLayers(10,8);
            //gridBuilder->addGrid(sphere, 5);


        
            //gridBuilder->addGeometry(SphereSTL);
            gridBuilder->addGeometry(sphere);

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!
            //////////////////////////////////////////////////////////////////////////
            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx , 0.0, 0.0);

            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vx, 0.0, 0.0);

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
            
            bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
            bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
            bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

            //////////////////////////////////////////////////////////////////////////
            SPtr<Grid> grid = gridBuilder->getGrid(gridBuilder->getNumberOfLevels() - 1);
            //////////////////////////////////////////////////////////////////////////

            //gridBuilder->writeGridsToVtk("F:/Work/Computations/out/Sphere/grid");
            //gridBuilder->writeArrows    ("F:/Work/Computations/out/Sphere/arrow");

            SimulationFileWriter::write("F:/Work/Computations/out/Sphere/grid/", gridBuilder, FILEFORMAT::BINARY);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( testcase == DrivAer )
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        {

            real dx = 0.2;
            real vx = 0.05;

            real L = 4.6;
            real Re = 1.0e6;

            para->setOutputPath( "F:/Work/Computations/out/DrivAerNew/" );
            para->setOutputPrefix( "DrivAer" );

            para->setPathAndFilename(para->getOutputPath() + "/" + para->getOutputPrefix());

            para->setPrintFiles(true);
    
            para->setVelocityLB( vx );
            para->setViscosityLB( ( vx * L / dx ) / Re );

            //para->setVelocityRatio(1.0 / velocityLB);
            para->setVelocityRatio(1.0);

            para->setTimestepOut( 10000 );
            para->setTimestepEnd( 100000 );

            para->setUseWale(false);

            para->setMainKernel("CumulantK20Comp");

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            TriangularMesh* DrivAerSTL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/DrivAer_Fastback_Coarse.stl");
            //TriangularMesh* triangularMesh = TriangularMesh::make("M:/TestGridGeneration/STL/DrivAer_NoSTLGroups.stl");
            //TriangularMesh* triangularMesh = TriangularMesh::make("M:/TestGridGeneration/STL/DrivAer_Coarse.stl");
            //TriangularMesh* DrivAerSTL = TriangularMesh::make("stl/DrivAer_Fastback_Coarse.stl");

            TriangularMesh* DrivAerRefBoxSTL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/DrivAer_REF_BOX_Adrea.stl");
            //TriangularMesh* DrivAerRefBoxSTL = TriangularMesh::make("stl/DrivAer_REF_BOX_Adrea.stl");

            real z0 = 0.318;

            gridBuilder->addCoarseGrid(- 5.0, -5.0, 0.0 - z0,
                                        15.0,  5.0, 5.0 - z0, dx);  // DrivAer

            //Object* floorBox = new Cuboid( -0.3, -1, -1, 4.0, 1, 0.2 );
            //Object* wakeBox  = new Cuboid(  3.5, -1, -1, 5.5, 1, 0.8 );

            //Conglomerate* refRegion = new Conglomerate();

            //refRegion->add(floorBox);
            //refRegion->add(wakeBox);
            //refRegion->add(DrivAerRefBoxSTL);

            gridBuilder->setNumberOfLayers(10,8);
            gridBuilder->addGrid(DrivAerRefBoxSTL, 4);
        
            gridBuilder->setNumberOfLayers(10,8);
            gridBuilder->addGrid(DrivAerSTL, 5);

            gridBuilder->addGeometry(DrivAerSTL);

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx , 0.0, 0.0);

            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vx, 0.0, 0.0);

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

            bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
            bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
            bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);

            //////////////////////////////////////////////////////////////////////////

            SPtr<Grid> grid = gridBuilder->getGrid(gridBuilder->getNumberOfLevels() - 1);

            gridBuilder->getGeometryBoundaryCondition(gridBuilder->getNumberOfLevels() - 1)->setTangentialVelocityForPatch( grid, 4, 0.0075, -2.0, 0.0,
                                                                                                                                     0.0075,  2.0, 0.0, -vx, 0.318);
            gridBuilder->getGeometryBoundaryCondition(gridBuilder->getNumberOfLevels() - 1)->setTangentialVelocityForPatch( grid, 3, 2.793 , -2.0, 0.0,
                                                                                                                                     2.793 ,  2.0, 0.0, -vx, 0.318);

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->writeGridsToVtk("F:/Work/Computations/out/DrivAerNew/DrivAer_Grid_");
            gridBuilder->writeArrows    ("F:/Work/Computations/out/DrivAerNew/DrivAer_Grid_arrow");

            //SimulationFileWriter::write("D:/GRIDGENERATION/files/", gridBuilder, FILEFORMAT::ASCII);
            //SimulationFileWriter::write("C:/Users/lenz/Desktop/Work/gridGenerator/grid/", gridBuilder, FILEFORMAT::ASCII);
            SimulationFileWriter::write("F:/Work/Computations/out/DrivAerNew/grid/", gridBuilder, FILEFORMAT::BINARY);
            //SimulationFileWriter::write("grid/", gridBuilder, FILEFORMAT::ASCII);

            //return;
            //gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( testcase == PaperPlane )
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        {

            real dx = 0.01;
            real vx = 0.05;

            real L = 0.3;
            real Re = 90000;

            para->setOutputPath( "F:/Work/Computations/out/PaperPlane/" );
            para->setOutputPrefix( "PaperPlaneK17winglet" );

            para->setPathAndFilename(para->getOutputPath() + "/" + para->getOutputPrefix());

            para->setPrintFiles(true);
    
            para->setVelocityLB( vx );
            para->setViscosityLB( ( vx * L / dx ) / Re );

            para->setVelocityRatio(1.0);

            para->setTimestepOut( 1000 );
            para->setTimestepEnd( 100000 );

            para->setUseWale(false);

            para->setMainKernel("CumulantAA2016CompSP27");

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            TriangularMesh* STL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/PaperPlane_1.stl");
            //TriangularMesh* STL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/PaperPlane_1_winglet.stl");
            
            TriangularMesh* RefBoxSTL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/PaperPlane_1_ref.stl");
            //TriangularMesh* RefBoxSTL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/PaperPlane_1_winglet_ref.stl");

            gridBuilder->addCoarseGrid(- 0.5, -0.3, -0.3,
                                         1.0,  0.3,  0.3, dx);

            gridBuilder->setNumberOfLayers(6,8);
            gridBuilder->addGrid(RefBoxSTL, 3);

            gridBuilder->addGeometry(STL);

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx , 0.0, 0.0);

            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vx, 0.0, 0.0);

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

            bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
            bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
            bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->writeGridsToVtk("F:/Work/Computations/out/PaperPlane/PaperPlane_Grid_");
            gridBuilder->writeArrows    ("F:/Work/Computations/out/PaperPlane/PaperPlane_Grid_arrow");

            //SimulationFileWriter::write("D:/GRIDGENERATION/files/", gridBuilder, FILEFORMAT::ASCII);
            //SimulationFileWriter::write("C:/Users/lenz/Desktop/Work/gridGenerator/grid/", gridBuilder, FILEFORMAT::ASCII);
            SimulationFileWriter::write("F:/Work/Computations/out/PaperPlane/grid/", gridBuilder, FILEFORMAT::BINARY);
            //SimulationFileWriter::write("grid/", gridBuilder, FILEFORMAT::ASCII);

            //return;
            //gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( testcase == StlGroupTest )
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        {

            real dx = 0.025;
            real vx = 0.05;

            real L = 1.0;
            real Re = 100;

            para->setOutputPath( "F:/Work/Computations/out/StlGroupTest/" );
            para->setOutputPrefix( "StlGroupTest" );

            para->setPathAndFilename(para->getOutputPath() + "/" + para->getOutputPrefix());

            para->setPrintFiles(true);
    
            para->setVelocityLB( vx );
            para->setViscosityLB( ( vx * L / dx ) / Re );

            para->setVelocityRatio(1.0);

            para->setTimestepOut( 1000 );
            para->setTimestepEnd( 100000 );

            para->setUseWale(false);

            para->setMainKernel("CumulantAA2016CompSP27");
            //para->setMainKernel(kernelMapper->getEnum("CumulantOneCompSP27"));

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            TriangularMesh* STL = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/STL_Group_Test_2_Cylinders.stl");

            gridBuilder->addCoarseGrid(- 2.0, -4.5, -2.0,
                                         4.0,  4.5,  2.0, dx);

            gridBuilder->addGeometry(STL);

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx , 0.0, 0.0);

            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vx, 0.0, 0.0);

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

            bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
            bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
            bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

            //////////////////////////////////////////////////////////////////////////

            SPtr<Grid> grid = gridBuilder->getGrid(gridBuilder->getNumberOfLevels() - 1);

            gridBuilder->getGeometryBoundaryCondition(gridBuilder->getNumberOfLevels() - 1)->setTangentialVelocityForPatch( grid, 1, 0.0, -2.0, 0.0,
                                                                                                                                     0.0,  2.0, 0.0, -vx, 0.5);

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->writeGridsToVtk("F:/Work/Computations/out/StlGroupTest/StlGroupTest_Grid_");
            gridBuilder->writeArrows    ("F:/Work/Computations/out/StlGroupTest/StlGroupTest_Grid_arrow");

            SimulationFileWriter::write("F:/Work/Computations/out/StlGroupTest/grid/", gridBuilder, FILEFORMAT::BINARY);

        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( testcase == DLC )
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        {
			real velocityRatio = 594.093427;

			real dx = 0.2;
			real vx = 0.065272188;

			real z0 = 0.24395 + 0.5*dx;

            std::vector<uint> ignorePatches = { 152, 153, 154 };

            //TriangularMesh* VW370_SERIE_STL = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/VW370_SERIE.stl", ignorePatches);
            TriangularMesh* VW370_SERIE_STL = TriangularMesh::make("stl/VW370_SERIE.stl", ignorePatches);

            //TriangularMesh* DLC_RefBox = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DLC_RefBox.stl");

            //TriangularMesh* DLC_RefBox_1 = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DLC_RefBox_withWake/DLC_RefBox_withWake_4m.stl");
            //TriangularMesh* DLC_RefBox_2 = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DLC_RefBox_withWake/DLC_RefBox_withWake_3m.stl");
            //TriangularMesh* DLC_RefBox_3 = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DLC_RefBox_withWake/DLC_RefBox_withWake_2m.stl");
            //TriangularMesh* DLC_RefBox_4 = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DLC_RefBox_withWake/DLC_RefBox_withWake_1m.stl");

            //TriangularMesh* DLC_RefBox_Level_3 = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DLC/DLC_RefBox_Level_3.stl");
            //TriangularMesh* DLC_RefBox_Level_4 = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DLC/DLC_RefBox_Level_4.stl");
            //TriangularMesh* DLC_RefBox_Level_5 = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DLC/DLC_RefBox_Level_5.stl");

            TriangularMesh* DLC_RefBox_Level_3 = TriangularMesh::make("stl/DLC/DLC_RefBox_Level_3.stl");
            TriangularMesh* DLC_RefBox_Level_4 = TriangularMesh::make("stl/DLC/DLC_RefBox_Level_4.stl");
            TriangularMesh* DLC_RefBox_Level_5 = TriangularMesh::make("stl/DLC/DLC_RefBox_Level_5.stl");

            //TriangularMesh* VW370_SERIE_STL = TriangularMesh::make("stl/VW370_SERIE.stl", ignorePatches);
            //TriangularMesh* DLC_RefBox = TriangularMesh::make("stl/DLC_RefBox.lnx.stl");
            //TriangularMesh* DLC_RefBox_4 = TriangularMesh::make("stl/DLC_RefBox_withWake/DLC_RefBox_withWake_1m.lnx.stl");

            gridBuilder->addCoarseGrid(-30.0, -20.0,  0.0 - z0,
                                        50.0,  20.0, 25.0 - z0, dx);
            
            gridBuilder->setNumberOfLayers(10,8);
            gridBuilder->addGrid( new Cuboid( - 6.6, -6, -0.7, 20.6 , 6, 5.3  ), 1 );
            gridBuilder->addGrid( new Cuboid( -3.75, -3, -0.7, 11.75, 3, 2.65 ), 2 );

            gridBuilder->setNumberOfLayers(10,8);
            gridBuilder->addGrid(DLC_RefBox_Level_3, 3);
            gridBuilder->addGrid(DLC_RefBox_Level_4, 4);
        
            Conglomerate* refinement = new Conglomerate();
            refinement->add(DLC_RefBox_Level_5);
            refinement->add(VW370_SERIE_STL);

            gridBuilder->setNumberOfLayers(10,8);
            gridBuilder->addGrid(refinement, 5);

            gridBuilder->addGeometry(VW370_SERIE_STL);

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx , 0.0, 0.0);

            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vx, 0.0, 0.0);

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

            bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
            bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
            bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);
 
            //////////////////////////////////////////////////////////////////////////

            SPtr<Grid> grid = gridBuilder->getGrid(gridBuilder->getNumberOfLevels() - 1);

            real wheelsFrontX = -0.081;
            real wheelsRearX  =  2.5486;

            real wheelsFrontZ =  0.0504;
            real wheelsRearZ  =  0.057;

            real wheelsRadius =  0.318;

			real wheelRotationFrequency = 1170.74376 / 60.0;

			real wheelTangentialVelocity = -2.0 * M_PI * wheelsRadius * wheelRotationFrequency / velocityRatio;

            std::vector<uint> frontWheelPatches = { 71, 86, 87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97, 159 };
            std::vector<uint> rearWheelPatches  = { 82, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 160 };

            for( uint patch : frontWheelPatches ){
                gridBuilder->getGeometryBoundaryCondition(gridBuilder->getNumberOfLevels() - 1)->setTangentialVelocityForPatch( grid, patch, wheelsFrontX, -2.0, wheelsFrontZ,
                                                                                                                                             wheelsFrontX,  2.0, wheelsFrontZ, 
					                                                                                                                         wheelTangentialVelocity, wheelsRadius);
            }

            for( uint patch : rearWheelPatches ){
                gridBuilder->getGeometryBoundaryCondition(gridBuilder->getNumberOfLevels() - 1)->setTangentialVelocityForPatch( grid, patch, wheelsRearX , -2.0, wheelsRearZ ,
                                                                                                                                             wheelsRearX ,  2.0, wheelsRearZ , 
					                                                                                                                         wheelTangentialVelocity, wheelsRadius);
            }

            //////////////////////////////////////////////////////////////////////////

            //gridBuilder->writeGridsToVtk("C:/Users/lenz/Desktop/Work/gridGenerator/grid/DLC_Grid");
            //gridBuilder->writeArrows    ("C:/Users/lenz/Desktop/Work/gridGenerator/grid/DLC_Grid_arrow");

            gridBuilder->writeGridsToVtk("grid/DLC_Grid");
            gridBuilder->writeArrows    ("grid/DLC_Grid_arrow");

            //SimulationFileWriter::write("D:/GRIDGENERATION/files/", gridBuilder, FILEFORMAT::ASCII);
            //SimulationFileWriter::write("C:/Users/lenz/Desktop/Work/gridGenerator/grid/", gridBuilder, FILEFORMAT::ASCII);
            SimulationFileWriter::write("grid/", gridBuilder, FILEFORMAT::ASCII);
            
            //gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( testcase == MultiGPU )
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        {

            real dx = 1.0 / 40.0;
            real vx = 0.05;

            real D = 1.0;
            real Re = 100;

            para->setOutputPath( "F:/Work/Computations/out/Sphere/" );
            para->setOutputPrefix( "Sphere" );

            para->setPathAndFilename(para->getOutputPath() + "/" + para->getOutputPrefix());

            para->setPrintFiles(true);
    
            para->setVelocityLB( vx );
            para->setViscosityLB( ( vx * D / dx ) / Re );

            para->setVelocityRatio(1.0);

            para->setTimestepOut( 1000 );
            para->setTimestepEnd( 100000 );

            para->setCalcDragLift(true);

            para->setUseWale(false);

            para->setMainKernel("CumulantK15Comp");

            para->setDevices( { 0, 1 } );
            para->setMaxDev(2);

            //const uint generatePart = 1;
            const uint generatePart = communicator.getPID();
            
            std::ofstream logFile2;
            
            if( generatePart == 0 )
                logFile2.open( "F:/Work/Computations/gridGenerator/grid/0/gridGeneratorLog.txt" );
                //logFile2.open( "grid/0/gridGeneratorLog.txt" );
            
            if( generatePart == 1 )
                logFile2.open( "F:/Work/Computations/gridGenerator/grid/1/gridGeneratorLog.txt" );
                //logFile2.open( "grid/1/gridGeneratorLog.txt" );

            logging::Logger::addStream(&logFile2);

            TriangularMesh* triangularMesh = TriangularMesh::make("F:/Work/Computations/gridGenerator/stl/Sphere/SphereNotOptimal.stl");
            //TriangularMesh* triangularMesh = TriangularMesh::make("stl/ShpereNotOptimal.lnx.stl");

            // all
            //gridBuilder->addCoarseGrid(-2, -2, -2,  
            //                            4,  2,  2, dx);

            real overlap = 10.0 * dx;

            if( generatePart == 0 )
                gridBuilder->addCoarseGrid(-2.0          , -2.0, -2.0,  
                                            0.5 + overlap,  2.0,  2.0, dx);

            if( generatePart == 1 )
                gridBuilder->addCoarseGrid( 0.5 - overlap, -2.0, -2.0,  
                                            4.0          ,  2.0,  2.0, dx);


            gridBuilder->setNumberOfLayers(10,8);
            gridBuilder->addGrid(triangularMesh, 1);

            gridBuilder->addGeometry(triangularMesh);
            
            if( generatePart == 0 )
                gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>( -2.0, 0.5, 
                                                                             -2.0, 2.0, 
                                                                             -2.0, 2.0 ) );
            
            if( generatePart == 1 )
                gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>(  0.5, 4.0, 
                                                                             -2.0, 2.0, 
                                                                             -2.0, 2.0 ) );

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!
            
            if( generatePart == 0 ){
                gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 1);
            }
            
            if( generatePart == 1 ){
                gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 0);
            }

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx , 0.0, 0.0);

            if (generatePart == 0) {
                gridBuilder->setVelocityBoundaryCondition(SideType::MX, vx, 0.0, 0.0);
            }
            if (generatePart == 1) {
                gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
            }

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

            bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
            bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
            bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

            //////////////////////////////////////////////////////////////////////////

            if (generatePart == 0) {
                //gridBuilder->writeGridsToVtk("F:/Work/Computations/gridGenerator/grid/0/Test_");
                //gridBuilder->writeArrows    ("F:/Work/Computations/gridGenerator/grid/0/Test_Arrow");
            }
            if (generatePart == 1) {
                //gridBuilder->writeGridsToVtk("F:/Work/Computations/gridGenerator/grid/1/Test_");
                //gridBuilder->writeArrows    ("F:/Work/Computations/gridGenerator/grid/1/Test_Arrow");
            }

            if (generatePart == 0)
                SimulationFileWriter::write("F:/Work/Computations/gridGenerator/grid/0/", gridBuilder, FILEFORMAT::ASCII);
                //SimulationFileWriter::write("grid/0/", gridBuilder, FILEFORMAT::ASCII);
            if (generatePart == 1)
                SimulationFileWriter::write("F:/Work/Computations/gridGenerator/grid/1/", gridBuilder, FILEFORMAT::ASCII);
                //SimulationFileWriter::write("grid/1/", gridBuilder, FILEFORMAT::ASCII);

            //return;
            
            //gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);
        }

    }
    else
    {
        //gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);
        //gridGenerator = GridProvider::makeGridReader(FILEFORMAT::BINARY, para, cudaMemoryManager);
    }

    logFile.close();

    //return;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    SPtr<GridProvider> gridGenerator;
    if( useGridGenerator ) gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);
    else                   gridGenerator = GridProvider::makeGridReader(FILEFORMAT::BINARY, para, cudaMemoryManager);

    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory);
    sim.run();
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

			multipleLevel(targetPath + "config.txt");

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

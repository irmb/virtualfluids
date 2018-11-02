//#define MPI_LOGGING


#include <mpi.h>
#if defined( MPI_LOGGING )
	#include <mpe.h>
#endif

#include <string>
#include <iostream>

#include "LBM/Simulation.h"

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "VirtualFluidsBasics/utilities/input/Input.h"
#include "VirtualFluidsBasics/utilities/StringUtil/StringUtil.h"
#include "grid/GridBuilder/LevelGridBuilder.h"
#include "utilities/transformator/TransformatorImp.h"
#include "io/GridVTKWriter/GridVTKWriter.h"

#include "io/SimulationFileWriter/SimulationFileWriter.h"
#include "grid/GridBuilder/LevelGridBuilder.h"
#include "grid/GridBuilder/ParallelGridBuilder.h"
#include "geometries/TriangularMesh/TriangularMesh.h"

#include "grid/GridFactory.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"
#include <grid/GridMocks.h>
#include "grid/GridStrategy/GridStrategyMocks.h"
#include "VirtualFluidsBasics/utilities/logger/Logger.h"
#include "geometries/Conglomerate/Conglomerate.h"
#include "io/STLReaderWriter/STLReader.h"
#include "io/STLReaderWriter/STLWriter.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"
#include "Output/FileWriter.h"


std::string getGridPath(std::shared_ptr<Parameter> para, std::string Gridpath)
{
    if (para->getNumprocs() == 1)
        return Gridpath + "/";
    
    return Gridpath + "/" + StringUtil::toString(para->getMyID()) + "/";
}

void setParameters(std::shared_ptr<Parameter> para, std::unique_ptr<input::Input> &input)
{
    std::string _path = input->getValue("Path");
    std::string _prefix = input->getValue("Prefix");
    std::string _gridpath = input->getValue("GridPath");
    para->setNumprocs(1);
    std::string gridPath = getGridPath(para, _gridpath);
    para->setMaxDev(StringUtil::toInt(input->getValue("NumberOfDevices")));
    para->setDevices(StringUtil::toVector(input->getValue("Devices")));
    para->setOutputPath(_path);
    para->setOutputPrefix(_prefix);
    para->setFName(_path + "/" + _prefix);
    para->setPrintFiles(false);
    para->setPrintFiles(StringUtil::toBool(input->getValue("WriteGrid")));
    para->setGeometryValues(StringUtil::toBool(input->getValue("GeometryValues")));
    para->setCalc2ndOrderMoments(StringUtil::toBool(input->getValue("calc2ndOrderMoments")));
    para->setCalc3rdOrderMoments(StringUtil::toBool(input->getValue("calc3rdOrderMoments")));
    para->setCalcHighOrderMoments(StringUtil::toBool(input->getValue("calcHigherOrderMoments")));
    para->setReadGeo(StringUtil::toBool(input->getValue("ReadGeometry")));
    para->setCalcMedian(StringUtil::toBool(input->getValue("calcMedian")));
    para->setConcFile(StringUtil::toBool(input->getValue("UseConcFile")));
    para->setUseMeasurePoints(StringUtil::toBool(input->getValue("UseMeasurePoints")));
    para->setUseWale(StringUtil::toBool(input->getValue("UseWale")));
    para->setSimulatePorousMedia(StringUtil::toBool(input->getValue("SimulatePorousMedia")));
    para->setD3Qxx(StringUtil::toInt(input->getValue("D3Qxx")));
    para->setTEnd(StringUtil::toInt(input->getValue("TimeEnd")));
    para->setTOut(StringUtil::toInt(input->getValue("TimeOut")));
    para->setTStartOut(StringUtil::toInt(input->getValue("TimeStartOut")));
    para->setTimeCalcMedStart(StringUtil::toInt(input->getValue("TimeStartCalcMedian")));
    para->setTimeCalcMedEnd(StringUtil::toInt(input->getValue("TimeEndCalcMedian")));
    para->setPressInID(StringUtil::toInt(input->getValue("PressInID")));
    para->setPressOutID(StringUtil::toInt(input->getValue("PressOutID")));
    para->setPressInZ(StringUtil::toInt(input->getValue("PressInZ")));
    para->setPressOutZ(StringUtil::toInt(input->getValue("PressOutZ")));
    //////////////////////////////////////////////////////////////////////////
	para->setCompOn(StringUtil::toBool(input->getValue("CompOn")));
	para->setDiffOn(StringUtil::toBool(input->getValue("DiffOn")));
    para->setDiffMod(StringUtil::toInt(input->getValue("DiffMod")));
    para->setDiffusivity(StringUtil::toFloat(input->getValue("Diffusivity")));
    para->setTemperatureInit(StringUtil::toFloat(input->getValue("Temp")));
    para->setTemperatureBC(StringUtil::toFloat(input->getValue("TempBC")));
    //////////////////////////////////////////////////////////////////////////
    para->setViscosity(StringUtil::toFloat(input->getValue("Viscosity_LB")));
    para->setVelocity(StringUtil::toFloat(input->getValue("Velocity_LB")));
    para->setViscosityRatio(StringUtil::toFloat(input->getValue("Viscosity_Ratio_World_to_LB")));
    para->setVelocityRatio(StringUtil::toFloat(input->getValue("Velocity_Ratio_World_to_LB")));
    para->setDensityRatio(StringUtil::toFloat(input->getValue("Density_Ratio_World_to_LB")));
    para->setPressRatio(StringUtil::toFloat(input->getValue("Delta_Press")));
    para->setRealX(StringUtil::toFloat(input->getValue("SliceRealX")));
    para->setRealY(StringUtil::toFloat(input->getValue("SliceRealY")));
    para->setFactorPressBC(StringUtil::toFloat(input->getValue("dfpbc")));
    para->setGeometryFileC(input->getValue("GeometryC"));
    para->setGeometryFileM(input->getValue("GeometryM"));
    para->setGeometryFileF(input->getValue("GeometryF"));
    //////////////////////////////////////////////////////////////////////////
    para->setgeoVec(gridPath + input->getValue("geoVec"));
    para->setcoordX(gridPath + input->getValue("coordX"));
    para->setcoordY(gridPath + input->getValue("coordY"));
    para->setcoordZ(gridPath + input->getValue("coordZ"));
    para->setneighborX(gridPath + input->getValue("neighborX"));
    para->setneighborY(gridPath + input->getValue("neighborY"));
    para->setneighborZ(gridPath + input->getValue("neighborZ"));
    para->setscaleCFC(gridPath + input->getValue("scaleCFC"));
    para->setscaleCFF(gridPath + input->getValue("scaleCFF"));
    para->setscaleFCC(gridPath + input->getValue("scaleFCC"));
    para->setscaleFCF(gridPath + input->getValue("scaleFCF"));
    para->setscaleOffsetCF(gridPath + input->getValue("scaleOffsetCF"));
    para->setscaleOffsetFC(gridPath + input->getValue("scaleOffsetFC"));
    para->setgeomBoundaryBcQs(gridPath + input->getValue("geomBoundaryBcQs"));
    para->setgeomBoundaryBcValues(gridPath + input->getValue("geomBoundaryBcValues"));
    para->setinletBcQs(gridPath + input->getValue("inletBcQs"));
    para->setinletBcValues(gridPath + input->getValue("inletBcValues"));
    para->setoutletBcQs(gridPath + input->getValue("outletBcQs"));
    para->setoutletBcValues(gridPath + input->getValue("outletBcValues"));
    para->settopBcQs(gridPath + input->getValue("topBcQs"));
    para->settopBcValues(gridPath + input->getValue("topBcValues"));
    para->setbottomBcQs(gridPath + input->getValue("bottomBcQs"));
    para->setbottomBcValues(gridPath + input->getValue("bottomBcValues"));
    para->setfrontBcQs(gridPath + input->getValue("frontBcQs"));
    para->setfrontBcValues(gridPath + input->getValue("frontBcValues"));
    para->setbackBcQs(gridPath + input->getValue("backBcQs"));
    para->setbackBcValues(gridPath + input->getValue("backBcValues"));
    para->setnumberNodes(gridPath + input->getValue("numberNodes"));
    para->setLBMvsSI(gridPath + input->getValue("LBMvsSI"));
    //////////////////////////////gridPath + ////////////////////////////////////////////
    para->setmeasurePoints(gridPath + input->getValue("measurePoints"));
    para->setpropellerValues(gridPath + input->getValue("propellerValues"));
    para->setclockCycleForMP(StringUtil::toFloat(input->getValue("measureClockCycle")));
    para->settimestepForMP(StringUtil::toInt(input->getValue("measureTimestep")));
    para->setcpTop(gridPath + input->getValue("cpTop"));
    para->setcpBottom(gridPath + input->getValue("cpBottom"));
    para->setcpBottom2(gridPath + input->getValue("cpBottom2"));
    para->setConcentration(gridPath + input->getValue("Concentration"));
    //////////////////////////////////////////////////////////////////////////
    //Normals - Geometry
    para->setgeomBoundaryNormalX(gridPath + input->getValue("geomBoundaryNormalX"));
    para->setgeomBoundaryNormalY(gridPath + input->getValue("geomBoundaryNormalY"));
    para->setgeomBoundaryNormalZ(gridPath + input->getValue("geomBoundaryNormalZ"));
    //Normals - Inlet
    para->setInflowBoundaryNormalX(gridPath + input->getValue("inletBoundaryNormalX"));
    para->setInflowBoundaryNormalY(gridPath + input->getValue("inletBoundaryNormalY"));
    para->setInflowBoundaryNormalZ(gridPath + input->getValue("inletBoundaryNormalZ"));
    //Normals - Outlet
    para->setOutflowBoundaryNormalX(gridPath + input->getValue("outletBoundaryNormalX"));
    para->setOutflowBoundaryNormalY(gridPath + input->getValue("outletBoundaryNormalY"));
    para->setOutflowBoundaryNormalZ(gridPath + input->getValue("outletBoundaryNormalZ"));
    //////////////////////////////////////////////////////////////////////////
    //Forcing
    para->setForcing(StringUtil::toFloat(input->getValue("ForcingX")), StringUtil::toFloat(input->getValue("ForcingY")), StringUtil::toFloat(input->getValue("ForcingZ")));
    //////////////////////////////////////////////////////////////////////////
    //Particles
    para->setCalcParticles(StringUtil::toBool(input->getValue("calcParticles")));
    para->setParticleBasicLevel(StringUtil::toInt(input->getValue("baseLevel")));
    para->setParticleInitLevel(StringUtil::toInt(input->getValue("initLevel")));
    para->setNumberOfParticles(StringUtil::toInt(input->getValue("numberOfParticles")));
    para->setneighborWSB(gridPath + input->getValue("neighborWSB"));
    para->setStartXHotWall(StringUtil::toDouble(input->getValue("startXHotWall")));
    para->setEndXHotWall(StringUtil::toDouble(input->getValue("endXHotWall")));
    //////////////////////////////////////////////////////////////////////////
    //for Multi GPU
    if (para->getNumprocs() > 1)
    {
        ////////////////////////////////////////////////////////////////////////////
        ////1D domain decomposition
        //std::vector<std::string> sendProcNeighbors;
        //std::vector<std::string> recvProcNeighbors;
        //for (int i = 0; i<para->getNumprocs();i++)
        //{
        // sendProcNeighbors.push_back(gridPath + StringUtil::toString(i) + "s.dat");
        // recvProcNeighbors.push_back(gridPath + StringUtil::toString(i) + "r.dat");
        //}
        //para->setPossNeighborFiles(sendProcNeighbors, "send");
        //para->setPossNeighborFiles(recvProcNeighbors, "recv");
        //////////////////////////////////////////////////////////////////////////
        //3D domain decomposition
        std::vector<std::string> sendProcNeighborsX, sendProcNeighborsY, sendProcNeighborsZ;
        std::vector<std::string> recvProcNeighborsX, recvProcNeighborsY, recvProcNeighborsZ;
        for (int i = 0; i < para->getNumprocs(); i++)
        {
            sendProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xs.dat");
            sendProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Ys.dat");
            sendProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zs.dat");
            recvProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xr.dat");
            recvProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Yr.dat");
            recvProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zr.dat");
        }
        para->setPossNeighborFilesX(sendProcNeighborsX, "send");
        para->setPossNeighborFilesY(sendProcNeighborsY, "send");
        para->setPossNeighborFilesZ(sendProcNeighborsZ, "send");
        para->setPossNeighborFilesX(recvProcNeighborsX, "recv");
        para->setPossNeighborFilesY(recvProcNeighborsY, "recv");
        para->setPossNeighborFilesZ(recvProcNeighborsZ, "recv");
    }
    //////////////////////////////////////////////////////////////////////////
    //para->setkFull(             input->getValue( "kFull" ));
    //para->setgeoFull(           input->getValue( "geoFull" ));
    //para->setnoSlipBcPos(       input->getValue( "noSlipBcPos" ));
    //para->setnoSlipBcQs(          input->getValue( "noSlipBcQs" ));
    //para->setnoSlipBcValues(      input->getValue( "noSlipBcValues" ));
    //para->setnoSlipBcValue(     input->getValue( "noSlipBcValue" ));
    //para->setslipBcPos(         input->getValue( "slipBcPos" ));
    //para->setslipBcQs(          input->getValue( "slipBcQs" ));
    //para->setslipBcValue(       input->getValue( "slipBcValue" ));
    //para->setpressBcPos(        input->getValue( "pressBcPos" ));
    //para->setpressBcQs(           input->getValue( "pressBcQs" ));
    //para->setpressBcValues(       input->getValue( "pressBcValues" ));
    //para->setpressBcValue(      input->getValue( "pressBcValue" ));
    //para->setvelBcQs(             input->getValue( "velBcQs" ));
    //para->setvelBcValues(         input->getValue( "velBcValues" ));
    //para->setpropellerCylinder( input->getValue( "propellerCylinder" ));
    //para->setpropellerQs(		 input->getValue( "propellerQs"      ));
    //para->setwallBcQs(            input->getValue( "wallBcQs"         ));
    //para->setwallBcValues(        input->getValue( "wallBcValues"     ));
    //para->setperiodicBcQs(        input->getValue( "periodicBcQs"     ));
    //para->setperiodicBcValues(    input->getValue( "periodicBcValues" ));
    //cout << "Try this: " << para->getgeomBoundaryBcValues() << endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Restart
    para->setTimeDoCheckPoint(StringUtil::toInt(input->getValue("TimeDoCheckPoint")));
    para->setTimeDoRestart(StringUtil::toInt(input->getValue("TimeDoRestart")));
    para->setDoCheckPoint(StringUtil::toBool(input->getValue("DoCheckPoint")));
    para->setDoRestart(StringUtil::toBool(input->getValue("DoRestart")));
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*    para->setMaxLevel(StringUtil::toInt(input->getValue("NOGL")));
    para->setGridX(StringUtil::toVector(input->getValue("GridX")));                           
    para->setGridY(StringUtil::toVector(input->getValue("GridY")));                           
    para->setGridZ(StringUtil::toVector(input->getValue("GridZ")));                  
    para->setDistX(StringUtil::toVector(input->getValue("DistX")));                  
    para->setDistY(StringUtil::toVector(input->getValue("DistY")));                  
    para->setDistZ(StringUtil::toVector(input->getValue("DistZ")));      */            

    para->setNeedInterface(std::vector<bool>{true, true, true, true, true, true});

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Kernel
	para->setMainKernel(input->getValue("MainKernelName"));
	para->setMultiKernelOn(StringUtil::toBool(input->getValue("multiKernelOn")));
	para->setMultiKernelLevel(StringUtil::toVector(input->getValue("multiKernelLevel")));
	para->setMultiKernelName(StringUtil::toStringVector(input->getValue("multiKernelName")));
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}



void multipleLevel(const std::string& configPath)
{
    logging::Logger::setStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    auto gridFactory = SPtr<GridFactory>(new GridFactory());
    gridFactory->setGridStrategy(SPtr<GridStrategy>(new GridCpuStrategy()));
    gridFactory->setGrid("grid");
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::RAYCASTING);

    //auto gridBuilderlevel = LevelGridBuilder::makeShared(Device::CPU, "D3Q27");
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    //Conglomerate* conglomerate = new Conglomerate();
    //conglomerate->add(new Cuboid(10, 10, 10, 30, 30, 30));
    //conglomerate->subtract(new Sphere(30, 20, 20, 4));
    //gridBuilder->addGrid(conglomerate, 2);


//    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 14, 10, 20, 0.25);
    //TriangularMesh* triangularMesh = TriangularMesh::make("D:/GRIDGENERATION/STL/quadarBinaer.stl", DiscretizationMethod::POINT_IN_OBJECT);


    gridBuilder->addCoarseGrid(-10, -8, -3, 50, 20, 20, 0.25);
    TriangularMesh* triangularMesh = TriangularMesh::make("D:/GRIDGENERATION/STL/input/local_input/bruecke.stl", DiscretizationMethod::RAYCASTING);


    //TriangleOffsetSurfaceGeneration::createOffsetTriangularMesh(triangularMesh, 5);

    //TriangularMesh* sphere = TriangularMesh::make("D:/GRIDGENERATION/STL/GTI.stl", DiscretizationMethod::RAYCASTING);
    //TransformatorImp trans(1.0, Vertex(5.5, 1, 12));
    //trans.transformWorldToGrid(*sphere);
    //STLWriter::writeSTL(sphere->triangleVec, "D:/GRIDGENERATION/STL/GTI2.stl", false);

    //gridBuilder->addGrid(new Sphere(20, 20, 20, 8));
    gridBuilder->addGrid(triangularMesh, 2);

    //gridBuilder->addFineGrid(new Cuboid(15, 15, 15, 25, 25, 25), 1);
    //gridBuilder->addFineGrid(new Cuboid(17, 17, 17, 23, 23, 23), 2);


    //gridBuilder->addFineGrid(17.0, 17.0, 17.0, 20.0, 20.0, 20.0, 3);
    //gridBuilder->addFineGrid(10.0, 10.0, 10.0, 20.0, 20.0, 20.0, 3);


    //gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTest_level_2", 2);

    gridBuilder->buildGrids();

    gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTestSphere_level_0", 0);
    gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTestSphere_level_1", 1);
    gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTestSphere_level_2", 2);

    //gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTestCuboid_level_0", 0);
    //gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTestCuboid_level_1", 1);

    //SimulationFileWriter::write("D:/GRIDGENERATION/couplingVF/test/simu/", gridBuilder, FILEFORMAT::ASCII);

    //const uint level = 2;
    //gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, level);
    //gridBuilderlevel->setGrids(gridBuilder->getGrids());


    //gridBuilder->addGrid(14.4921875, 14.4921875, 14.4921875, 16.5078125, 16.5078125, 16.5078125, 0.015625, "cpu", "D3Q27", false, false, false);
    //gridBuilder->addGrid(13.984375, 13.984375, 13.984375, 17.015625, 17.015625, 17.015625, 0.03125, "cpu", "D3Q27", false, false, false);
    //gridBuilder->addGrid(13.46875, 13.46875, 13.46875, 17.53125, 17.53125, 17.53125, 0.0625, "cpu", "D3Q27", false, false, false);
    //gridBuilder->addGrid(12.4375, 12.4375, 12.4375, 18.5625, 18.5625, 18.5625, 0.125, "gpu", "D3Q27", false, false, false);
    //gridBuilder->addGrid(10.375, 10.375, 10.375, 20.625, 20.625, 20.625, 0.25, "gpu", "D3Q27", false, false, false);
    //gridBuilder->addGrid(5.25, 5.25, 5.25, 24.75, 24.75, 24.75, 0.5, "gpu", "D3Q27", false, false, false);
    //gridBuilder->addGrid(0.0, 0.0, 0.0, 30.0, 30.0, 30.0, 1.0, "gpu", "D3Q27", true, true, true);


    //gridBuilder->copyDataFromGpu();

    //gridBuilder->meshGeometry("D:/GRIDGENERATION/STL/circleBinaer.stl", 1);
    //gridBuilder->meshGeometry("D:/GRIDGENERATION/STL/circleBinaer.stl", 0);
    //gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTest_level_1", 1);
    //gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTest_level_0", 0);
    //gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTest_level_2", 2);

    SPtr<Parameter> para = Parameter::make();
    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para);
    //SPtr<GridProvider> gridGenerator = GridProvider::makeGridReader(false, para);

    std::ifstream stream;
    stream.open(configPath.c_str(), std::ios::in);
    if (stream.fail())
        throw "can not open config file!\n";

    UPtr<input::Input> input = input::Input::makeInput(stream, "config");

    setParameters(para, input);

    Simulation sim;
    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    sim.init(para, gridGenerator, fileWriter);
    sim.run();
}


int main( int argc, char* argv[])
{
   MPI_Init(&argc, &argv);
   std::string str, str2; 
   if ( argv != NULL )
   {
      str = static_cast<std::string>(argv[0]);
      if (argc > 1)
      {
         str2 = static_cast<std::string>(argv[1]);
         try
         {
             multipleLevel(str2);
         }
         catch (std::exception e)
         {
             std::cout << e.what() << std::flush;
             //MPI_Abort(MPI_COMM_WORLD, -1);
         }
      }
      else
      {
          std::cout << "Configuration file must be set!: lbmgm <config file>" << std::endl << std::flush;
         //MPI_Abort(MPI_COMM_WORLD, -1);
      }
   }
   /*
   MPE_Init_log() & MPE_Finish_log() are NOT needed when
   liblmpe.a is linked with this program.  In that case,
   MPI_Init() would have called MPE_Init_log() already.
   */
#if defined( MPI_LOGGING )
   MPE_Init_log();
#endif

   

#if defined( MPI_LOGGING )
   if ( argv != NULL )
      MPE_Finish_log( argv[0] );
   if ( str != "" )
      MPE_Finish_log( str.c_str() );
   else
      MPE_Finish_log( "TestLog" );
#endif

   MPI_Finalize();
   return 0;
}

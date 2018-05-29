//#define MPI_LOGGING


#include <mpi.h>
#if defined( MPI_LOGGING )
	#include <mpe.h>
#endif

#include <string>
#include <iostream>

#include "LBM/Simulation.h"

#include "DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "grid/GridBuilder/LevelGridBuilder.h"

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "VirtualFluidsBasics/utilities/input/Input.h"
#include "VirtualFluidsBasics/utilities/StringUtil/StringUtil.h"
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
//#include "DataStructureInitializer/GridReaderFiles/GridReader.h"

#include "grid/BoundaryConditions/Side.h"

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
    para->setMaxLevel(StringUtil::toInt(input->getValue("NOGL")));
    para->setGridX(StringUtil::toVector(input->getValue("GridX")));                           
    para->setGridY(StringUtil::toVector(input->getValue("GridY")));                           
    para->setGridZ(StringUtil::toVector(input->getValue("GridZ")));                  
    para->setDistX(StringUtil::toVector(input->getValue("DistX")));                  
    para->setDistY(StringUtil::toVector(input->getValue("DistY")));                  
    para->setDistZ(StringUtil::toVector(input->getValue("DistZ")));                

    para->setNeedInterface(std::vector<bool>{true, true, true, true, true, true});
}



void multipleLevel(const std::string& configPath)
{
    logging::Logger::setStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::RAYCASTING);

    //auto gridBuilderlevel = LevelGridBuilder::makeShared(Device::CPU, "D3Q27");
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////Conglomerate* conglomerate = new Conglomerate();
    ////conglomerate->add(new Cuboid(10, 10, 10, 30, 30, 30));
    ////conglomerate->subtract(new Sphere(30, 20, 20, 4));
    ////gridBuilder->addGrid(conglomerate, 2);

    ////gridBuilder->addCoarseGrid(0, 0, 0, 40, 40, 40, 1.0);
    ////Object* cuboid = new Cuboid(15, 15, 15, 20, 20, 20);
    ////gridBuilder->addGrid(cuboid, 2);

    ////gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 14, 10, 20, 0.25);
    ////TriangularMesh* triangularMesh = TriangularMesh::make("D:/GRIDGENERATION/STL/quadarBinaer.stl", DiscretizationMethod::POINT_IN_OBJECT);



    //gridBuilder->addGrid(new Cuboid(-10, -8, -8, 50, 22, 22), 1);

    //gridBuilder->addCoarseGrid(-16, -14, -14, 59, 28, 29, 1.0);
    //TriangularMesh* triangularMesh = TriangularMesh::make("D:/GRIDGENERATION/STL/input/local_input/bruecke.stl");

    gridBuilder->addCoarseGrid(-10, -10, -10, 9, 9, 9, 0.25);
    TriangularMesh* triangularMesh = TriangularMesh::make("D:/GRIDGENERATION/STL/cubeBinaer1x1.stl");

    gridBuilder->addGeometry(triangularMesh);

    gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.001, 0.0, 0.0);
    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.001);

    gridBuilder->setPeriodicBoundaryCondition(false, false, true);

    gridBuilder->setNoSlipBoundaryCondition(SideType::MY);

    //gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.001, 0.0, 0.0);
    //gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.001, 0.0, 0.0);
    //gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.001, 0.0, 0.0);
    //gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.001, 0.0, 0.0);
    //gridBuilder->setVelocityBoundaryCondition(SideType::PZ, 0.001, 0.0, 0.0);


    //TriangleOffsetSurfaceGeneration::createOffsetTriangularMesh(triangularMesh, 5);

    //TriangularMesh* sphere = TriangularMesh::make("D:/GRIDGENERATION/STL/GTI.stl", DiscretizationMethod::RAYCASTING);
    //TransformatorImp trans(1.0, Vertex(5.5, 1, 12));
    //trans.transformWorldToGrid(*sphere);
    //STLWriter::writeSTL(sphere->triangleVec, "D:/GRIDGENERATION/STL/GTI2.stl", false);

    //gridBuilder->addGrid(new Sphere(20, 20, 20, 8));
    //gridBuilder->addGeometry(triangularMesh, 1);

    //gridBuilder->addFineGrid(new Cuboid(15, 15, 15, 25, 25, 25), 1);
    //gridBuilder->addFineGrid(new Cuboid(17, 17, 17, 23, 23, 23), 2);


    //gridBuilder->addFineGrid(17.0, 17.0, 17.0, 20.0, 20.0, 20.0, 3);
    //gridBuilder->addFineGrid(10.0, 10.0, 10.0, 20.0, 20.0, 20.0, 3);


    //gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTest_level_2", 2);

    gridBuilder->buildGrids();
    //SimulationFileWriter::write("D:/GRIDGENERATION/files/", gridBuilder, FILEFORMAT::ASCII);

    gridBuilder->writeGridsToVtk("D:/GRIDGENERATION/");
    gridBuilder->writeArrows("D:/arrows");


    SPtr<Parameter> para = Parameter::make();
    SPtr<GridProvider> gridGenerator = GridGenerator::make(gridBuilder, para);
    //SPtr<GridProvider> gridGenerator = GridReader::make(FileFormat::BINARY, para);

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

//#define MPI_LOGGING


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

#include "core/LbmOrGks.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"

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

#include "geometries/Sphere/Sphere.h"
#include "geometries/VerticalCylinder/VerticalCylinder.h"
#include "geometries/Cuboid/Cuboid.h"
#include "geometries/TriangularMesh/TriangularMesh.h"
#include "geometries/Conglomerate/Conglomerate.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"

#include "grid/GridFactory.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"
//#include "grid/GridMocks.h"
//#include "grid/GridStrategy/GridStrategyMocks.h"
#include "VirtualFluidsBasics/utilities/logger/Logger.h"
//#include "VirtualFluidsBasics/basics/utilities/UbLogger.h"
#include "io/STLReaderWriter/STLReader.h"
#include "io/STLReaderWriter/STLWriter.h"
#include "Output/FileWriter.h"
#include "DataStructureInitializer/GridReaderFiles/GridReader.h"

#include "utilities/math/Math.h"

#include "grid/BoundaryConditions/Side.h"
#include "grid/BoundaryConditions/BoundaryCondition.h"

#include "utilities/communication.h"

std::string getGridPath(std::shared_ptr<Parameter> para, std::string Gridpath)
{
    if (para->getNumprocs() == 1)
        return Gridpath + "/";
    
    return Gridpath + "/" + StringUtil::toString(para->getMyID()) + "/";
}

void setParameters(std::shared_ptr<Parameter> para, std::unique_ptr<input::Input> &input)
{
	Communicator* comm = Communicator::getInstanz();

	para->setMaxDev(StringUtil::toInt(input->getValue("NumberOfDevices")));
	para->setNumprocs(comm->getNummberOfProcess());
	para->setDevices(StringUtil::toVector(input->getValue("Devices")));
	para->setMyID(comm->getPID());
	
	std::string _path = input->getValue("Path");
    std::string _prefix = input->getValue("Prefix");
    std::string _gridpath = input->getValue("GridPath");
    std::string gridPath = getGridPath(para, _gridpath);
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
    //std::ofstream logFile( "C:/Users/lenz/Desktop/Work/gridGenerator/grid/gridGeneratorLog.txt" );
    std::ofstream logFile( "grid/gridGeneratorLog.txt" );
    logging::Logger::addStream(&logFile);

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    //UbLog::reportingLevel() = UbLog::logLevelFromString("DEBUG5");

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    //gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::RAYCASTING);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    //gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_UNDER_TRIANGLE);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    
    SPtr<Parameter> para = Parameter::make();
    SPtr<GridProvider> gridGenerator;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool useGridGenerator = false;

    if(useGridGenerator){

        enum testCase{ 
            DrivAer,
            DLC,
            MultiGPU
        };

        int testcase = DrivAer;
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( testcase == DrivAer )
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        {
            real dx = 0.2;
            real vx = 0.05;

            TriangularMesh* DrivAerSTL = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DrivAer_Fastback_Coarse.stl");
            //TriangularMesh* triangularMesh = TriangularMesh::make("M:/TestGridGeneration/STL/DrivAer_NoSTLGroups.stl");
            //TriangularMesh* triangularMesh = TriangularMesh::make("M:/TestGridGeneration/STL/DrivAer_Coarse.stl");
            //TriangularMesh* DrivAerSTL = TriangularMesh::make("stl/DrivAer_Fastback_Coarse.stl");

            TriangularMesh* DrivAerRefBoxSTL = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/DrivAer_REF_BOX_Adrea.stl");
            //TriangularMesh* DrivAerRefBoxSTL = TriangularMesh::make("stl/DrivAer_REF_BOX_Adrea.stl");

            real z0 = 0.318+0.5*dx;

            gridBuilder->addCoarseGrid(- 5.0, -5.0, 0.0 - z0,
                                        15.0,  5.0, 5.0 - z0, dx);  // DrivAer

            //Object* floorBox = new Cuboid( -0.3, -1, -1, 4.0, 1, 0.2 );
            //Object* wakeBox  = new Cuboid(  3.5, -1, -1, 5.5, 1, 0.8 );

            //Conglomerate* refRegion = new Conglomerate();

            //refRegion->add(floorBox);
            //refRegion->add(wakeBox);
            //refRegion->add(DrivAerRefBoxSTL);

            gridBuilder->setNumberOfLayers(10,8);
            gridBuilder->addGrid(DrivAerRefBoxSTL, 2);
        
            //gridBuilder->setNumberOfLayers(10,8);
            //gridBuilder->addGrid(DrivAerSTL, 5);

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
            
            //////////////////////////////////////////////////////////////////////////

            SPtr<Grid> grid = gridBuilder->getGrid(gridBuilder->getNumberOfLevels() - 1);

            gridBuilder->getGeometryBoundaryCondition(gridBuilder->getNumberOfLevels() - 1)->setTangentialVelocityForPatch( grid, 4, 0.0075, -2.0, 0.0,
                                                                                                                                     0.0075,  2.0, 0.0, -vx, 0.318);
            gridBuilder->getGeometryBoundaryCondition(gridBuilder->getNumberOfLevels() - 1)->setTangentialVelocityForPatch( grid, 3, 2.793 , -2.0, 0.0,
                                                                                                                                     2.793 ,  2.0, 0.0, -vx, 0.318);

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->writeGridsToVtk("C:/Users/lenz/Desktop/Work/gridGenerator/grid/DrivAer_Grid");
            gridBuilder->writeArrows    ("C:/Users/lenz/Desktop/Work/gridGenerator/grid/DrivAer_Grid_arrow");

            //SimulationFileWriter::write("D:/GRIDGENERATION/files/", gridBuilder, FILEFORMAT::ASCII);
            //SimulationFileWriter::write("C:/Users/lenz/Desktop/Work/gridGenerator/grid/", gridBuilder, FILEFORMAT::ASCII);
            SimulationFileWriter::write("C:/Users/lenz/Desktop/Work/gridGenerator/grid/", gridBuilder, FILEFORMAT::BINARY);
            //SimulationFileWriter::write("grid/", gridBuilder, FILEFORMAT::ASCII);

            return;

            gridGenerator = GridGenerator::make(gridBuilder, para);
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

            gridGenerator = GridGenerator::make(gridBuilder, para);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( testcase == MultiGPU )
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        {
            //const uint generatePart = 1;
            const uint generatePart = Communicator::getInstanz()->getPID();
            
            std::ofstream logFile2;
            
            if( generatePart == 0 )
                logFile2.open( "C:/Users/lenz/Desktop/Work/gridGenerator/grid/0/gridGeneratorLog.txt" );
            
            if( generatePart == 1 )
                logFile2.open( "C:/Users/lenz/Desktop/Work/gridGenerator/grid/1/gridGeneratorLog.txt" );

            logging::Logger::addStream(&logFile2);

            real dx = 1.0 / 40.0;
            real vx = 0.05;

            TriangularMesh* triangularMesh = TriangularMesh::make("C:/Users/lenz/Desktop/Work/gridGenerator/stl/ShpereNotOptimal.stl");
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
                gridBuilder->findCommunicationIndices(CommunicationDirections::PX);
                gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 1);
            }
            
            if( generatePart == 1 ){
                gridBuilder->findCommunicationIndices(CommunicationDirections::MX);
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
        
            //////////////////////////////////////////////////////////////////////////

            if (generatePart == 0) {
                gridBuilder->writeGridsToVtk("C:/Users/lenz/Desktop/Work/gridGenerator/grid/0/Test_");
                gridBuilder->writeArrows    ("C:/Users/lenz/Desktop/Work/gridGenerator/grid/0/Test_Arrow");
            }
            if (generatePart == 1) {
                gridBuilder->writeGridsToVtk("C:/Users/lenz/Desktop/Work/gridGenerator/grid/1/Test_");
                gridBuilder->writeArrows    ("C:/Users/lenz/Desktop/Work/gridGenerator/grid/1/Test_Arrow");
            }

            if (generatePart == 0)
                SimulationFileWriter::write("C:/Users/lenz/Desktop/Work/gridGenerator/grid/0/", gridBuilder, FILEFORMAT::ASCII);
            if (generatePart == 1)
                SimulationFileWriter::write("C:/Users/lenz/Desktop/Work/gridGenerator/grid/1/", gridBuilder, FILEFORMAT::ASCII);

            //return;

            gridGenerator = GridGenerator::make(gridBuilder, para);
        }
    }
    else
    {
        gridGenerator = GridReader::make(FileFormat::BINARY, para);
        //gridGenerator = GridReader::make(FileFormat::ASCII, para);
    }

    logFile.close();

    //return;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::ifstream stream;
    stream.open(configPath.c_str(), std::ios::in);
    if (stream.fail())
        throw std::runtime_error("can not open config file!");

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
            catch (const std::exception& e)
            {
                *logging::out << logging::Logger::ERROR << e.what() << "\n";
                //MPI_Abort(MPI_COMM_WORLD, -1);
            }
            catch (...)
            {
                std::cout << "unknown exeption" << std::endl;
            }
        }
        else
        {
            try
            {
                multipleLevel("C:/Users/lenz/Desktop/Work/gridGenerator/inp/configTest.txt");
            }
            catch (const std::exception& e)
            {
                
                *logging::out << logging::Logger::ERROR << e.what() << "\n";
                //std::cout << e.what() << std::flush;
                //MPI_Abort(MPI_COMM_WORLD, -1);
            }
            catch (const std::bad_alloc e)
            {
                
                *logging::out << logging::Logger::ERROR << "Bad Alloc:" << e.what() << "\n";
                //std::cout << e.what() << std::flush;
                //MPI_Abort(MPI_COMM_WORLD, -1);
            }
            catch (...)
            {
                *logging::out << logging::Logger::ERROR << "Unknown exception!\n";
                //std::cout << "unknown exeption" << std::endl;
            }

            std::cout << "\nConfiguration file must be set!: lbmgm <config file>" << std::endl << std::flush;
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

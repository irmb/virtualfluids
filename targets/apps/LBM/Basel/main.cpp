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

#include "metis.h"

#include "Core/LbmOrGks.h"
#include "Core/Input/Input.h"
#include "Core/StringUtilities/StringUtil.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"

#include "global.h"

#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"

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
    //std::ofstream logFile( "F:/Work/Computations/gridGenerator/grid/gridGeneratorLog.txt" );
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

            real dx = 1.2;
            real vx = 0.05;

            TriangularMesh* BaselSTL = TriangularMesh::make("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/stl/BaselUrbanProfile_066_deg_bridge_All_CLOSED.stl");

            gridBuilder->addCoarseGrid(-256.0, -256.0, -  8.0,
                                        256.0,  256.0,  160.0, dx);  

            gridBuilder->addGeometry(BaselSTL);

            gridBuilder->setPeriodicBoundaryCondition(true, true, false);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx , 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx , 0.0, 0.0);

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

            //////////////////////////////////////////////////////////////////////////

            gridBuilder->writeGridsToVtk("C:/Users/hiwi/Desktop/Basel_Ergebnisse/grids/BaselUni/Basel_Grid");

            SimulationFileWriter::write("C:/Users/hiwi/Desktop/Basel_Ergebnisse/grids/BaselUni/", gridBuilder, FILEFORMAT::BINARY);

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			StreetPointFinder finder;

			finder.readStreets("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/targets/apps/LBM/streetTest/resources/ExampleStreets.txt");

			finder.writeVTK("C:/Users/hiwi/Desktop/Basel_Ergebnisse/ExampleStreets.vtk");

			finder.findIndicesLB(gridBuilder->getGrid(0));

			finder.writeConnectionVTK("C:/Users/hiwi/Desktop/Basel_Ergebnisse/grids/BaselUni/Basel_Grid/ExampleStreetsConnection.vtk", gridBuilder->getGrid(0));

			finder.writeSimulationFile("C:/Users/hiwi/Desktop/Basel_Ergebnisse/grids/BaselUni/", 1.0, gridBuilder->getNumberOfLevels(), 0);

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			return;

            gridGenerator = GridGenerator::make(gridBuilder, para);

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
    //if (stream.fail())
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
				multipleLevel("C:/Users/hiwi/Desktop/MS/configBasel.txt");
				//multipleLevel("F:/Work/Computations/gridGenerator/inp/configTest.txt");
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

            //std::cout << "\nConfiguration file must be set!: lbmgm <config file>" << std::endl << std::flush;
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

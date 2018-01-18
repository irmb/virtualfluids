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
#include "geometries/Geometry/Geometry.cuh"

#include "grid/GridFactory.h"

std::string getGridPath(std::shared_ptr<Parameter> para, std::string Gridpath)
{
    if (para->getNumprocs() == 1)
        return Gridpath + "/";
    else
        return Gridpath + "/" + StringUtil::toString(para->getMyID()) + "/";
}

void setParameters(std::shared_ptr<Parameter> para, std::unique_ptr<input::Input> &input)
{
    std::string _path = input->getValue("Path");
    std::string _prefix = input->getValue("Prefix");
    std::string _gridpath = input->getValue("GridPath");
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
    para->setMaxLevel(StringUtil::toInt(input->getValue("NOGL")));
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
    para->setGridX(StringUtil::toVector(input->getValue("GridX")));                      //GridX = StringUtil::toVector<int>(input->getValue( "GridX" ));          
    para->setGridY(StringUtil::toVector(input->getValue("GridY")));                      //GridY = StringUtil::toVector<int>(input->getValue( "GridY" ));          
    para->setGridZ(StringUtil::toVector(input->getValue("GridZ")));                      //GridZ = StringUtil::toVector<int>(input->getValue( "GridZ" ));
    para->setDistX(StringUtil::toVector(input->getValue("DistX")));                      //DistX = StringUtil::toVector<int>(input->getValue( "DistX" ));
    para->setDistY(StringUtil::toVector(input->getValue("DistY")));                      //DistY = StringUtil::toVector<int>(input->getValue( "DistY" ));
    para->setDistZ(StringUtil::toVector(input->getValue("DistZ")));                      //DistZ = StringUtil::toVector<int>(input->getValue( "DistZ" )); 
                                                                                          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    para->setNeedInterface(std::vector<bool>{true, true, true, true, true, true});
}



void multipleLevel(const std::string& configPath)
{
    

    SPtr<LevelGridBuilder> gridBuilder(new LevelGridBuilder());
    gridBuilder->addGrid(20, 0, -10, 80, 60, 50, 1.0, "cpu", "D3Q27");
    gridBuilder->addGrid(35.25, 3.25, -4.75, 70.75, 40.75, 25.75, 0.5, "cpu", "D3Q27");


    gridBuilder->meshGeometry("D:/GRIDGENERATION/STL/circleBinaer.stl", 1);
    gridBuilder->meshGeometry("D:/GRIDGENERATION/STL/circleBinaer.stl", 0);
    gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTest3", 1);
    gridBuilder->writeGridToVTK("D:/GRIDGENERATION/gridTest2", 0);

    //SPtr<Parameter> para = Parameter::make();
    //SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(builder, para);

    //std::ifstream stream;
    //stream.open(configPath.c_str(), std::ios::in);
    //if (stream.fail())
    //    throw "can not open config file!\n";

    //UPtr<input::Input> input = input::Input::makeInput(stream, "config");

    //setParameters(para, input);

    //Simulation sim;
    //sim.init(para, gridGenerator);
    //sim.run();
}

void simulate(const std::string& configPath)
{
    SPtr<LevelGridBuilder> builder(new LevelGridBuilder());


    SPtr<Parameter> para = Parameter::make();
    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(builder, para);
    //std::shared_ptr<GridProvider> reader = GridProvider::makeGridReader(true, para);

    std::ifstream stream;
    stream.open(configPath.c_str(), std::ios::in);
    if (stream.fail())
        throw "can not open config file!\n";

    UPtr<input::Input> input = input::Input::makeInput(stream, "config");

    setParameters(para, input);

    SPtr<Transformator> trans(new TransformatorImp());
    //builder->addGrid(para->getGridX()[0], para->getGridY()[0], para->getGridZ()[0], 1.0, "D3Q27", trans);

    SPtr<Transformator> transRefine1(new TransformatorImp(para->getDistX()[1], para->getDistY()[1], para->getDistZ()[1], 0.5));
    //builder->addGrid(para->getGridX()[1], para->getGridY()[1], para->getGridZ()[1], 1.0, "D3Q27", transRefine1);

    //builder->getGridWrapper(0, 0)->copyDataFromGPU();
    //builder->getGridWrapper(1, 0)->copyDataFromGPU();

    //GridVTKWriter::writeSparseGridToVTK(builder->getGridWrapper(0, 0)->grid, "D:/GRIDGENERATION/couplingVF/periodicTaylor/testFile", trans);
    SimulationFileWriter::writeSimulationFiles("D:/GRIDGENERATION/couplingVF/periodicTaylor/simuFiles/", builder, false, trans);

    Simulation sim;
    sim.init(para, gridGenerator);
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
         catch (std::string e)
         {
             std::cout << e << std::flush;
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

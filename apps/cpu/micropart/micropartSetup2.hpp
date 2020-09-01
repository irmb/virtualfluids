#include <iostream>
#include <string>
#include <map>
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "OverlapBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include "geometry3d/GbSystem3D.h"
#include "geometry3d/GbCuboid3D.h"
#include "geometry3d/GbCylinder3D.h"
#include "geometry3d/GbSphere3D.h"
#include "BlocksPostprocessor.h"
#include "Grid3D.h"
#include "Patch3D.h"
#include "Patch3DSystem.h"
#include "Block3D.h"
#include "LBMKernelETD3Q27Cascaded.h"
#include "LBMKernelETD3Q27BGK.h"
#include "CalculationManager.h" 
#include "D3Q27SetConnectorsBlockVisitor.h" 
#include "D3Q27ETInitDistributionsBlockVisitor.h"
#include "D3Q27Interactor.h"
#include "D3Q27NoSlipBCAdapter.h"
#include "D3Q27VelocityBCAdapter.h"
#include "D3Q27DensityBCAdapter.h"
#include "SimulationParameters.h"
#include "Communicator.h"
#include "MPICommunicator.h"
#include "SimpleGeometricPartitioner.h"
#include "D3Q27MacroscopicQuantitiesPostprocessor.h"
#include "D3Q27ETBCProcessor.h"
#include "D3Q27TriFaceMeshInteractor.h"
#include "ConfigFileReader.h"
#include "StringUtil.hpp"
#include "D3Q27PressureDifferencePostprocessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "NUPSCounterPostprocessor.h"
#include "PQueuePartitioningGridVisitor.h"
#include "SetKernelBlockVisitor.h"
#include "GenBlocksGridVisitor.h"
#include "D3Q27PathLinePostprocessor.h"
#include "D3Q27SetUndefinedNodesBlockVisitor.h"
   //
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "geometry3d/creator/GbTriFaceMesh3DCreator.h"
#include "geometry3d/GbTriFaceMesh3D.h"
#include "D3Q27System.h"
#include <basics/transmitter/TbTransmitterMpiPool.h>
#include "MathUtil.hpp"
#include "D3Q27OffsetInterpolationProcessor.h"
#include "SolidBlocksHelper.h"
#include "MetisPartitioningGridVisitor.h"
#include "RestartPostprocessor.h"
#include "D3Q27IncompressibleOffsetInterpolationProcessor.h"
#include "LBMKernelETD3Q27CCLB.h"

using namespace std;

void runSetup2(const char *cstr)
{
   try
   {
      CommunicatorPtr comm(new MPICommunicator());
      int myid = comm->getProcessID();
      int numprocs = comm->getNumberOfProcesses();

      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
      double availMem = 0;
      string geoFile;
      int numOfThreads = 1;

      if(machine == "BOMBADIL") 
      {
         pathname = "c:/temp/micropart";
         availMem = 3.0e9;
         //geoFile = "c:/Data/micropart/DK19_7_02_Martin.stl";
         //geoFile = "c:/Data/micropart/ktoolcav.stl";
         //geoFile = "c:/Data/micropart/boxN.stl";
         //geoFile = "c:/Data/bananas/Banana_boxD.stl";
      }
      else if(machine == "M01" || machine == "M02")      
      {
         pathname = "/work/koskuche/scratch/micropart2";
         availMem = 12.0e9;
         //geoFile = "/home/koskuche/data/micropart/DK19_7_02_Martin.stl";

         numOfThreads = 1;
         if(myid ==0)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      UbLog::reportingLevel() = logINFO;
      //UbLog::reportingLevel() = logDEBUG1;

      int nodePerBlockX1 = 16; //Anzahl an Knoten pro Block
      int nodePerBlockX2 = 16;//(int)16;
      int nodePerBlockX3 = 16;//(int)16;

      double bH = nodePerBlockX1;    //gewuenschte Rand- und Blockbreite

      //Simulation Parameters

      //length [m]
      double lSI = 0.067;
      //length [LB]
      double lLB = 30;

      double dx = 0.0134*0.5;//lSI/lLB;

      double left_offset = 0.5;
      double right_offset  = 0.5;//2*0.5
      double front_offset = 0.15;
      double back_offset  = 0.15;
      double top_offset = 0.0;
      double bottom_offset  = 0.07;

      LBMReal vLB = 0.016103;
      LBMReal Re;
      LBMReal rhoLB = 0.0;
      LBMReal nueLB = 0.0000249*2.0;//(vLB*lLB)/Re;
      Re = (vLB*(0.303/dx))/nueLB;
      const int baseLevel = 0;
      const int refineLevel = 5;
      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      //////////////////////////////////////////////////////////////////////////
      GbObject3DPtr refineCube1(new  GbCuboid3D(-0.2222890,-0.52993, -0.141754, 0.578916113,0.6089970,0.0446053));
      if(myid == 0) GbSystem3D::writeGeoObject(refineCube1.get(), pathname+"/geo/refineCube1", WbWriterVtkXmlASCII::getInstance());

      GbObject3DPtr refineCube2(new  GbCuboid3D(-0.16,-0.05, -0.141754, 0.2,0.05,0.0446053));
      if(myid == 0) GbSystem3D::writeGeoObject(refineCube2.get(), pathname+"/geo/refineCube2", WbWriterVtkXmlASCII::getInstance());
      //////////////////////////////////////////////////////////////////////////

      Grid3DPtr grid(new Grid3D());

      UbSchedulerPtr rSch(new UbScheduler(1000, 1000));
      //RestartPostprocessorPtr rp(new RestartPostprocessor(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY));

      std::string opt;

      if(cstr!= NULL)
         opt = std::string(cstr);

      if/*(cstr== NULL)*/(cstr!= NULL)
      {
         opt = std::string(cstr);

         if(myid==0) UBLOG(logINFO,"Restart step: " << opt);

         //grid = rp->restart(UbSystem::stringTo<int>(opt));
         //rp->reconnect();
         grid->setTimeStep(UbSystem::stringTo<int>(opt));

         if(myid ==0) UBLOG(logINFO,"TimeStep = " <<grid->getTimeStep());

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );

         //domain decomposition
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);
      }
      else
      {
         if(myid ==0)
         {
            UBLOG(logINFO,"L = " <<lLB );
            UBLOG(logINFO,"v = " <<vLB );
            UBLOG(logINFO,"rho = " <<rhoLB );
            UBLOG(logINFO,"nue = " << nueLB );
            UBLOG(logINFO,"Re = " << Re );
            UBLOG(logINFO,"Preprozess - start");
         }


         ////////////////////////////////////////////////////////////////////////
         //Grid
         //////////////////////////////////////////////////////////////////////////
         grid->setDeltaX(dx);
         grid->setBlockNX(nodePerBlockX1, nodePerBlockX2, nodePerBlockX2);

         ////////////////////////////////////////////////////////////////////////////
         //// Geometrie
         ////////////////////////////////////////////////////////////////////////////
         //GbTriFaceMesh3DPtr geo (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(geoFile,"geo"));

         //if(myid == 0) GbSystem3D::writeGeoObject(geo.get(), pathname+"/geo/geo", WbWriterVtkXmlASCII::getInstance());

         ////////////////////////////////////////////////////////////////////////////
         //// Randgeometrien erstellen
         ////////////////////////////////////////////////////////////////////////////
         double shiftForMG=grid->getDeltaX(refineLevel)*nodePerBlockX1 / 3.0*2.0;
         GbCuboid3DPtr plate1  = GbCuboid3DPtr( new GbCuboid3D( -7.5, -1.515e-1, -6.831e-2, 7.5, 1.515e-1, 0.0 ));

         GbCuboid3DPtr plate2  = GbCuboid3DPtr( new GbCuboid3D( -1.5e-1, -16.51e-1, -16.831e-2, 1.5e-1, -1.6e-2, 1.0 ));
         GbCuboid3DPtr plate3  = GbCuboid3DPtr( new GbCuboid3D( -1.5e-1, 1.6e-2, -16.831e-2, 1.5e-1, 16.515e-1, 1.0 ));

         GbCuboid3DPtr plate1_1  = GbCuboid3DPtr( new GbCuboid3D( -7.5, -2.515e-1, -1.0e-1, 7.5, 2.515e-1, -6.831e-2 ));
         GbCuboid3DPtr plate1_2  = GbCuboid3DPtr( new GbCuboid3D( -7.5, -2.515e-1, -0.0000001, 7.5, 2.515e-1, 1.0e-1 ));
         GbCuboid3DPtr plate1_3  = GbCuboid3DPtr( new GbCuboid3D( -7.5, 1.515e-1, -6.831e-2, 7.5, 2.515e-1, 0.0  ));
         GbCuboid3DPtr plate1_4  = GbCuboid3DPtr( new GbCuboid3D( -7.5, -2.515e-1, 0.0, 7.5, -1.515e-1, -1.0e-1 ));

         GbCuboid3DPtr inflow  = GbCuboid3DPtr( new GbCuboid3D( -8.0, -1.0, -1.0, -7.5, 1.0, 1.0 ));
         GbCuboid3DPtr outflow = GbCuboid3DPtr( new GbCuboid3D( 7.5, -1.0, -1.0, 8.0, 1.0, 1.0 ));

         if(myid == 0)
         {
            GbSystem3D::writeGeoObject(plate1.get(),pathname+"/geo/plate1", WbWriterVtkXmlASCII::getInstance());
            GbSystem3D::writeGeoObject(plate2.get(),pathname+"/geo/plate2", WbWriterVtkXmlASCII::getInstance());
            GbSystem3D::writeGeoObject(plate3.get(),pathname+"/geo/plate3", WbWriterVtkXmlASCII::getInstance());

            GbSystem3D::writeGeoObject(plate1_1.get(),pathname+"/geo/plate1_1", WbWriterVtkXmlASCII::getInstance());
            GbSystem3D::writeGeoObject(plate1_2.get(),pathname+"/geo/plate1_2", WbWriterVtkXmlASCII::getInstance());
            GbSystem3D::writeGeoObject(plate1_3.get(),pathname+"/geo/plate1_3", WbWriterVtkXmlASCII::getInstance());
            GbSystem3D::writeGeoObject(plate1_4.get(),pathname+"/geo/plate1_4", WbWriterVtkXmlASCII::getInstance());

            GbSystem3D::writeGeoObject(inflow.get(),pathname+"/geo/inflow", WbWriterVtkXmlASCII::getInstance());
            GbSystem3D::writeGeoObject(outflow.get(),pathname+"/geo/outflow", WbWriterVtkXmlASCII::getInstance());
         }

         GbObject3DPtr gridCube(new GbCuboid3D(plate1->getX1Minimum()-shiftForMG, plate1->getX2Minimum()-shiftForMG, plate1->getX3Minimum()-shiftForMG,
                                                plate1->getX1Maximum()+shiftForMG, 
                                                plate1->getX2Maximum()+shiftForMG, 
                                                plate1->getX3Maximum()+shiftForMG));

         GenBlocksGridVisitor genBlocks;
         genBlocks.addGeoObject(gridCube);
         grid->accept(genBlocks);

         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
            RefineCrossAndInsideGbObjectBlockVisitor refVisitor1(refineCube1, baseLevel, refineLevel-3);
            grid->accept(refVisitor1);

            RefineCrossAndInsideGbObjectBlockVisitor refVisitor2(refineCube2, baseLevel, refineLevel-1);
            grid->accept(refVisitor2);

            RatioBlockVisitor ratioVisitor(refineLevel);
            grid->accept(ratioVisitor);

            RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
            grid->accept(ratioSmoothVisitor);

            OverlapBlockVisitor overlapVisitor(refineLevel);
            grid->accept(overlapVisitor);

            std::vector<int> dirs;
            D3Q27System::getLBMDirections(dirs);
            SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
            grid->accept(interDirsVisitor);
            if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         }

         //////////////////////////////////////////////////////////////////////////
         //INTERAKTOREN SETZEN (=Randbedingungen)
         //////////////////////////////////////////////////////////////////////////
         //oben/unten = Haftrand
         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
         //D3Q27TriFaceMeshInteractorPtr geoInt = D3Q27TriFaceMeshInteractorPtr( new D3Q27TriFaceMeshInteractor(geo, grid, D3Q27BoundaryConditionAdapterPtr(new D3Q27NoSlipBCAdapter(bbOption)),Interactor3D::SOLID));
         //geoInt->setUseHalfSpaceCheck(true);
         //geoInt->setRegardPointInObjectTest(true);

         //D3Q27InteractorPtr plate1Int(new D3Q27Interactor(plate1, grid, bcObst,Interactor3D::INVERSESOLID));
         D3Q27InteractorPtr plate2Int(new D3Q27Interactor(plate2, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr plate3Int(new D3Q27Interactor(plate3, grid, bcObst,Interactor3D::SOLID));

         D3Q27InteractorPtr plate1_1Int(new D3Q27Interactor(plate1_1, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr plate1_2Int(new D3Q27Interactor(plate1_2, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr plate1_3Int(new D3Q27Interactor(plate1_3, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr plate1_4Int(new D3Q27Interactor(plate1_4, grid, bcObst,Interactor3D::SOLID));

         //links: geschwindigkeits-einfluss
         //Velocity-BC
         //////////////////////////////////////////////////////////////////////////
         mu::Parser fct;
         fct.DefineConst("vx1"  , vLB*9.0/4.0 );
         fct = Utilities::getDuctParaboloidX(plate1->getX2Centroid(), plate1->getX2Maximum() - plate1->getX2Minimum(), plate1->getX3Centroid(), plate1->getX3Minimum() - plate1->getX3Maximum(), vLB*9.0/4.0);
         //fct.SetExpr("vx1");
         //////////////////////////////////////////////////////////////////////////

         //////////////////////////////////////////////////////////////////////////
         D3Q27BoundaryConditionAdapterPtr velBCAdapter = D3Q27BoundaryConditionAdapterPtr(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
         velBCAdapter->setSecondaryBcOption(2);
         D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(inflow, grid, velBCAdapter, Interactor3D::SOLID));

         //rechts: druckrand
         //Density-BC
         //fuer Kompressibles Modell  rho = 1.0
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         denBCAdapter->setSecondaryBcOption(1);
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(outflow, grid, denBCAdapter,Interactor3D::SOLID));

         MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B);
         grid->accept( metisVisitor );

         SolidBlocksHelper sd(grid, comm);
         //sd.addInteractor(geoInt);
         sd.addInteractor(inflowInt);
         sd.addInteractor(outflowInt);
         sd.addInteractor(plate1_1Int);
         sd.addInteractor(plate1_2Int);
         sd.addInteractor(plate1_3Int);
         sd.addInteractor(plate1_4Int);
         sd.addInteractor(plate2Int);
         sd.addInteractor(plate3Int);
         sd.deleteSolidBlocks();     

         grid->accept( metisVisitor );

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
         ppblocks->update(0);
         ppblocks.reset();

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nod_temp = nob * (nodePerBlockX1+gl) * (nodePerBlockX2+gl) * (nodePerBlockX3+gl);
         unsigned long nod = nob * (nodePerBlockX1) * (nodePerBlockX2) * (nodePerBlockX3);
         double needMemAll  = double(nod_temp*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());

         if(myid == 0)
         {
            UBLOG(logINFO,"Number of blocks = " << nob);
            UBLOG(logINFO,"Number of nodes  = " << nod);
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
         }  

         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(nodePerBlockX1, nodePerBlockX2, nodePerBlockX2));
         LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(nodePerBlockX1, nodePerBlockX2, nodePerBlockX2));

         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //canal
         //grid->addAndInitInteractor(geoInt);
         grid->addAndInitInteractor(plate1_1Int);
         grid->addAndInitInteractor(plate1_2Int);
         grid->addAndInitInteractor(plate1_3Int);
         grid->addAndInitInteractor(plate1_4Int);
         grid->addAndInitInteractor(plate2Int);
         grid->addAndInitInteractor(plate3Int);

         //inflow
         grid->addAndInitInteractor(inflowInt);

         //outflow
         grid->addAndInitInteractor(outflowInt);

         //////////////////////////////////////////////////////////////////////////
         //connectoren setzen:
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
         //////////////////////////////////////////////////////////////////////////
         //domain decomposition
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);
         //////////////////////////////////////////////////////////////////////////
         //Stroemungsfeld initialisieren
         //////////////////////////////////////////////////////////////////////////
         D3Q27ETInitDistributionsBlockVisitor initVisitor(rhoLB); //1.0
         //initVisitor.setVx1(0.0); 
         grid->accept(initVisitor);

         //if(myid == 0)
         //{
         //   //Abst�nde "q" als Linien rausschreiben
         //   std::vector< UbTupleFloat3 > nodes;
         //   std::vector< UbTupleInt2 >   lines;
         //   geoInt->addQsLineSet(nodes, lines);
         //   WbWriterVtkXmlBinary::getInstance()->writeLines(pathname+"/grid/qs",nodes,lines);
         //}

         if(myid == 0) UBLOG(logINFO,"Preprozess - end");

      }
      //////////////////////////////////////////////////////////////////////////
      //Set Postprozessors
      //////////////////////////////////////////////////////////////////////////
      {
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessor ppgeo(grid,geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv,  comm, true);
         grid->doPostProcess(0);
      }

      UbSchedulerPtr nupsSch(new UbScheduler(1, 5, 10));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      double outTime = 2000;
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid,stepSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv,  comm);
      //////////////////////////////////////////////////////////////////////////
      //PathLine
      //UbSchedulerPtr plSch(new UbScheduler(10, 1500));
      //D3Q27PathLinePostprocessor pathLine(grid, pathname + "/pathLine", WbWriterVtkXmlASCII::getInstance(), conv, plSch, comm, -0.3285474538, 0.09692341,-0.0376166666, nueLB, iProcessor);
      //////////////////////////////////////////////////////////////////////////
      //Simulation
      //////////////////////////////////////////////////////////////////////////
      double endTime = 1000000;
      UbSchedulerPtr visSch(stepSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      if(myid == 0) UBLOG(logINFO,"Simulation-start");
      calculation->calculate();
      if(myid == 0) UBLOG(logINFO,"Simulation-end");

   }
   catch(std::exception& e)
   {
      cerr << e.what() << endl;
   }
   catch(std::string& s)
   {
      cerr << s << endl;
   }
   catch(...)
   {
      cerr << "unknown exception" << endl;
   }

}

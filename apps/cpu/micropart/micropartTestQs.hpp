#include <iostream>
#include <string>
#include <map>
#include <vfluids.h>


using namespace std;

void micropartTestQs(const char *cstr)
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
      int numOfThreads = 3;

      if(machine == "BOMBADIL") 
      {
         pathname = "c:/temp/micropart";
         availMem = 3.0e9;
         //geoFile = "c:/Data/micropart/DK19_7_02_Martin.stl";
         geoFile = "d:/Data/micropart/E0019B_mit_Radien.stl";
         //geoFile = "c:/Data/micropart/boxN.stl";
         //geoFile = "c:/Data/bananas/Banana_boxD.stl";
      }
      else if(machine == "M01" || machine == "M02")      
      {
         pathname = "/work/koskuche/scratch/micropart3";
         //pathname = "/work/koskuche/scratch/micropart2";
         availMem = 12.0e9;
         geoFile = "/home/koskuche/data/micropart/E0019B_mit_Radien_Inv_new_Box.stl";

         numOfThreads = 8;
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

      int nodePerBlockX1 = 8; //Anzahl an Knoten pro Block
      int nodePerBlockX2 = 8;//(int)16;
      int nodePerBlockX3 = 8;//(int)16;

      double bH = nodePerBlockX1;    //gewuenschte Rand- und Blockbreite

      //Simulation Parameters

      //length [m]
      double lSI = 0.067;
      //length [LB]
      double lLB = 30;

      double dx = 5;//0.0134;//lSI/lLB;

      double left_offset = 0.5;
      double right_offset  = 0.5;//2*0.5
      double front_offset = 0.15;
      double back_offset  = 0.15;
      double top_offset = 0.0;
      double bottom_offset  = 0.07;

      LBMReal vLB = 0.016103;
      LBMReal Re;
      LBMReal rhoLB = 0.0;
      LBMReal nueLB = 0.0000249;//(vLB*lLB)/Re;
      Re = (vLB*(0.303/dx))/nueLB;
      const int baseLevel = 0;
      const int refineLevel = 2;
      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      double ft=1000.0;
      //////////////////////////////////////////////////////////////////////////
      GbObject3DPtr refineCube1(new  GbCuboid3D(-0.2222890*ft,-0.52993*ft, -0.141754*ft, /*0.578916113*ft*/275.0,0.6089970*ft,0.0446053*ft));
      if(myid == 0) GbSystem3D::writeGeoObject(refineCube1.get(), pathname+"/geo/refineCube1", WbWriterVtkXmlASCII::getInstance());

      GbObject3DPtr refineCube2(new  GbCuboid3D(-0.16*ft-10.0,-0.05*ft, -0.141754*ft, 0.2*ft+10.0,0.05*ft,0.0446053*ft));
      if(myid == 0) GbSystem3D::writeGeoObject(refineCube2.get(), pathname+"/geo/refineCube2", WbWriterVtkXmlASCII::getInstance());
      //////////////////////////////////////////////////////////////////////////

      Grid3DPtr grid(new Grid3D());

      UbSchedulerPtr rSch(new UbScheduler(1000, 1000));
      RestartPostprocessorPtr rp(new RestartPostprocessor(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY));

      std::string opt;

      if(cstr!= NULL)
         opt = std::string(cstr);

      if/*(cstr== NULL)*/(cstr!= NULL)
      {
         opt = std::string(cstr);

         if(myid==0) UBLOG(logINFO,"Restart step: " << opt);

         grid = rp->restart(UbSystem::stringTo<int>(opt));
         rp->reconnect(grid);
         grid->setTimeStep(UbSystem::stringTo<int>(opt));

         if(myid ==0) UBLOG(logINFO,"TimeStep = " <<grid->getTimeStep());

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
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
         UBLOG(logINFO,"Read geometry: start");
         GbTriFaceMesh3DPtr geo (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(geoFile,"geo"));
         UBLOG(logINFO,"Read geometry: end");
         if(myid == 0) GbSystem3D::writeGeoObject(geo.get(), pathname+"/geo/geo", WbWriterVtkXmlASCII::getInstance());


         ////////////////////////////////////////////////////////////////////////////
         //// Randgeometrien erstellen
         ////////////////////////////////////////////////////////////////////////////
         double shiftForMG=grid->getDeltaX(refineLevel)*nodePerBlockX1 / 3.0*2.0;

         GbCuboid3DPtr inflow  = GbCuboid3DPtr( new GbCuboid3D(geo->getX1Minimum()+9100.0,  geo->getX2Minimum()-200.0, geo->getX3Minimum()-200.0,
                                                               geo->getX1Minimum()+10000.0, geo->getX2Maximum()+200.0, geo->getX3Maximum()+200.0));
         GbCuboid3DPtr outflow  = GbCuboid3DPtr( new GbCuboid3D(geo->getX1Maximum()-10000.0,  geo->getX2Minimum()-200.0, geo->getX3Minimum()-200.0,
                                                                geo->getX1Maximum()-9100.0, geo->getX2Maximum()+200.0, geo->getX3Maximum()+200.0));

         if(myid == 0)
         {
            GbSystem3D::writeGeoObject(inflow.get(),pathname+"/geo/inflow", WbWriterVtkXmlASCII::getInstance());
            GbSystem3D::writeGeoObject(outflow.get(),pathname+"/geo/outflow", WbWriterVtkXmlASCII::getInstance());
         }

         //GbObject3DPtr gridCube(new GbCuboid3D(geo->getX1Minimum()-(double)nodePerBlockX1*dx, geo->getX2Minimum()-(double)nodePerBlockX1*dx, geo->getX3Minimum()-(double)nodePerBlockX1*dx,
         //   geo->getX1Maximum()+(double)nodePerBlockX1*dx, 
         //   geo->getX2Maximum()+(double)nodePerBlockX1*dx, 
         //   geo->getX3Maximum()+(double)nodePerBlockX1*dx));


         shiftForMG=0.0;
         GbObject3DPtr gridCube(new GbCuboid3D(geo->getX1Minimum()+10000.0, geo->getX2Minimum()-shiftForMG, -0.141754*ft/2.0/*geo->getX3Minimum()-shiftForMG*/,
            geo->getX1Maximum()-10000.0, 
            geo->getX2Maximum()+shiftForMG, 
            geo->getX3Maximum()+shiftForMG));

         if(myid == 0) GbSystem3D::writeGeoObject(gridCube.get(),pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());

         GenBlocksGridVisitor genBlocks;
         genBlocks.addGeoObject(gridCube);
         grid->accept(genBlocks);

         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
            //RefineCrossAndInsideGbObjectBlockVisitor refVisitor1(refineCube1, baseLevel, refineLevel-3);
            //grid->accept(refVisitor1);

            RefineCrossAndInsideGbObjectBlockVisitor refVisitor2(refineCube1, refineLevel);
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
         int bbOption = 2; //0=simple Bounce Back, 1=quadr. BB, 2=quadr. BB 2nd choice 
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
         D3Q27TriFaceMeshInteractorPtr geoInt = D3Q27TriFaceMeshInteractorPtr( new D3Q27TriFaceMeshInteractor(geo, grid, D3Q27BoundaryConditionAdapterPtr(new D3Q27NoSlipBCAdapter(bbOption)),Interactor3D::INVERSESOLID));
         geoInt->setUseHalfSpaceCheck(true);
         geoInt->setRegardPointInObjectTest(true);

         //links: geschwindigkeits-einfluss
         //Velocity-BC
         //////////////////////////////////////////////////////////////////////////
         mu::Parser fct;
         fct.DefineConst("vx1"  , vLB           );
         //fct = MathUtil::getDuctParaboloidX(plate1->getX2Centroid(), plate1->getX2Maximum() - plate1->getX2Minimum(), plate1->getX3Centroid(), plate1->getX3Minimum() - plate1->getX3Maximum(), vLB*9.0/4.0);
         fct.SetExpr("vx1");
         //////////////////////////////////////////////////////////////////////////
         //////////////////////////////////////////////////////////////////////////
         D3Q27BoundaryConditionAdapterPtr velBCAdapter = D3Q27BoundaryConditionAdapterPtr(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
         velBCAdapter->setSecondaryBcOption(2);
         D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(inflow, grid, velBCAdapter, Interactor3D::SOLID));
         //D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(inflow, grid, bcObst, Interactor3D::SOLID));

         //rechts: druckrand
         //Density-BC
         //fuer Kompressibles Modell  rho = 1.0
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(outflow, grid, D3Q27BoundaryConditionAdapterPtr(new D3Q27DensityBCAdapter(rhoLB)),Interactor3D::SOLID));
         //D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(outflow, grid, bcObst,Interactor3D::SOLID));

         MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B);
         grid->accept( metisVisitor );

         SolidBlocksHelper sd(grid, comm);
         sd.addInteractor(geoInt);
         sd.addInteractor(inflowInt);
         sd.addInteractor(outflowInt);
         sd.deleteSolidBlocks();     

         grid->accept( metisVisitor );

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
         ppblocks->update(0);
         ppblocks.reset();

         UBLOG(logINFO,grid->getBlock(10,12,0,1)->toString());
         vector<Block3DPtr> blocks;
         //grid->getNeighborBlocksForDirection(D3Q27System::W,10,12,0,1,3,blocks);
         grid->getNeighborBlocksForDirection(D3Q27System::E,4,6,0,0,2,blocks);
         BOOST_FOREACH(Block3DPtr b, blocks)
            UBLOG(logINFO, b->toString());



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
         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27BGK(nodePerBlockX1, nodePerBlockX2, nodePerBlockX2, true));
         LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(nodePerBlockX1, nodePerBlockX2, nodePerBlockX2,0));

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
         grid->addAndInitInteractor(geoInt);

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

         if(myid == 0)
         {
            //Abstände "q" als Linien rausschreiben
            std::vector< UbTupleFloat3 > nodes;
            std::vector< UbTupleInt2 >   lines;
            geoInt->addQsLineSet(nodes, lines);
            WbWriterVtkXmlBinary::getInstance()->writeLines(pathname+"/grid/qs",nodes,lines);
         }

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

      double outTime = 100;
      UbSchedulerPtr visSch(new UbScheduler(outTime));
      //visSch->addSchedule(20, 1010, 1100);
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid,visSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv,  comm);
      //////////////////////////////////////////////////////////////////////////
      //PathLine
      //UbSchedulerPtr plSch(new UbScheduler(10, 1500));
      //D3Q27PathLinePostprocessor pathLine(grid, pathname + "/pathLine", WbWriterVtkXmlASCII::getInstance(), conv, plSch, comm, -0.3285474538, 0.09692341,-0.0376166666, nueLB, iProcessor);
      //////////////////////////////////////////////////////////////////////////
      //Simulation
      //////////////////////////////////////////////////////////////////////////
      double endTime = 1000;
      UbSchedulerPtr visSch1(new UbScheduler(1));
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




#include <iostream>
#include <string>
#include <math.h> 

#include <sys/types.h> //mkdir rights
#include <sys/stat.h> //mkdir
#include <vfluids.h>

using namespace std;


void run(const char *cstr)
{
   try
   {
      string pathname; 
      string pathnameRestart;
      int numOfThreads =1;
      bool logfile = false;
      stringstream logFilename;
      double availMem = 0;

      UbLog::reportingLevel() = logINFO;

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      string machine = string(cstr);

      if(machine == "my") 
      {
         pathname = "d:/temp/BKanal";
         numOfThreads = 4;
         logfile = false;
         availMem = 10.0e9;
      }
      else if(machine == "Ludwig")      
      {
         pathname =        "/work/koskuche/SFB880/BKanal";
         pathnameRestart = "/work/koskuche/SFB880/BKanal";
         numOfThreads = 8;
         availMem = 1.0e9;
         logfile = true;

         if (myid==0) 
         {
            const char* str = pathname.c_str();
#if defined(__unix__)
            int status=mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif 
         }

         if(myid ==0)
         {
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(myid)+".txt";
         }

      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      if(myid ==0 && logfile)
      {
         UbLog::output_policy::setStream(logFilename.str());
      }

      UBLOG(logINFO,"Testcase BreugemChannel");

      //////////////////////////////////////////////////////////////////////////
      //physik
      //////////////////////////////////////////////////////////////////////////
      double Re    = 5500;
      double uLB   = 0.1;  
      double rhoLB = 0.0;

      int blocknx[3];
      blocknx[0] = 20;//10;//5;
      blocknx[1] = 20;//10;//5;
      blocknx[2] = 20;//10;//5;

      int nx[3];
      nx[0] = 15;
      nx[1] = 10;
      nx[2] = 10;

      double coarseNodeDx = 1.0;
      double H     = (double)(nx[2]*blocknx[2])/2.0;
      double hLB   = H/coarseNodeDx;
      double nuLB  = (uLB*hLB)/Re;

      int baseLevel   = 0;
      int refineLevel = 0;//2;//1;////3;//3 soll 1 test;

      ///////////////Weltabmessungen:
      double blockLengthx1 = blocknx[0]*coarseNodeDx; //geowerte
      double blockLengthx2 = blockLengthx1;
      double blockLengthx3 = blockLengthx1;

      double originX1 = 0.0;
      double originX2 = 0.0;
      double originX3 = 0.0;

      //bounding box
      double g_minX1 = originX1;
      double g_minX2 = originX2;
      double g_minX3 = originX3;

      double g_maxX1 = originX1 + 3.0*H;
      double g_maxX2 = originX2 + 2.0*H;
      double g_maxX3 = originX3 + 2.0*H;

      //double geoLength[]   = {  nx[0]*blockLengthx1, nx[1]*blockLengthx2, nx[2]*blockLengthx3}; 

      bool periodicx1 = true;
      bool periodicx2 = true;
      bool periodicx3 = false;

      LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter());

      Grid3DPtr grid(new Grid3D(comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(10000,10000,10000000));
      RestartPostprocessor rp(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY);
      grid = rp.restart(-1);
      //////////////////////////////////////////////////////////////////////////

       if (grid->getTimeStep() == 0)
       {

         //set grid
         grid->setDeltaX(coarseNodeDx);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
         grid->setPeriodicX1(periodicx1);
         grid->setPeriodicX2(periodicx2);
         grid->setPeriodicX3(periodicx3);


         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if(myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);



         //bottom and top solid bc
         //iteractors
         int bbOption1 = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption1));
         double geoOverlap = coarseNodeDx;
         GbCuboid3DPtr bottomBCCuboid(new GbCuboid3D(g_minX1-blockLengthx1, g_minX2-blockLengthx1, g_minX3-blockLengthx1, g_maxX1+blockLengthx1, g_maxX2+blockLengthx1, g_minX3));
         if(myid == 0) GbSystem3D::writeGeoObject(bottomBCCuboid.get(), pathname+"/geo/bottomBCCuboid", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr bottomBCInteractor(new D3Q27Interactor(bottomBCCuboid,grid,bcObst,Interactor3D::SOLID)); 

         GbCuboid3DPtr topBCCuboid(new GbCuboid3D(g_minX1-blockLengthx1, g_minX2-blockLengthx1, g_maxX3, g_maxX1+blockLengthx1, g_maxX2+blockLengthx1, g_maxX3+blockLengthx1));
         if(myid == 0) GbSystem3D::writeGeoObject(topBCCuboid.get(), pathname+"/geo/topBCCuboid", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr topBCInteractor(new D3Q27Interactor(topBCCuboid,grid,bcObst,Interactor3D::SOLID)); 
         //grid->addAndInitInteractor( bottomBCInteractor ); 
         // grid->addAndInitInteractor( topBCInteractor ); 
         //////////////////////////////////////////////////////////////////////////
         if(myid == 0)
         {
            UBLOG(logINFO, "*****************************************");
            UBLOG(logINFO, "* Parameters                            *");
            UBLOG(logINFO, "* Re            ="<<Re);
            //UBLOG(logINFO, "* Ma            ="<<Ma);
            //UBLOG(logINFO, "* uReal         ="<<uReal);
            //UBLOG(logINFO, "* nueReal       ="<<nueReal);
            UBLOG(logINFO, "* nue           ="<<nuLB);
            UBLOG(logINFO, "* velocity      ="<<uLB);
            // UBLOG(logINFO, "* LX1 (world/LB)="<<kanallaengeSI<<"/"<<kanallaengeSI/coarseNodeDx);
            //  UBLOG(logINFO, "* LX2 (world/LB)="<<kanalbreiteSI<<"/"<<kanalbreiteSI/coarseNodeDx);
            //UBLOG(logINFO, "* LX3 (world/LB)="<<kanalhoeheSI<<"/"<<kanalhoeheSI/coarseNodeDx);
            UBLOG(logINFO, "* cdx           ="<<coarseNodeDx);
            //UBLOG(logINFO, "* fdx           ="<<fineNodeDx);
            UBLOG(logINFO, "* dx_base       ="<<coarseNodeDx<<" == "<<coarseNodeDx);
            //UBLOG(logINFO, "* dx_refine     ="<<fineNodeDx<<" == "<<fineNodeDx );
            //UBLOG(logINFO, "* nx1/2/3       ="<<nx[0]<<"/"<<nx[1]<<"/"<<nx[2]);
            UBLOG(logINFO, "* blocknx1/2/3  ="<<blocknx[0]<<"/"<<blocknx[1]<<"/"<<blocknx[2]);
            UBLOG(logINFO, "* x2Periodic    ="<<periodicx2);
            UBLOG(logINFO, "* x3Periodic    ="<<periodicx3);
            UBLOG(logINFO, "* number of threads ="<<numOfThreads);
            UBLOG(logINFO, "* number of processes ="<<comm->getNumberOfProcesses());
            UBLOG(logINFO, "*****************************************");
/*            UBLOGML(logINFO, "UnitConverter:"<<unitConverter->toString());
            UBLOG(logINFO, "*****************************************");  */   
         }

         if(myid == 0) UBLOG(logINFO,"Refinement - start");	

         //////////////////////////////////////////////////////////////////////////
         // refine
         //////////////////////////////////////////////////////////////////////////
         //GbCuboid3DPtr wallsX1X2minRef3(new GbCuboid3D(  originX1-3.0*geoOverlap   , originX2-3.0*geoOverlap  , originX3-3.0*geoOverlap
         //   , originX1+geoLength[0]+geoOverlap, originX2+geoOverlap+geoLength[1], kanalhoeheSI*0.6/*0.55*/));

         //GbCuboid3DPtr wallsX1X2minRef4(new GbCuboid3D(  originX1-3.0*geoOverlap   , originX2-3.0*geoOverlap  ,kanalhoeheSI*0.49
         //   , originX1+geoLength[0]+geoOverlap, originX2+geoOverlap+geoLength[1], kanalhoeheSI*0.53));

         //GbCuboid3DPtr wallsX1X2maxRef2(new GbCuboid3D(  originX1-3.0*geoOverlap   , originX2-3.0*geoOverlap  ,kanalhoeheSI*0.9
         //   , originX1+geoLength[0]+geoOverlap, originX2+geoOverlap+geoLength[1], originX3+geoOverlap+geoLength[2]));

         //GbCuboid3DPtr wallsX1X2maxRef1(new GbCuboid3D(  originX1-3.0*geoOverlap   , originX2-3.0*geoOverlap  ,kanalhoeheSI*0.95
         //   , originX1+geoLength[0]+geoOverlap, originX2+geoOverlap+geoLength[1], originX3+geoOverlap+geoLength[2]));

         //if (refineLevel > 0)
         //{

         //   RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
         //   refineHelper.addGbObject(wallsX1X2minRef3, refineLevel-1);
         //   refineHelper.addGbObject(wallsX1X2minRef4, refineLevel);
         //   refineHelper.addGbObject(wallsX1X2maxRef2, refineLevel-1);
         //   refineHelper.addGbObject(wallsX1X2maxRef1, refineLevel);

         //   refineHelper.refine();
         //   if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         //}

         ///interactoren
         //int bbOption1 = 0; //0=simple Bounce Back, 1=quadr. BB
         //D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption1));
         ///////würfel unten version ende
         ////////////////////////////////////////////////////////////////////////////////
         ////////PM grid
         //Temporär:
         //double  H=1.0;

         vector<D3Q27InteractorPtr> D3Q27InteractorPtrarray;
         ////////////////////////////////////////////////////////////////////////////////
         double dpCubes=(double)H/20.0;
         double distanceXY=dpCubes/2.0-coarseNodeDx*0.5;
         double distanceZ=0;
   
         for (int x = 0; x<30; x++)
            for (int y = 0; y<20; y++)
               for (int z = 0; z<9; z++)
               {
                  double xminCubes = originX1+distanceXY+2.0*dpCubes*x;
                  double yminCubes = originX2+distanceXY+2.0*dpCubes*y;
                  double zminCubes = originX3+distanceZ+2.0*dpCubes*z;
                  double xmaxCubes = xminCubes+dpCubes;
                  double ymaxCubes = yminCubes+dpCubes;
                  double zmaxCubes = zminCubes+dpCubes;
                  GbCuboid3DPtr rectTemp(new GbCuboid3D(xminCubes, yminCubes, zminCubes, xmaxCubes, ymaxCubes, zmaxCubes));
                  D3Q27BoundaryConditionAdapterPtr cubeBCAdapter(new D3Q27NoSlipBCAdapter());                   
                  D3Q27InteractorPtr cubeInteractor( new D3Q27Interactor(rectTemp,grid,cubeBCAdapter,Interactor3D::SOLID));
                  D3Q27InteractorPtrarray.push_back(cubeInteractor); 
               }

         ////////////////
         //ende cubes
         //////////
         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));

         ////////////////////////////////////////////
         /////delete solid blocks
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(topBCInteractor);
         intHelper.addInteractor(bottomBCInteractor);
         for(size_t i=0; i<D3Q27InteractorPtrarray.size(); ++i)
         {
            intHelper.addInteractor(D3Q27InteractorPtrarray[i]);
         }
         intHelper.selectBlocks();
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	 
         //////////////////////////////////////

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         if(myid == 0) UBLOG(logINFO,"Write blocks - start");
         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
         if(myid == 0) ppblocks->update(0);
         if(myid == 0) UBLOG(logINFO,"Write blocks - end");

         unsigned long nob = grid->getNumberOfBlocks();
         unsigned long nod = nob * blocknx[0]*blocknx[1]*blocknx[2];
         unsigned long nod_real = nob * (blocknx[0]+3)*(blocknx[1]+3)*(blocknx[2]+3);

         double needMemAll  = double(nod_real*(27*sizeof(double) + sizeof(int)));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());

         if(myid == 0)
         {
            UBLOG(logINFO,"Number of blocks = " << nob);
            UBLOG(logINFO,"Number of nodes  = " << nod);
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
         }

         LBMKernel3DPtr kernel;
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));
         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);
     
         mu::Parser fctForcingX1;
         fctForcingX1.SetExpr("Fx1*dx");
         fctForcingX1.DefineConst("Fx1", 0.6*5.0e-6);//9.99685e-7);

         kernel->setForcingX1(fctForcingX1);
         kernel->setWithForcing(true); 

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         //////////////////////////////////
         //undef nodes
         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }
         //////////////////////////////////////////
         intHelper.setBC();

         for(size_t i=0; i<D3Q27InteractorPtrarray.size(); ++i)
         {
            grid->addAndInitInteractor( D3Q27InteractorPtrarray[i] ); 
            char numstr[21];
            sprintf(numstr, "%f", (double)i);
            std::string pathObstCube = pathname+"/geo/obstBCCuboid"+ numstr;
            if(myid == 0) GbSystem3D::writeGeoObject(D3Q27InteractorPtrarray[i]->getGbObject3D().get(),
               /* rectTemp.get(),*/ pathObstCube, WbWriterVtkXmlASCII::getInstance());
         }


         ppblocks.reset();

         //inflow
         mu::Parser inflowProfile;
         inflowProfile.SetExpr("uLB*0.9"); 
         inflowProfile.DefineConst("uLB",uLB);

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept( setConnsVisitor );

         //initialization of decompositions
         D3Q27ETInitDistributionsBlockVisitor initVisitor( nuLB,rhoLB);
         initVisitor.setVx1(inflowProfile);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), 
            unitConverter,true));

         ppgeo->update(0);
         ppgeo.reset();
         geoSch.reset();

         if(myid == 0) UBLOG(logINFO,"Preprozess - end");      

      }
            else
            {
               //set forcing
               mu::Parser fctForcingX1, fctForcingX2, fctForcingX3;
               fctForcingX1.SetExpr("Fx1*dx");
               fctForcingX1.DefineConst("Fx1", 0.6*5.0e-6);
               fctForcingX2.SetExpr("0.0");
               fctForcingX3.SetExpr("0.0");
               SetForcingBlockVisitor forcingVisitor(fctForcingX1, fctForcingX2, fctForcingX3);
               grid->accept(forcingVisitor);

               //set connectors
               D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
               D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
               grid->accept( setConnsVisitor );
               if(myid == 0) UBLOG(logINFO,"Restart - end"); 
      }


      UbSchedulerPtr visSch(new UbScheduler());
      visSch->addSchedule(100,100,1000);
      visSch->addSchedule(1000,1000,10000);
      visSch->addSchedule(10000,10000,100000);
      //visSch->addSchedule(20000,20000,800000);
      //visSch->addSchedule(50,350000,350500);
      //visSch->addSchedule(50,420000,420500);
      //visSch->addSchedule(50000,420500,10000000);
      //visSch->addSchedule(2250,140000,450001);
      //UbSchedulerPtr resSch(new UbScheduler());
      //resSch->addSchedule(20000,20,10000000);
      //UbSchedulerPtr resSchRMS(new UbScheduler());
      //resSchRMS->addSchedule(40000,420000,10000000);
      //UbSchedulerPtr resSchMeans(new UbScheduler());
      //resSchMeans->addSchedule(40000,0,10000000);
      //UbSchedulerPtr stepAvSch(new UbScheduler());
      //stepAvSch->addSchedule(20,0,10000000);
      //AverageValuesPostprocessor Avpp(grid, pathname + "/steps/stepAV", WbWriterVtkXmlBinary::getInstance(), 
      //                                visSch/*wann wird rausgeschrieben*/, stepAvSch/*wann wird gemittelt*/, resSchMeans, resSchRMS/*wann wird resettet*/);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), unitConverter);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);


      UbSchedulerPtr AdjForcSch(new UbScheduler());
      AdjForcSch->addSchedule(100,100,20000000);
      //D3Q27IntegrateValuesHelperPtr IntValHelp(new D3Q27IntegrateValuesHelper(grid, comm, 
      //   originX1, originX2, kanalhoeheSI*0.55/*0.501*/, 
      //   nx[0]*blockLengthx1, nx[1]*blockLengthx2, kanalhoeheSI*0.999));
      D3Q27IntegrateValuesHelperPtr IntValHelp(new D3Q27IntegrateValuesHelper(grid, comm, 
         originX1, originX2, g_maxX3*0.55/*0.501*/, 
         g_maxX1, g_maxX2, g_maxX3*0.999));

      double vxZiel=uLB;
      //D3Q27AdjustForcingPostprocessor AdjForcPPPtr(grid, AdjForcSch,IntValHelp, vxZiel*0.6, comm);//da am 11.3.2013 velo um 1/0.6 zu hoch
      D3Q27AdjustForcingPostprocessor AdjForcPPPtr(grid, AdjForcSch,IntValHelp, vxZiel, comm);//dies sollte zu Re=5500 fuehren..

      UbSchedulerPtr visQSch(new UbScheduler());
      visQSch->addSchedule(10,90100,90130);
      QCriterionPostprocessor QKritPtr(grid,pathname+"/steps/Qq",WbWriterVtkXmlBinary::getInstance(),visQSch, comm);

      //mu::Parser decrViscFunc;
      //decrViscFunc.SetExpr("nue0+c0/(t+1)/(t+1)");
      //decrViscFunc.DefineConst("nue0", nuLB);
      //decrViscFunc.DefineConst("c0", 0.1);
      //UbSchedulerPtr DecrViscSch(new UbScheduler());
      //DecrViscSch->addSchedule(10,10,1000);
      //DecreaseViscosityPostprocessor decrViscPPPtr(grid, DecrViscSch,&decrViscFunc, comm);

      cout << "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()<<endl;
      cout << "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed()<<endl;
      cout << "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()<<endl;

      double endTime = 2000000;//20000001;
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      if(myid == 0) UBLOG(logINFO,"Simulation-start");
      calculation->calculate();
      if(myid == 0) UBLOG(logINFO,"Simulation-end");
   }
   catch(std::exception& e)
   {
      cerr << e.what() << endl << flush;
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
int main(int argc, char* argv[])
{

   run(argv[1]);

   return 0;
}


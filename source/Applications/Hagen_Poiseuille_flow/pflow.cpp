#include <iostream>
#include <string>

#include <vfluids.h>

using namespace std;


void pflowForcing(const char *cstr)
{
   try
   {

      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
      int numOfThreads = 4;
      double availMem = 0;

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if(machine == "BOMBADIL") 
      {
         pathname = "d:/temp/Hagen_Poiseuille_flow_forcing2";
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         pathname = "/work/koskuche/scratch/Hagen_Poiseuille_flow_forcing";
         availMem = 12.0e9;

//         if(myid ==0)
//         {
//            stringstream logFilename;
//            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
//            UbLog::output_policy::setStream(logFilename.str());
//         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      double dx = 1;

      //const int blocknx1 = 16;
      //const int blocknx2 = 16;
      //const int blocknx3 = 16;

      const int blocknx1 = 8;
      const int blocknx2 = 8;
      const int blocknx3 = 8;

      LBMReal uLB = 0.05;
      LBMReal rhoLB = 0.0;
      LBMReal nueLB = 0.05842;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = 1;

      //bounding box
      double g_minX1 = 0;
      double g_minX2 = 0;
      double g_minX3 = 0;

      //double g_maxX1 = 16;
      //double g_maxX2 = 128;
      //double g_maxX3 = 16;

      double g_maxX1 = 32;
      double g_maxX2 = 32;
      double g_maxX3 = 32;

      double blockLength = blocknx1*dx;

      Grid3DPtr grid(new Grid3D(comm));
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(true);
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);

      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());      

      //GbObject3DPtr refineCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2/3.0, g_maxX3));
      //if(myid ==0) GbSystem3D::writeGeoObject(refineCube.get(),pathname + "/geo/refineCube", WbWriterVtkXmlBinary::getInstance());

      GbObject3DPtr refineCube(new GbCuboid3D(g_minX1+blockLength+2*dx, g_minX2+blockLength+2*dx, g_minX3+blockLength+2*dx, g_maxX1-blockLength-2*dx, g_maxX2-blockLength-2*dx, g_maxX3-blockLength-2*dx));
      if(myid ==0) GbSystem3D::writeGeoObject(refineCube.get(),pathname + "/geo/refineCube", WbWriterVtkXmlBinary::getInstance());
      
      //UbSchedulerPtr rSch(new UbScheduler());
      //rSch->addSchedule(50, 50, 50);
      //RestartPostprocessorPtr rp(new RestartPostprocessor(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY));

      std::string opt;

      if(cstr!= NULL)
         opt = std::string(cstr);

      if/*(cstr== NULL)*/(cstr!= NULL)
      {
         //if(myid==0) UBLOG(logINFO,"Restart step: " << opt);
         //grid = rp->restart(UbSystem::stringTo<int>(opt));
         //rp->reconnect(grid);

         mu::Parser fctForcingX1;
         mu::Parser fctForcingX2;
         mu::Parser fctForcingX3;
         fctForcingX1.SetExpr("Fx1*dx");
         fctForcingX1.DefineConst("Fx1", 9.99685e-7);
         fctForcingX2.SetExpr("0.0");
         fctForcingX3.SetExpr("0.0");

         SetForcingBlockVisitor forcingVisitor(fctForcingX1, fctForcingX2, fctForcingX3);
         grid->accept(forcingVisitor);

         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
      }
      else
      {

      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      if(myid ==0)
      {
         UBLOG(logINFO,"Parameters:");
         //UBLOG(logINFO,"L = " << L2/dx );
         UBLOG(logINFO,"v = " << uLB );
         UBLOG(logINFO,"rho = " << rhoLB );
         UBLOG(logINFO,"nue = " << nueLB );
         //UBLOG(logINFO,"Re = " << Re );
         UBLOG(logINFO,"dx = " << dx );
         UBLOG(logINFO,"number of levels = " << refineLevel+1 );
         UBLOG(logINFO,"numOfThreads = " << numOfThreads );
         UBLOG(logINFO,"Preprozess - start");
      }

      //walls
      GbCuboid3DPtr addWallYmin (new GbCuboid3D(g_minX1-1.0*blockLength, g_minX2-1.0*blockLength, g_minX3-1.0*blockLength, g_maxX1+1.0*blockLength, g_minX2, g_maxX3+1.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallYmax (new GbCuboid3D(g_minX1-1.0*blockLength, g_maxX2, g_minX3-1.0*blockLength, g_maxX1+1.0*blockLength, g_maxX2+1.0*blockLength, g_maxX3+1.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));

      if (refineLevel > 0)
      {
         if(myid == 0) UBLOG(logINFO,"Refinement - start");	
         RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
         refineHelper.addGbObject(refineCube, 1);
         refineHelper.refine();
         if(myid == 0) UBLOG(logINFO,"Refinement - end");	
      }

      int bbOption = 0; //0=simple Bounce Back, 1=quadr. BB
      D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
      D3Q27BoundaryConditionAdapterPtr bcObst2(new D3Q27SlipBCAdapter(bbOption));

      //walls
      D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcObst,Interactor3D::SOLID));
      D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcObst,Interactor3D::SOLID));

      ////////////////////////////////////////////
      //METIS
      Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));   
      ////////////////////////////////////////////
      /////delete solid blocks
      if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
      InteractorsHelper intHelper(grid, metisVisitor);
      intHelper.addInteractor(addWallYminInt);
      intHelper.addInteractor(addWallYmaxInt);
      intHelper.selectBlocks();
      if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	 
      //////////////////////////////////////

      //set connectors
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
      grid->accept( setConnsVisitor );

      //domain decomposition for threads
      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      grid->accept(pqPartVisitor);

      ppblocks->update(0);
      ppblocks.reset();

      unsigned long nob = grid->getNumberOfBlocks();
      int gl = 3;
      unsigned long nodb = (blocknx1) * (blocknx2) * (blocknx3);
      unsigned long nod = nob * (blocknx1) * (blocknx2) * (blocknx3);
      unsigned long nodg = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);
      double needMemAll  = double(nodg*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
      double needMem  = needMemAll / double(comm->getNumberOfProcesses());

      if(myid == 0)
      {
         UBLOG(logINFO,"Number of blocks = " << nob);
         UBLOG(logINFO,"Number of nodes  = " << nod);
         int minInitLevel = grid->getCoarsestInitializedLevel();
         int maxInitLevel = grid->getFinestInitializedLevel();
         for(int level = minInitLevel; level<=maxInitLevel; level++)
         {
            int nobl = grid->getNumberOfBlocks(level);
            UBLOG(logINFO,"Number of blocks for level " << level <<" = " << nob);
            UBLOG(logINFO,"Number of nodes for level " << level <<" = " << nob*nodb);
         }
         UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
         UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
         UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
      }            

      int kernelType = 2;
      LBMKernel3DPtr kernel;
      if (kernelType == 0)
      {
         rhoLB = 1.0;
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));
      }
      else if (kernelType == 1)
      {
         rhoLB = 1.0;
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3));
      }
      else if (kernelType == 2)
      {
         rhoLB = 0.0;
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));
         //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLBex2(blocknx1, blocknx2, blocknx3, 0, grid));
      }

      mu::Parser fctForcingX1;
      fctForcingX1.SetExpr("Fx1*dx");
      fctForcingX1.DefineConst("Fx1", 9.99685e-7);

      kernel->setForcingX1(fctForcingX1);
      kernel->setWithForcing(true);
      
      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
      grid->accept(kernelVisitor);

      if (refineLevel > 0)
      {
         D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
         grid->accept(undefNodesVisitor);
      }

      //walls
      intHelper.setBC();

      //initialization of distributions
      D3Q27ETInitDistributionsBlockVisitor initVisitor(nueLB, rhoLB);
      grid->accept(initVisitor);

      //Postrozess
      UbSchedulerPtr geoSch(new UbScheduler(1));
      D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
         new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, true));
      ppgeo->update(0);
      ppgeo.reset();

      if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 
}
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      double outTime = 1;
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      //UbSchedulerPtr stepSch(new UbScheduler());
      //stepSch->addSchedule(10, 100, 1000);
      //nodeSch->addSchedule(1000, 1000, 10000);
      //nodeSch->addSchedule(10000, 10000, 50000);
      //stepSch->addSchedule(100, 100, 1000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv);
      
      double fdx = grid->getDeltaX(grid->getFinestInitializedLevel());

      D3Q27IntegrateValuesHelperPtr h1(new D3Q27IntegrateValuesHelper(grid, comm, 
                                       g_minX1, g_minX2, g_minX3, 
                                       g_minX1+1.0*fdx, g_maxX2, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(h1->getBoundingBox().get(),pathname + "/geo/iv1", WbWriterVtkXmlBinary::getInstance());
      D3Q27IntegrateValuesHelperPtr h2(new D3Q27IntegrateValuesHelper(grid, comm, 
                                       g_maxX1-1.0*fdx, g_minX2, g_minX3, 
                                       g_maxX1, g_maxX2, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(h2->getBoundingBox().get(),pathname + "/geo/iv2", WbWriterVtkXmlBinary::getInstance());
      LBMReal rhoReal = rhoLB;
      LBMReal uReal = uLB;
      D3Q27PressureDifferencePostprocessor rhopp(grid, stepSch, pathname + "/results/rho_diff.txt", h1, h2, rhoReal, uReal, uLB, comm);

      double endTime =10000;//10001.0;

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
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
//////////////////////////////////////////////////////////////////////////
void pflowdp(const char *cstr)
{
   try
   {

      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
      int numOfThreads = 1;
      double availMem = 0;

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if(machine == "BOMBADIL") 
      {
         pathname = "d:/temp/Hagen_Poiseuille_flow_tube_testC";
         availMem = 4.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         pathname = "/work/koskuche/scratch/Hagen_Poiseuille_flow_tube5";
         availMem = 10.0e9;
         numOfThreads = 8;

         if(myid ==0)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      double dx = 1;

      const int blocknx1 = 16;
      const int blocknx2 = 16;
      const int blocknx3 = 16;

      LBMReal rhoLB = 0.0;//1.4e-7; //0.0;
      LBMReal nuLB = 0.168666666667;

      double coord[6];

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = 0;

      double Ri = 16;
      double f = 1;
      double r = f*Ri;

      //bounding box
      double g_minX1 = 0;
      double g_minX2 = 0;
      double g_minX3 = 0;

      double g_maxX1 = 2.0*r*4.0;
      double g_maxX2 = 2.0*r;
      double g_maxX3 = 2.0*r;

      double blockLength = (double)blocknx1*dx;

      double h = (g_maxX2)/2.0;

      double rhoLBinflow = 0.003;//9e-7;
      
    //  rhoLBinflow/=1000.0;

      double dpLB = (rhoLBinflow - rhoLB)/3.0;

      
      double dex = g_maxX1;
      double Umax = (1.0/(4.0*nuLB))*(dpLB/dex)*(h*h);

      double Re = (4*h*Umax)/(3*nuLB);

      Grid3DPtr grid(new Grid3D(comm));
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(false);
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);

      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());      

      double k1 = 4;
      double k2 = 8;

      GbObject3DPtr refineCube1_1(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2/k1-1.0, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(refineCube1_1.get(),pathname + "/geo/refineCube1_1", WbWriterVtkXmlBinary::getInstance());

      GbObject3DPtr refineCube1_2(new GbCuboid3D(g_minX1, g_maxX2-g_maxX2/k1+1.0, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(refineCube1_2.get(),pathname + "/geo/refineCube1_2", WbWriterVtkXmlBinary::getInstance());

      GbObject3DPtr refineCube2_1(new GbCuboid3D(g_minX1+2*blockLength+2*dx, g_minX2, g_minX3, g_maxX1-2*blockLength-2*dx, g_maxX2/k2-1.0, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(refineCube2_1.get(),pathname + "/geo/refineCube2_1", WbWriterVtkXmlBinary::getInstance());

      GbObject3DPtr refineCube2_2(new GbCuboid3D(g_minX1+2*blockLength+2*dx, g_maxX2-g_maxX2/k2+1.0, g_minX3, g_maxX1-2*blockLength-2*dx, g_maxX2, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(refineCube2_2.get(),pathname + "/geo/refineCube2_2", WbWriterVtkXmlBinary::getInstance());
      
      GbObject3DPtr refineCube2_3(new GbCuboid3D(g_minX1+blockLength+2*dx, g_minX2+blockLength+2*dx, g_minX3+blockLength+2*dx, g_maxX1-blockLength-2*dx, g_maxX2-blockLength-2*dx, g_maxX3-blockLength-2*dx));
      if(myid ==0) GbSystem3D::writeGeoObject(refineCube2_3.get(),pathname + "/geo/refineCube2_3", WbWriterVtkXmlBinary::getInstance());

      //UbSchedulerPtr rSch(new UbScheduler());
      //rSch->addSchedule(50, 50, 50);
      //RestartPostprocessorPtr rp(new RestartPostprocessor(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY));

      std::string opt;

      if(cstr!= NULL)
         opt = std::string(cstr);

      //////////////////////////////////////////////////////////////////////////
      //restart
      double restartStep = 300000;
      double restartStart = 200000;
      UbSchedulerPtr rSch(new UbScheduler(restartStep));
      RestartPostprocessor rp(grid, rSch, comm, pathname, RestartPostprocessor::TXT);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep() == 0)
      {
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         if(myid ==0)
         {
            UBLOG(logINFO,"Parameters:");
            UBLOG(logINFO,"h = " << h );
            UBLOG(logINFO,"rho = " << rhoLB );
            UBLOG(logINFO,"nue = " << nuLB );
            UBLOG(logINFO,"Re = " << Re );
            UBLOG(logINFO,"dx = " << dx );
            UBLOG(logINFO,"dpLB = " << dpLB );
            UBLOG(logINFO,"Umax = " << Umax );
            UBLOG(logINFO,"number of levels = " << refineLevel+1 );
            UBLOG(logINFO,"numOfThreads = " << numOfThreads );
            UBLOG(logINFO,"path = " << pathname );
            UBLOG(logINFO,"Preprozess - start");
         }

         //tube
         GbCylinder3DPtr tube(new GbCylinder3D(gridCube->getX1Minimum(),gridCube->getX2Centroid(),gridCube->getX3Centroid(),
                                               gridCube->getX1Maximum(),gridCube->getX2Centroid(),gridCube->getX3Centroid(), g_maxX3/2.0));

         coord[0] = tube->getX1Minimum();
         coord[1] = tube->getX2Minimum();
         coord[2] = tube->getX3Minimum();
         coord[3] = tube->getX1Maximum();
         coord[4] = tube->getX2Maximum();
         coord[5] = tube->getX3Maximum();

         ////walls
         //GbCuboid3DPtr addWallYmin (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_minX3-4.0*blockLength, g_maxX1+4.0*blockLength, g_minX2, g_maxX3+4.0*blockLength));
         //if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         //GbCuboid3DPtr addWallYmax (new GbCuboid3D(g_minX1-4.0*blockLength, g_maxX2, g_minX3-4.0*blockLength, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_maxX3+4.0*blockLength));
         //if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_minX3-4.0*blockLength, g_minX1, g_maxX2+4.0*blockLength, g_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow (new GbCuboid3D(g_maxX1, g_minX2-4.0*blockLength, g_minX3-4.0*blockLength, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         ////inflow
         //GbCuboid3DPtr geoInflow (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_minX3-4.0*blockLength, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_minX3));
         //if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         ////outflow
         //GbCuboid3DPtr geoOutflow (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_maxX3, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_maxX3+4.0*blockLength));
         //if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            //refineHelper.addGbObject(refineCube1_1, 1);
            //refineHelper.addGbObject(refineCube1_2, 1);
            //refineHelper.addGbObject(refineCube2_1, 2);
            //refineHelper.addGbObject(refineCube2_2, 2);
            refineHelper.addGbObject(refineCube2_3, refineLevel);
            refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         }

         int bbOption = 0; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
         D3Q27BoundaryConditionAdapterPtr bcObst2(new D3Q27SlipBCAdapter(bbOption));

         mu::Parser fct0;
         fct0.SetExpr("0.0001");
         D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct0, 0, D3Q27BCFunction::INFCONST));
         velBCAdapter->setSecondaryBcOption(2);

         D3Q27InteractorPtr tubeInt(new D3Q27Interactor(tube, grid, bcObst,Interactor3D::INVERSESOLID));

         //walls
         //D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcObst,Interactor3D::SOLID));
         //D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcObst,Interactor3D::SOLID));
         //D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, velBCAdapter,Interactor3D::SOLID));

         //inflow
         //double dp_Ph=0.1*10000.0;//dp in Bar
         //double dp_lb=dp_Ph*0.001*(nueLB*dx)*(nueLB*dx);//nue_ph=10e-6
         //if(myid == 0) UBLOG(logINFO,"dp_lb = " << dp_lb );
         //double rhoLBinflow = 3.0*(dp_lb-rhoLB);

         D3Q27BoundaryConditionAdapterPtr denBCAdapterFront(new D3Q27DensityBCAdapter(rhoLBinflow));
         denBCAdapterFront->setSecondaryBcOption(0);
         D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, denBCAdapterFront, Interactor3D::SOLID));

         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         denBCAdapter->setSecondaryBcOption(0);
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));   
         ////////////////////////////////////////////
         /////delete solid blocks
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(tubeInt);
         //intHelper.addInteractor(addWallYminInt);
         //intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	 
         //////////////////////////////////////

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept( setConnsVisitor );

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         ppblocks->update(0);
         ppblocks.reset();

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx1) * (blocknx2) * (blocknx3);
         unsigned long nod = nob * (blocknx1) * (blocknx2) * (blocknx3);
         unsigned long nodg = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);
         double needMemAll  = double(nodg*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());

         if(myid == 0)
         {
            UBLOG(logINFO,"Number of blocks = " << nob);
            UBLOG(logINFO,"Number of nodes  = " << nod);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for(int level = minInitLevel; level<=maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO,"Number of blocks for level " << level <<" = " << nobl);
               UBLOG(logINFO,"Number of nodes for level " << level <<" = " << nobl*nodb);
            }
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
         }            

         int kernelType = 2;
         LBMKernel3DPtr kernel;
         //if (kernelType == 0)
         //{
         //   rhoLB = 1.0;
         //   kernel = LBMKernel3DPtr(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));
         //}
         //else if (kernelType == 1)
         //{
         //   rhoLB = 1.0;
         //   kernel = LBMKernel3DPtr(new LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3));
         //}
         //else if (kernelType == 2)
         //{
         //   rhoLB = 0.0;
            kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));
         //}

         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());

         BoundaryConditionPtr densityBC(new NonEqDensityBoundaryCondition());
         //BoundaryConditionPtr densityBC(new NonReflectingDensityBoundaryCondition());
         //BoundaryConditionPtr noSlipBC(new HighViscosityNoSlipBoundaryCondition());
         BoundaryConditionPtr noSlipBC(new NoSlipBoundaryCondition());

         bcProc->addBC(densityBC);
         bcProc->addBC(noSlipBC);

         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //walls
         intHelper.setBC();

         BoundaryConditionBlockVisitor bcVisitor;
         grid->accept(bcVisitor);

         //initialization of distributions
         mu::Parser fct;
         fct.SetExpr("-(1.0/(2.0*nu))*(dp/dx)*((x2-h)^2 - h^2)");
         fct.DefineConst("dp", dpLB);
         fct.DefineConst("dx", dex);
         fct.DefineConst("h", h);
         fct.DefineConst("nu", nuLB);

         D3Q27ETInitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         //initVisitor.setVx1(fct);
         initVisitor.setVx1(0.0);
         //initVisitor.setVx3(fct);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
         ppgeo->update(0);
         ppgeo.reset();

         //grid->addInteractor(inflowInt);
         //grid->addInteractor(outflowInt);

         if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 
      }
      else
      {
         Grid3D::Interactor3DSet interactors = grid->getInteractors();
         interactors[0]->setGrid3D(grid);
         boost::dynamic_pointer_cast<D3Q27Interactor>(interactors[0])->deleteBCAdapter();
         D3Q27BoundaryConditionAdapterPtr denBCAdapterFront(new D3Q27DensityBCAdapter(rhoLBinflow));
         boost::dynamic_pointer_cast<D3Q27Interactor>(interactors[0])->addBCAdapter(denBCAdapterFront);
         interactors[0]->updateInteractor();

         //BOOST_FOREACH(Interactor3DPtr i, interactors)
         //{
         //   i->setGrid3D(grid);
         //   i->updateInteractor();
         //}

         BoundaryConditionBlockVisitor bcVisitor;
         grid->accept(bcVisitor);

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch,numOfThreads, comm);

      D3Q27IntegrateValuesHelperPtr ih1(new D3Q27IntegrateValuesHelper(grid, comm, coord[0], coord[1], coord[2], coord[3], coord[4], coord[5]));
      if (myid == 0) GbSystem3D::writeGeoObject(ih1->getBoundingBox().get(), pathname + "/geo/ih1", WbWriterVtkXmlBinary::getInstance());

      UbSchedulerPtr stepMV(new UbScheduler(5));

      TimeseriesPostprocessor tsp(grid, stepMV, ih1, pathname+"/ts/ts1", comm);

      double outTime = 100;
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      stepSch->addSchedule(1000, 1000, 1000);
      stepSch->addSchedule(384000, 384000, 384000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv);

      double endTime =400000;//10001.0;

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
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
//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

   if ( argv != NULL )
   {
//      if (argc > 1)
//      {
         //pflowForcing(argv[1]);
         pflowdp(argv[1]);
//      }
//      else
//      {
//         cout << "Configuration file must be set!: " <<  argv[0] << " <config file>" << endl << std::flush;
//      }
   }

   return 0;
}

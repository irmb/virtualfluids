#include <iostream>
#include <string>

#include <boost/pointer_cast.hpp>

#include "vfluids.h"

using namespace std;

#include <omp.h>

void block_test_incompressible(const char *cstr1, const char *cstr2)
{

   try
   {

      //Sleep(30000);

      ConfigFileReader cf(cstr1);
      if ( !cf.read() )
      {
         std::string exceptionText = "Unable to read configuration file\n";
         throw exceptionText;
      }

      string machine = QUOTEME(CAB_MACHINE);
      string pathname = cf.getValue("path"); 
      int numOfThreads = UbSystem::stringTo<int>(cf.getValue("numOfThreads"));
      double availMem = 0;

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if(machine == "BOMBADIL") 
      {
         //pathname = "c:/temp/block_test";
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         //pathname = "/work/koskuche/scratch/block_test";
         availMem = 12.0e9;
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");


      //if(myid ==0)
      //{
      //   UbLog::reportingLevel() = logDEBUG5;
      //   stringstream logFilename;
      //   logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(myid)+".txt";
      //   UbLog::output_policy::setStream(logFilename.str());
      //}

      double dx = 1.0;

      const int blocknx1 = UbSystem::stringTo<int>(cf.getValue("blocknx1")); //16;
      const int blocknx2 = UbSystem::stringTo<int>(cf.getValue("blocknx2"));//16;
      const int blocknx3 = UbSystem::stringTo<int>(cf.getValue("blocknx3"));//16;

      const int gridNx1 = UbSystem::stringTo<int>(cf.getValue("gridNx1"));//3;
      const int gridNx2 = UbSystem::stringTo<int>(cf.getValue("gridNx2"));//3;
      const int gridNx3 = UbSystem::stringTo<int>(cf.getValue("gridNx3"));//3;


      double L1 = gridNx1*blocknx1*dx;
      double L2, L3, H;
      L2 = L3 = H = gridNx2*blocknx1*dx;

      LBMReal radius = 3.0*dx;
      LBMReal uLB = 0.01;
      LBMReal Re = 0.5;
      LBMReal rhoLB = 0.0;
      LBMReal l = L2 / dx;
      LBMReal nuLB = (((4.0/9.0)*uLB)*2.0*(radius/dx))/Re;
      //LBMReal nueLB = 0.005842;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = UbSystem::stringTo<int>(cf.getValue("refineLevel"));

      //bounding box
      double d_minX1 = 0.0;
      double d_minX2 = 0.0;
      double d_minX3 = 0.0;

      double d_maxX1 = L1;
      double d_maxX2 = L2;
      double d_maxX3 = L3;

      double offs = 0.0;

      //double g_minX1 = d_minX1-offs-0.499999*dx;
      double g_minX1 = d_minX1-offs;
      double g_minX2 = d_minX2-offs;
      double g_minX3 = d_minX3-offs;

      double g_maxX1 = d_maxX1+offs;
      double g_maxX2 = d_maxX2+offs;
      double g_maxX3 = d_maxX3+offs;

      double blockLength = blocknx1*dx;

      //obstacle
      GbObject3DPtr cylinder(new GbCylinder3D(L1*0.5-2*blockLength, L2*0.5+dx, -1.0*dx, L1*0.5-2*blockLength, L2*0.5+dx, L3+1.0*dx, radius));
      GbSystem3D::writeGeoObject(cylinder.get(),pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

      D3Q27InteractorPtr cylinderInt;
      D3Q27InteractorPtr addWallZminInt;

      //refinement area
      double off = dx;
      GbObject3DPtr refineCube(new  GbCuboid3D(cylinder->getX1Minimum()-off, cylinder->getX2Minimum()-off, cylinder->getX3Minimum(), 
         cylinder->getX1Maximum()+off, cylinder->getX2Maximum()+off, cylinder->getX3Maximum()));

      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance()); 

      Grid3DPtr grid(new Grid3D(comm));
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);

      //Grid3DPtr grid(new Grid3D(comm, blocknx1, blocknx2, blocknx3, gridNx1, gridNx2, gridNx3));
      //grid->setPeriodicX1(true);
      //grid->setPeriodicX2(true);
      //grid->setPeriodicX3(true);

      double outTime = 1.0;
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      //PostprocessorPtr pp; //(new D3Q27MacroscopicQuantitiesPostprocessor(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlASCII::getInstance(), conv, comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(1000,10,10000));
      //RestartPostprocessor rp(grid, rSch, comm, pathname, RestartPostprocessor::BINARY);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep() == 0)
      {
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         //rp->addPostprocessor(pp);
         if(myid ==0)
         {
            UBLOG(logINFO,"Parameters:");
            UBLOG(logINFO,"L = " << L2/dx );
            UBLOG(logINFO,"v = " << uLB );
            UBLOG(logINFO,"rho = " << rhoLB );
            UBLOG(logINFO,"nue = " << nuLB );
            UBLOG(logINFO,"Re = " << Re );
            UBLOG(logINFO,"dx = " << dx );
            UBLOG(logINFO,"number of levels = " << refineLevel+1 );
            UBLOG(logINFO,"numOfThreads = " << numOfThreads );
            UBLOG(logINFO,"Preprozess - start");
         }

         if(myid ==0) GbSystem3D::writeGeoObject(refineCube.get(),pathname + "/geo/refineCube", WbWriterVtkXmlBinary::getInstance());

         //walls
         GbCuboid3DPtr addWallYmin (new GbCuboid3D(d_minX1-blockLength, d_minX2-blockLength, d_minX3-blockLength, d_maxX1+blockLength, d_minX2, d_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin (new GbCuboid3D(d_minX1-blockLength, d_minX2-blockLength, d_minX3-blockLength, d_maxX1+blockLength, d_maxX2+blockLength, d_minX3));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax (new GbCuboid3D(d_minX1-blockLength, d_maxX2, d_minX3-blockLength, d_maxX1+blockLength, d_maxX2+blockLength, d_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax (new GbCuboid3D(d_minX1-blockLength, d_minX2-blockLength, d_maxX3, d_maxX1+blockLength, d_maxX2+blockLength, d_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow (new GbCuboid3D(d_minX1-blockLength, d_minX2-blockLength, d_minX3-blockLength, d_minX1, d_maxX2+blockLength, d_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow (new GbCuboid3D(d_maxX1, d_minX2-blockLength, d_minX3-blockLength, d_maxX1+blockLength, d_maxX2+blockLength, d_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         //GbCuboid3DPtr addWallYmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_maxX2-2*dx, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         //if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         //GbCuboid3DPtr addWallZmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_maxX3-2*dx, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         //if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         ////inflow
         //GbCuboid3DPtr geoInflow (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_minX1+2*dx, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         //if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         ////outflow
         //GbCuboid3DPtr geoOutflow (new GbCuboid3D(d_maxX1-2*dx, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         //if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            refineHelper.addGbObject(refineCube, refineLevel);
            refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         }


         int bbOptionC = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObstC(new D3Q27NoSlipBCAdapter(bbOptionC));
         cylinderInt = D3Q27InteractorPtr ( new D3Q27Interactor(cylinder, grid, bcObstC,Interactor3D::SOLID));

         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
         //walls
         D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcObst,Interactor3D::SOLID));
         addWallZminInt = D3Q27InteractorPtr(new D3Q27Interactor(addWallZmin, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcObst,Interactor3D::SOLID));

         mu::Parser fct;
         //fct.SetExpr("16*U*x2*x3*(H-x2)*(H-x3)/H^4");
         //fct.SetExpr("-4*U*(x2^2-H*x2)/H^2");
         //fct.DefineConst("U", 3/2*uLB);
         //fct.DefineConst("H", H);

         fct.SetExpr("U");
         fct.DefineConst("U", uLB);

         //inflow
         D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
         velBCAdapter->setSecondaryBcOption(2);
         //double dp_Ph=0.1*10000.0;//
         //double dp_lb=dp_Ph*0.001*(nueLB*dx)*(nueLB*dx);//nue_ph=10e-6
         //if(myid == 0) UBLOG(logINFO,"dp_lb = " << dp_lb );
         //D3Q27BoundaryConditionAdapterPtr denBCAdapterFront(new D3Q27DensityBCAdapter(3.0*(dp_lb-rhoLB)));
         //denBCAdapterFront->setSecondaryBcOption(0);
         D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));
         //D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, denBCAdapterFront, Interactor3D::SOLID));

         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         denBCAdapter->setSecondaryBcOption(0);
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));
         //D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, bcObst,Interactor3D::SOLID));

         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(cylinderInt);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();

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
               UBLOG(logINFO,"Number of blocks for level " << level <<" = " << nob);
               UBLOG(logINFO,"Number of nodes for level " << level <<" = " << nob*nodb);
            }
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
         }            

         //int kernelType = UbSystem::stringTo<int>(cf.getValue("kernel"));
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
            rhoLB = 0.0;
            kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));
            //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLBWithSpongeLayer(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));
            //int nx[4];
            //nx[0]=gridNx1*blocknx1*(1<<refineLevel);
            //nx[1]=gridNx2*blocknx2*(1<<refineLevel);
            //nx[2]=gridNx3*blocknx3*(1<<refineLevel);
            //nx[3]=refineLevel+1;
            //EsoTwistD3Q27SparseData::setSize(nx);
            //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLBSparse(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLBSparse::NORMAL));
            //kernel = LBMKernel3DPtr(new LBMKernelESD3Q27CCLB(blocknx1, blocknx2, blocknx3, grid));
            //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLBex(blocknx1, blocknx2, blocknx3, 0, grid));
         //}

         //int sizeSP=2;
         //mu::Parser spongeLayer;
         //spongeLayer.SetExpr("x1>=(sizeX-sizeSP)/dx ? (sizeX-(x1+1))/sizeSP/2.0 + 0.5 : 1.0");
         //spongeLayer.DefineConst("sizeX", gridNx1*blocknx1);
         //spongeLayer.DefineConst("sizeSP", sizeSP*blocknx1);
         //kernel->setWithSpongeLayer(true);
         //kernel->setSpongeLayer(spongeLayer);



         //mu::Parser fctForcingX1;
         //fctForcingX1.SetExpr("Fx1");
         //fctForcingX1.DefineConst("Fx1", 9.99685e-7);

         //kernel->setForcingX1(fctForcingX1);
         //kernel->setWithForcing(true);
         //
         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         //BCProcessorPtr bcProc(new D3Q27ETForThinWallBCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);

         grid->accept(kernelVisitor);

         //////////////////////////////////////////////////////////////////////////
         //Experemintel
         //////////////////////////////////////////////////////////////////////////
         //int minInitLevel = grid->getCoarsestInitializedLevel();

         //for(int level = minInitLevel; level<=maxInitLevel; level++)
         //{
         //   vector<Block3DPtr> blockVector;
         //   grid->getBlocks(level, blockVector);
         //   BOOST_FOREACH(Block3DPtr block, blockVector)
         //   {
         //      if (block)
         //      {
         //         boost::dynamic_pointer_cast<LBMKernelESD3Q27CCLB>(block->getKernel())->initNeighbours();
         //      }
         //   }
         //}
         //////////////////////////////////////////////////////////////////////////


         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //UbSchedulerPtr geoSch(new UbScheduler(1));
         //D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
         // new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, comm, true));
         //ppgeo->update(0);
         //ppgeo.reset();

         //return;

         intHelper.setBC();

         //initialization of distributions
         D3Q27ETInitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         //initVisitor.setVx1(fct);
         //initVisitor.setNu(nueLB);
         //initVisitor.setVx1(0.01);
         //initVisitor.setVx2(0.02);
         //initVisitor.setVx3(0.03);
         grid->accept(initVisitor);

         //Postrozess
         //UbSchedulerPtr geoSch(new UbScheduler(1));
         //D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
         //   new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, comm, true));
         //ppgeo->update(0);
         //ppgeo.reset();

         {
            UbSchedulerPtr geoSch(new UbScheduler(1));
            //D3Q27MacroscopicQuantitiesPostprocessor ppgeo(grid,geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv,  comm, true);
            D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
               new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
            //grid->addObserver(ppgeo);
            grid->doPostProcess(0);
            //grid->notifyObservers(0);
            //grid->removeObserver(ppgeo);
         }

         //grid->notifyObservers(0);

         //UbSchedulerPtr stepSch(new UbScheduler(outTime));
         //D3Q27MacroscopicQuantitiesPostprocessorPtr pp(new D3Q27MacroscopicQuantitiesPostprocessor(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlASCII::getInstance(), conv));
         //rp->addPostprocessor(pp);

         //for (int i=0; i < 10; i++)
         //{
         //   grid->doPostProcess(i);
         //}

         //return;

         //UbSchedulerPtr rs(new UbScheduler(3));
         //D3Q27ShearStressPostprocessorPtr shsPp(new D3Q27ShearStressPostprocessor(grid,pathname + "/shs/shs", WbWriterVtkXmlASCII::getInstance(), stepSch, rs));
         //shsPp->addInteractor(boost::dynamic_pointer_cast<D3Q27Interactor>(addWallZminInt));
         //rp->addPostprocessor(shsPp);

         if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 
      }

      else
      {
         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept( setConnsVisitor );
         if(myid == 0) UBLOG(logINFO,"Restart - end"); 
      }
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      //UbSchedulerPtr visSch(new UbScheduler());
      //visSch->addSchedule(1,1,3);
      //visSch->addSchedule(100,100,1000);
      //visSch->addSchedule(1000,1000,5000);
      //visSch->addSchedule(5000,5000,100000);
      //visSch->addSchedule(100000,100000,10000000);

      
      D3Q27IntegrateValuesHelperPtr ih1(new D3Q27IntegrateValuesHelper(grid, comm, gridCube->getX1Minimum(), gridCube->getX2Minimum(), gridCube->getX3Minimum(), 
         gridCube->getX1Maximum(), gridCube->getX2Maximum(), gridCube->getX3Maximum()));
      if (myid == 0) GbSystem3D::writeGeoObject(ih1->getBoundingBox().get(), pathname + "/geo/ih1", WbWriterVtkXmlBinary::getInstance());

      double factorp = 1; // dp_real / dp_LB;
      double factorv = 1;// dx / dt;
      UbSchedulerPtr stepMV(new UbScheduler(500));

      TimeseriesPostprocessor tsp(grid, stepMV, ih1, pathname+cf.getValue("timeSeriesOut"), comm);

      UbSchedulerPtr visSch(stepSch);

      //UbSchedulerPtr avSch(new UbScheduler());
      //avSch->addSchedule(100,100,10000);
      //
      //double startStep = 32000;
      //UbSchedulerPtr resSchRMS(new UbScheduler());
      //resSchRMS->addSchedule(100000, startStep, 10000000);
      //UbSchedulerPtr resSchMeans(new UbScheduler());
      //resSchMeans->addSchedule(100000, startStep, 10000000);
      //UbSchedulerPtr stepAvSch(new UbScheduler());
      //int averageInterval=100;
      //stepAvSch->addSchedule(averageInterval,0,10000000);

      //AverageValuesPostprocessor Avpp(grid, pathname, WbWriterVtkXmlBinary::getInstance(), visSch/*wann wird rausgeschrieben*/, stepAvSch/*wann wird gemittelt*/, resSchMeans,resSchRMS/*wann wird resettet*/);

      UbSchedulerPtr emSch(new UbScheduler(100));
      //EmergencyExitPostprocessor empr(grid, emSch, pathname, RestartPostprocessorPtr(&rp), comm);

      //rp->addPostprocessor(avPp);

      //D3Q27ShearStressPostprocessor shs(grid,pathname, WbWriterVtkXmlASCII::getInstance(), stepSch, resSchMeans);
      //shs.addInteractor(boost::dynamic_pointer_cast<D3Q27Interactor>(addWallZminInt));

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname, WbWriterVtkXmlASCII::getInstance(), conv);

      //UbSchedulerPtr visSch(new UbScheduler(1));

      double endTime = UbSystem::stringTo<int>(cf.getValue("endTime"));//10001.0;

      cout << "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()<<endl;
      cout << "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed()<<endl;
      cout << "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()<<endl;

      //#pragma omp parallel num_threads(4)
      //      {
      //         int i = omp_get_thread_num();
      //         printf_s("Hello from thread %d\n", i);
      //      }

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




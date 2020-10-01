#include <vfluids.h>

#include <fbond.h>
#include <Version.h>

#include <JM.h>

using namespace std;

using namespace fetol;

////////////////////////////////////////////////////////////////////////
void chanel(const char *cstr)
{
   try
   {
      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
      int numOfThreads = 1;
      double availMem = 0;

      UBLOG(logINFO,"Communicator-init::strat");
      CommunicatorPtr comm = FETOLCommunicator::getInstance();
      UBLOG(logINFO,"Communicator-init::end");

      int myid = comm->getProcessID();
      int mybundle = comm->getBundleID();
      int root = comm->getRoot();
      int myrank = boost::dynamic_pointer_cast<FETOLCommunicator>(comm)->getMPIRank();

      UBLOG(logINFO,"myid = " << myid );
      UBLOG(logINFO,"mybundle = " << mybundle );
      UBLOG(logINFO,"myrank = " << myrank );

      JM::init(mybundle, myrank);

      //UbLog::reportingLevel() = logDEBUG5;

      if(machine == "BOMBADIL") 
      {
         pathname = "d:/temp/fetol_demo";
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         //pathname = "/work/koskuche/scratch/fetol_demo";

         pathname = string(cstr);
         availMem = 1.5e9;

        if(myid==root && mybundle==root)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/log/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(myid)+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      double dx = 1;

      const int blocknx1 = 40;
      const int blocknx2 = 40;
      const int blocknx3 = 40;

      const int gridNx1 = 4;
      const int gridNx2 = 2;
      const int gridNx3 = 2;

      double L1 = gridNx1*blocknx1;
      double L2, L3, H;
      L2 = H = gridNx2*blocknx1;
      L3 = gridNx3*blocknx1;

      LBMReal radius = 7;
      LBMReal uLB = 0.1;
      LBMReal Re = 3000.0;
      LBMReal rhoLB = 0.0;
      LBMReal l = L2 / dx;
      LBMReal nueLB = (uLB*2.0*radius)/Re;

      Grid3DPtr grid(new Grid3D(comm));
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);

      //UBLOG(logINFO,"Restart:start");
      UbSchedulerPtr restartSch(new UbScheduler(1000, 1000, 100000));
      RestartPostprocessor rp(grid, restartSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY);
      grid = rp.restart(-1);
      //UBLOG(logINFO,"Restart:end");

      if (grid->getTimeStep() == 0)
      {

         const int baseLevel = 0;
         const int refineLevel = 0;

         //obstacle
         GbObject3DPtr sphere(new GbSphere3D(L1/4.0, L2*0.5, L3*0.5, radius));
         GbSystem3D::writeGeoObject(sphere.get(),pathname + "/geo/sphere", WbWriterVtkXmlBinary::getInstance());

         D3Q27InteractorPtr sphereInt;

         //bounding box
         double d_minX1 = 0.0;
         double d_minX2 = 0.0;
         double d_minX3 = 0.0;

         double d_maxX1 = L1;
         double d_maxX2 = L2;
         double d_maxX3 = L3;

         double blockLength = blocknx1*dx;

         GbObject3DPtr gridCube(new GbCuboid3D(d_minX1, d_minX2, d_minX3, d_maxX1, d_maxX2, d_maxX3));
         if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance()); 

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         if(myid ==0)
         {
            UBLOG(logINFO,"Parameters:");
            UBLOG(logINFO,"L = " << L2/dx );
            UBLOG(logINFO,"v = " << uLB );
            UBLOG(logINFO,"rho = " << rhoLB );
            UBLOG(logINFO,"nue = " << nueLB );
            UBLOG(logINFO,"Re = " << Re );
            UBLOG(logINFO,"dx = " << dx );
            UBLOG(logINFO,"number of levels = " << refineLevel+1 );
            UBLOG(logINFO,"numOfThreads = " << numOfThreads );
            UBLOG(logINFO,"Preprozess - start");
         }

         //walls
         GbCuboid3DPtr addWallYmin (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_minX2, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_minX3));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_maxX2, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_maxX3, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_minX1, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow (new GbCuboid3D(d_maxX1, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));

         //sphere
         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
         sphereInt = D3Q27InteractorPtr ( new D3Q27Interactor(sphere, grid, bcObst,Interactor3D::SOLID));

         //walls
         D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcObst,Interactor3D::SOLID));

         mu::Parser fct;
         fct.SetExpr("16*U*x2*x3*(H-x2)*(H-x3)/H^4");
         fct.DefineConst("U", uLB);
         fct.DefineConst("H", H);

         //inflow
         D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
         velBCAdapter->setSecondaryBcOption(2);
         D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));

         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(sphereInt);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();

         ppblocks->update(0);
         ppblocks.reset();

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
         FETOLSetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nod = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);

         double needMemAll  = double(nod*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());

         if(myid == 0)
         {
            UBLOG(logINFO,"Number of blocks = " << nob);
            UBLOG(logINFO,"Number of nodes  = " << nod);
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
         }            

         LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));

         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
         grid->accept(kernelVisitor);

         intHelper.setBC();

         //initialization of distributions
         D3Q27ETInitDistributionsBlockVisitor initVisitor(nueLB, rhoLB);
         initVisitor.setVx1(fct);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, true));
         ppgeo->update(0);
         ppgeo.reset();

         if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 
      }
      else
      {
         UBLOG(logINFO,"SetConnectors - start, id="<<myid);

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
         FETOLSetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
         
         UBLOG(logINFO,"SetConnectors - end, id="<<myid); 
      }

      UbSchedulerPtr stepSch(new UbScheduler());
      stepSch->addSchedule(1000, 0, 1000000);
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv);

      //UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      //NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      double endTime = 100000;
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch, CalculationManager::FETOL));
      
      if(myid == 0) 
         UBLOG(logINFO,"Simulation-start");
      
      calculation->calculate();
      
      if(myid == 0) 
         UBLOG(logINFO,"Simulation-end");

      JM::finalize();
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
      if (argc > 1)
      {
         chanel(argv[1]);
      }
      else
      {
         cout << "Configuration file must be set!: " <<  argv[0] << " <config file>" << endl << std::flush;
      }
   }
}


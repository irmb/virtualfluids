
#include <vfluids.h>

#include "fbond.h"
#include "Version.h"


using namespace std;


//////////////////////////////////////////////////////////////////////////
void chanel(const char *cstr)
{
   try
   {
      ConfigFileReader cf(cstr);
      if ( !cf.read() )
      {
         std::string exceptionText = "Unable to read configuration file\n";
         throw exceptionText;
      }

      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
      int numOfThreads = UbSystem::stringTo<int>(cf.getValue("numOfThreads"));
      double availMem = 0;

      CommunicatorPtr comm;

      string comm_type = cf.getValue("comm");
      if(comm_type == "MPI")
         comm = vf::mpi::MPICommunicator::getInstance();
      else if(comm_type == "BOND")
         comm = BondCommunicator::getInstance();
      
      int myid = comm->getProcessID();
      int mybundle = comm->getBundleID();
      int root = comm->getRoot();

      //UbLog::reportingLevel() = logDEBUG5;
      
      
      pathname = cf.getValue("path");

      if(machine == "BOMBADIL") 
      {
         //pathname = "c:/temp/bond_test";
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         //pathname = "/work/koskuche/scratch/bond_test";
         availMem = 1.5e9;

         if(myid==root && mybundle==root)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(mybundle)+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      
      double dx = 1;

      const int blocknx1 = UbSystem::stringTo<int>(cf.getValue("blocknx1")); //16;
      const int blocknx2 = UbSystem::stringTo<int>(cf.getValue("blocknx2"));//16;
      const int blocknx3 = UbSystem::stringTo<int>(cf.getValue("blocknx3"));//16;

      const int gridNx1 = UbSystem::stringTo<int>(cf.getValue("gridNx1"));//3;
      const int gridNx2 = UbSystem::stringTo<int>(cf.getValue("gridNx2"));//3;
      const int gridNx3 = UbSystem::stringTo<int>(cf.getValue("gridNx3"));//3;


      double L1 = gridNx1*blocknx1;
      double L2, L3, H;
      L2 = L3 = H = gridNx2*blocknx1;

      LBMReal radius = 7;
      LBMReal uLB = 0.05;
      LBMReal Re = 300.0;
      LBMReal rhoLB = 1.0;
      LBMReal l = L2 / dx;
      LBMReal nueLB = (((4.0/9.0)*uLB)*2.0*(radius/dx))/Re;


      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = UbSystem::stringTo<int>(cf.getValue("refineLevel"));

      //obstacle
      GbObject3DPtr cylinder(new GbCylinder3D(L1*0.5, L2*0.5, 0, L1*0.5, L2*0.5, L3, radius));
      GbSystem3D::writeGeoObject(cylinder.get(),pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

      D3Q27InteractorPtr cylinderInt;

      //bounding box
      double d_minX1 = 0.0;
      double d_minX2 = 0.0;
      double d_minX3 = 0.0;

      double d_maxX1 = L1;
      double d_maxX2 = L2;
      double d_maxX3 = L3;

      double offs = dx;

      //double g_minX1 = d_minX1-offs-0.499999*dx;
      double g_minX1 = d_minX1-offs;
      double g_minX2 = d_minX2-offs;
      double g_minX3 = d_minX3-offs;

      double g_maxX1 = d_maxX1+offs;
      double g_maxX2 = d_maxX2+offs;
      double g_maxX3 = d_maxX3+offs;

      double blockLength = blocknx1*dx;

      //refinement area
      double off = 1;
      GbObject3DPtr refineCube(new  GbCuboid3D(cylinder->getX1Minimum()-off, cylinder->getX2Minimum()-off, cylinder->getX3Minimum(), 
         cylinder->getX1Maximum()+off, cylinder->getX2Maximum()+off, cylinder->getX3Maximum()));

      Grid3DPtr grid(new Grid3D(comm, blocknx1, blocknx2, blocknx3, gridNx1, gridNx2, gridNx3));


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

      if(myid ==0) GbSystem3D::writeGeoObject(refineCube.get(),pathname + "/geo/refineCube", WbWriterVtkXmlBinary::getInstance());

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

      if (refineLevel > 0)
      {
         if(myid == 0) UBLOG(logINFO,"Refinement - start");	
         RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
         refineHelper.addGbObject(refineCube, 1);
         refineHelper.refine();
         if(myid == 0) UBLOG(logINFO,"Refinement - end");	
      }

      
      if(comm_type == "MPI")
      {
         MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, true, numOfThreads);
         grid->accept( metisVisitor );
      }
      else if(comm_type == "BOND")
      {
         MetisPartitioningWithBundlesGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, true, numOfThreads);
         grid->accept( metisVisitor );
      }
      
      int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
      D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
      cylinderInt = D3Q27InteractorPtr ( new D3Q27Interactor(cylinder, grid, bcObst,Interactor3D::SOLID));

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

      ppblocks->update(0);
      ppblocks.reset();

      //set connectors
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());

      if(comm_type == "MPI")
      {
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
      }
      else if(comm_type == "BOND")
      {
         D3Q27BondSetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
      }

      //domain decomposition for threads
      //PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      //grid->accept(pqPartVisitor);

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

      //LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3));
      LBMKernel3DPtr kernel(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));
      //option = 0 - ohne param., option = 1 - mit param.
      //int option = 0;
      //LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, option));

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
      grid->addAndInitInteractor(addWallYminInt);
      grid->addAndInitInteractor(addWallZminInt);
      grid->addAndInitInteractor(addWallYmaxInt);
      grid->addAndInitInteractor(addWallZmaxInt);

      //obstacle
      //grid->addAndInitInteractor(cylinderInt);

      //inflow
      grid->addAndInitInteractor(inflowInt);

      //outflow
      grid->addAndInitInteractor(outflowInt);

      //initialization of distributions
      D3Q27ETInitDistributionsBlockVisitor initVisitor(nueLB, rhoLB);
      //initVisitor.setVx1(fct);
      grid->accept(initVisitor);

      //Postrozess
      UbSchedulerPtr geoSch(new UbScheduler(1));
      D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
         new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, true));
      ppgeo->update(0);
      ppgeo.reset();

      if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 

            //double outTime = 100.0;
      UbSchedulerPtr stepSch(new UbScheduler());
      //stepSch->addSchedule(10, 0, 1000);
      //nodeSch->addSchedule(1000, 1000, 10000);
      //nodeSch->addSchedule(10000, 10000, 50000);
      //nodeSch->addSchedule(100, 100, 10000);

      //D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv, comm);

      double step = UbSystem::stringTo<double>(cf.getValue("step"));
      double begin = UbSystem::stringTo<double>(cf.getValue("begin"));
      double end = UbSystem::stringTo<double>(cf.getValue("end"));
      
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      UbSchedulerPtr visSch(new UbScheduler(10.0));
      double endTime = UbSystem::stringTo<double>(cf.getValue("endTime"));;
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
      //if(myid == 0) 
      UBLOG(logINFO,"Simulation-start");
      calculation->calculate();
      //if(myid == 0) 
      UBLOG(logINFO,"Simulation-end");
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
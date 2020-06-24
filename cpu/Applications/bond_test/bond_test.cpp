#ifdef VF_BOND

#include <iostream>
#include <string>

#include <vfluids.h>

#include "fbond.h"
#include "Version.h"


using namespace std;

int agent_main();
void simulation(const char *cstr);

int main(int argc, char* argv[])
{
   int returnval = 0;
   try
   {
      bond::init();
      returnval = agent_main();
      bond::finalize();

      //CommunicatorPtr comm(new BondCommunicator());
      //cout<<"Bundle ID = "<<comm->getBundleID()<<", MPI rank = "<<comm->getProcessID()<<", root = "<<comm->getRoot()<<endl;

      //if ( argv != NULL )
      //{
      //   if (argc > 1)
      //   {
      //      simulation(argv[1]);
      //   }
      //   else
      //   {
      //      cout << "Configuration file must be set!: " <<  argv[0] << " <config file>" << endl << std::flush;
      //   }
      //}
   }
   catch(std::runtime_error& e)
   {
      std::cerr<<"\nRUNTIME ERROR: "<<e.what()<<"\n"<<std::endl;
   }
   catch(...)
   {
      std::cerr<<"unknown error"<<std::endl;
   }
   return returnval;
}


int agent_main()
{
   cout<<"\n=== bond lib info:\n"<<bond::Version::info()<<"\n===\n\n";

   // try to work around a bug in mpich (at least mpich2-1.4.1p1 and mpich2-1.5a1)
   int _mpiInitialized = (int)false;
   MPI_Initialized(&_mpiInitialized);
   if(!_mpiInitialized)
      MPI_Init(0, 0);	

   int mpi_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
   int mpi_size;
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
   cout<<"I am process "<<bond::processID()<<" of "<<bond::processCount()
      <<", bundle ID "<<bond::bundleID()<<" of "<<bond::bundleCount()
      <<", MPI rank "<<mpi_rank<<" of "<<mpi_size<<"\n";

   if(bond::processID() == 0)
   {
      try
      {
         Sleep(10000);
         // send
         vector<double> data(42);
         data[0] = 123.1;
         data[data.size()-1] = -999.1;

         int dst_rank = 1;
         int msg_tag = 42;
         cout<<"["<<bond::processID()<<"] nonblocking send ... "<<data[0]<<"..."<<data[data.size()-1]<<"\n";
         std::tr1::shared_ptr<bond::FutureSend> fus = bond::sendFuture(&data[0], data.size(), MPI_DOUBLE, dst_rank, msg_tag);

         vector<double> data2(42);
         data2[0] = 123.2;
         data2[data2.size()-1] = -999.2;
         cout<<"["<<bond::processID()<<"] blocking send ... "<<data2[0]<<"..."<<data2[data.size()-1]<<"\n";
         bond::sendComplete(&data2[0], data2.size(), MPI_DOUBLE, dst_rank, msg_tag);

         //Sleep(10000);

         fus->complete();
      }
      catch(std::runtime_error& e)
      {
         std::cerr<<"\nSEND ERROR: "<<e.what()<<"\n"<<std::endl;
      }
   }
   else
   {
      try
      {
         // receive
         vector<double> data(42);
         int src_rank = 0;
         cout<<"["<<bond::processID()<<"] nonblocking receive ...\n";
         int msg_tag = 42;
         std::tr1::shared_ptr<bond::FutureReceive> fur = bond::receiveFuture(&data[0], data.size(), MPI_DOUBLE, src_rank, msg_tag);


         //Sleep(10000);

         cout<<"["<<bond::processID()<<"] blocking receive ...\n";
         vector<double> data2(42);
         bond::receiveComplete(&data2[0], data2.size(), MPI_DOUBLE, src_rank, msg_tag);
         cout<<"received blocking "<<data2[0]<<"..."<<data2[data.size()-1]<<"\n";

         

         fur->complete();
         cout<<"received nonblocking "<<data[0]<<"..."<<data[data.size()-1]<<"\n";
      }
      catch(std::runtime_error& e)
      {
         std::cerr<<"\nRECEIVE ERROR: "<<e.what()<<"\n"<<std::endl;
      }
   }

   cout<<"process "<<bond::processID()<<" done\n";
   return 0;
}
//////////////////////////////////////////////////////////////////////////
void simulation(const char *cstr)
{
   try
   {
      ConfigFileReader cf(cstr);
      if ( !cf.read() )
      {
         std::string exceptionText = "Unable to read configuration file\n";
         throw exceptionText;
      }

      //UbLog::reportingLevel() = logDEBUG5;

      string machine = QUOTEME(CAB_MACHINE);
      string pathname = cf.getValue("path"); 
      int numOfThreads = UbSystem::stringTo<int>(cf.getValue("numOfThreads"));
      double availMem = 0;

      CommunicatorPtr comm;
      string comm_type = cf.getValue("comm");
      if(comm_type == "MPI")
         comm = MPICommunicator::getInstance();
      else if(comm_type == "BOND")
         comm = BondCommunicator::getInstance();

      int myid = comm->getProcessID();
      int mybundle = comm->getBundleID();
      int root = comm->getRoot();

      //UbLog::reportingLevel() = logDEBUG5;

      if(machine == "BOMBADIL") 
      {
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         availMem = 12.0e9;

         if(myid==root && mybundle==root)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
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
      LBMReal rhoLB = 0.0;
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

      //grid->setPeriodicX1(true);
      //grid->setPeriodicX2(true);
      //grid->setPeriodicX3(true);


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
         RefineCrossAndInsideGbObjectBlockVisitor refVisitor(refineCube, refineLevel);
         grid->accept(refVisitor);

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


      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());

      if(comm_type == "MPI")
      {
         MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, true, numOfThreads);
         grid->accept( metisVisitor );
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
      }
      else if(comm_type == "BOND")
      {
         MetisPartitioningWithBundlesGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, true, numOfThreads);
         grid->accept( metisVisitor );
         D3Q27BondSetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
         setConnsVisitor.activate();
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

      //SolidBlocksHelper sd(grid, comm);
      //sd.addInteractor(cylinderInt);
      //sd.addInteractor(addWallYminInt);
      //sd.addInteractor(addWallZminInt);
      //sd.addInteractor(addWallYmaxInt);
      //sd.addInteractor(addWallZmaxInt);
      //sd.addInteractor(inflowInt);
      //sd.addInteractor(outflowInt);

      //sd.deleteSolidBlocks();

      //grid->accept( metisVisitor );

      //grid->getBlock(0)->setBundle(0);
      //grid->getBlock(0)->setRank(0);
      //grid->getBlock(1)->setBundle(1);
      //grid->getBlock(1)->setRank(1);
      //grid->getBlock(2)->setBundle(0);
      //grid->getBlock(2)->setRank(0);
      //grid->getBlock(3)->setBundle(1);
      //grid->getBlock(3)->setRank(1);
      //grid->getBlock(4)->setBundle(0);
      //grid->getBlock(4)->setRank(0);
      //grid->getBlock(5)->setBundle(1);
      //grid->getBlock(5)->setRank(1);
      //grid->getBlock(6)->setBundle(1);
      //grid->getBlock(6)->setRank(1);
      //grid->getBlock(7)->setBundle(0);
      //grid->getBlock(7)->setRank(0);

      ppblocks->update(0);
      ppblocks.reset();

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
      //LBMKernel3DPtr kernel(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));
      //option = 0 - ohne param., option = 1 - mit param.
      int option = 0;
      LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));

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
      initVisitor.setVx1(fct);
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
      stepSch->addSchedule(100, 0, 1000);
      //nodeSch->addSchedule(1000, 1000, 10000);
      //nodeSch->addSchedule(10000, 10000, 50000);
      //nodeSch->addSchedule(100, 100, 10000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv, comm);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      UbSchedulerPtr visSch(new UbScheduler(10.0));
      double endTime = UbSystem::stringTo<int>(cf.getValue("endTime"));
      //CalculatorPtr calc = CalculatorPtr(new FETOLCalculator());
      //CalculatorPtr calc = CalculatorPtr(new Calculator());
      //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch, calc));
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
#endif












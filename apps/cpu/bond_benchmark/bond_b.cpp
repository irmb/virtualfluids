#ifdef VF_BOND

#include <iostream>
#include <string>

#include <vfluids.h>

#include "fbond.h"
#include "Version.h"

#include <stdlib.h>

using namespace std;


//////////////////////////////////////////////////////////////////////////
void periodic(const char *cstr1, const char *cstr2)
{
   try
   {
      //Sleep(10000);
      ConfigFileReader cf(cstr1);
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

      UbLog::reportingLevel() = logDEBUG5;
      system("hostname");
      
////////////////////////////////////////////// 
//       char hostname[1024];
//       hostname[1023] = '\0';
//       gethostname(hostname, 1023);
//       puts(hostname);
//       UBLOG(logINFO,"hostname = " << string(hostname) );
//////////////////////////////////////////////      

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

         if(myid==root /*&& mybundle==root*/)
         {
            //UBLOG(logINFO,"bundle = " << mybundle);
            //UBLOG(logINFO,"process ID = " << myid);
            stringstream logFilename;
            //logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(mybundle)+"_"+UbSystem::toString(myid)+"_"+comm_type+".txt";
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(comm->getNumberOfProcesses())+"p_"+comm_type+".txt";
            UbLog::output_policy::setStream(logFilename.str());
            //UbLog::reportingLevel() = logDEBUG5;
         }
      }
      else if(machine == "HICEGATE0")      
      {
         //pathname = "/work/koskuche/scratch/block_test";
         availMem = 6.0e9;

         if(myid ==0)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      } 
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      //UBLOG(logINFO,"bundle = " << mybundle);
      //UBLOG(logINFO,"process ID = " << myid);

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

      LBMReal uLB = 0.05;
      LBMReal Re = 20.0;
      LBMReal rhoLB = 0.0;
      LBMReal nueLB = 0.05842;

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

      double offs = dx;

      double g_minX1 = d_minX1-offs;
      double g_minX2 = d_minX2-offs;
      double g_minX3 = d_minX3-offs;

      double g_maxX1 = d_maxX1+offs;
      double g_maxX2 = d_maxX2+offs;
      double g_maxX3 = d_maxX3+offs;

      double blockLength = blocknx1*dx;

      //refinement area
      //double off = 1;
      //GbObject3DPtr refineCube(new  GbCuboid3D(cylinder->getX1Minimum()-off, cylinder->getX2Minimum()-off, cylinder->getX3Minimum(), 
      //   cylinder->getX1Maximum()+off, cylinder->getX2Maximum()+off, cylinder->getX3Maximum()));

      Grid3DPtr grid(new Grid3D(comm, blocknx1, blocknx2, blocknx3, gridNx1, gridNx2, gridNx3));
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(true);


      if(myid ==0)
      {
         //UBLOG(logINFO,"bundle = " << mybundle);
         //UBLOG(logINFO,"process ID = " << myid);
         UBLOG(logINFO,"Parameters:");
         UBLOG(logINFO,"Communicator =  " << comm_type);
         UBLOG(logINFO,"Grid size =  " << gridNx1);
         UBLOG(logINFO,"L = " << L1 );
         UBLOG(logINFO,"v = " << uLB );
         UBLOG(logINFO,"rho = " << rhoLB );
         UBLOG(logINFO,"nue = " << nueLB );
         UBLOG(logINFO,"Re = " << Re );
         UBLOG(logINFO,"dx = " << dx );
         UBLOG(logINFO,"number of levels = " << refineLevel+1 );
         UBLOG(logINFO,"number of threads = " << numOfThreads );
         UBLOG(logINFO,"number of processes = " << comm->getNumberOfProcesses() );
         UBLOG(logINFO,"Preprocess - start");
      }


      //if (refineLevel > 0)
      //{
      //   if(myid == 0) UBLOG(logINFO,"Refinement - start");	
      //   RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
      //   refineHelper.addGbObject(refineCube, 1);
      //   refineHelper.refine();
      //   if(myid == 0) UBLOG(logINFO,"Refinement - end");	
      //}

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
         //MetisPartitioningWithBundlesGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, true, numOfThreads);
	 MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, true, numOfThreads);
         grid->accept( metisVisitor );
         D3Q27BondSetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
         setConnsVisitor.activate();
      }

      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
      ppblocks->update(0);
      ppblocks.reset();

      unsigned long nob = grid->getNumberOfBlocks();
      int gl = 3;

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
         unsigned long nodb = (blocknx1) * (blocknx2) * (blocknx3);
         for(int level = minInitLevel; level<=maxInitLevel; level++)
         {
            unsigned long nobl = grid->getNumberOfBlocks(level);
            UBLOG(logINFO,"Number of blocks for level " << level <<" = " << nobl);
            UBLOG(logINFO,"Number of nodes for level " << level <<" = " << nobl*nodb);
         }
         UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
         UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
         UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
      }            

      LBMKernel3DPtr kernel;
      rhoLB = 0.0;
      kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));

      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);


      SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
      grid->accept(kernelVisitor);

      //if (refineLevel > 0)
      //{
      //   D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
      //   grid->accept(undefNodesVisitor);
      //}

      //initialization of distributions
      D3Q27ETInitDistributionsBlockVisitor initVisitor(nueLB, rhoLB);
      initVisitor.setVx1(0.0);
      grid->accept(initVisitor);

      if(myid == 0) UBLOG(logINFO,"Preprocess - end"); 
      
      UbSchedulerPtr stepSch(new UbScheduler());
      
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      UbSchedulerPtr visSch(new UbScheduler());
      double endTime = UbSystem::stringTo<int>(cf.getValue("endTime"));//10001.0;

      if(myid == 0)
      {
         UBLOG(logINFO,"//////////////////////////////////////////////////////////////////////////");
         UBLOG(logINFO,"System information:");
         UBLOG(logINFO,"Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()/1073741824.0<< " GB");
         UBLOG(logINFO,"Physical Memory currently used: " << Utilities::getPhysMemUsed()/1073741824.0<<" GB");
         UBLOG(logINFO,"Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         UBLOG(logINFO,"//////////////////////////////////////////////////////////////////////////");
      }

      CalculationManagerPtr calculation;
      if(comm_type == "MPI")
         calculation = CalculationManagerPtr(new CalculationManager(grid, numOfThreads, endTime, stepSch, CalculationManager::MPI));
      else if(comm_type == "BOND")
         calculation = CalculationManagerPtr(new CalculationManager(grid, numOfThreads, endTime, stepSch, CalculationManager::MPI));
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
   int returnval = 0;
   try
   {
      if ( argv != NULL )
      {
         if (argc > 1)
         {
            //chanel(argv[1]);
            periodic(argv[1], argv[2]);
         }
         else
         {
            cout << "Configuration file must be set!: " <<  argv[0] << " <config file>" << endl << std::flush;
         }
      }
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

#endif










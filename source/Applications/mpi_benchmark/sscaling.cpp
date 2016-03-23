#include <iostream>
#include <string>

#include "vfluids.h"

using namespace std;

int main(int argc, char* argv[])
{
   //// Initialize the MPI environment
   //MPI_Init(NULL, NULL);

   //// Get the number of processes
   //int world_size;
   //MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   //// Get the rank of the process
   //int world_rank;
   //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   //// Get the name of the processor
   //char processor_name[MPI_MAX_PROCESSOR_NAME];
   //int name_len;
   //MPI_Get_processor_name(processor_name, &name_len);

   //// Print off a hello world message
   //printf("Hello world from processor %s, rank %d"
   //   " out of %d processors\n",
   //   processor_name, world_rank, world_size);

   //// Finalize the MPI environment.
   //MPI_Finalize();

   //Sleep(30000);

   CommunicatorPtr comm = MPICommunicator::getInstance();
   int myid = comm->getProcessID();
   
   // Get the name of the processor
   char machinename[MPI_MAX_PROCESSOR_NAME];
   int name_len;
   MPI_Get_processor_name(machinename, &name_len);

   try
   {
      double availMem = 1.2e9;
      int numOfThreads = 1;

      //UbLog::reportingLevel() = UbLog::logLevelFromString("DEBUG5");

      stringstream logFilename;
      logFilename <<  "logfile_"+UbSystem::toString(machinename)+"_PID_"+UbSystem::toString(myid)+".txt";
      UbLog::output_policy::setStream(logFilename.str());

      UBLOG(logINFO, "MPI benchmark");
      UBLOG(logINFO, "1. PID = " << myid << " host name: " << machinename);
      UBLOG(logINFO, "1. PID = " << myid << " Number of processes = " << comm->getNumberOfProcesses());
      UBLOG(logINFO, "1. PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()/1073741824.0<< " GB");
      UBLOG(logINFO, "1. PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed()/1073741824.0<< " GB");
      UBLOG(logINFO, "1. PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()/1073741824.0<< " GB");

      double dx = 1;

      const int blocknx1 = 64;
      const int blocknx2 = 64;
      const int blocknx3 = 64;

      int gs = 60; // 30;
      const int gridNx1 = gs; // *comm->getNumberOfProcesses();
      const int gridNx2 = gs;
      const int gridNx3 = gs;


      double L1 = gridNx1*blocknx1;
      double L2, L3, H;
      L2 = L3 = H = gridNx2*blocknx1;

      LBMReal uLB = 0.05;
      LBMReal Re = 20.0;
      LBMReal rhoLB = 0.0;
      LBMReal l = L2 / dx;

      LBMReal nueLB = 0.05842;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      double blockLength = blocknx1*dx;

      Grid3DPtr grid(new Grid3D(comm, blocknx1, blocknx2, blocknx3, gridNx1, gridNx2, gridNx3));
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(true);

      UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
      UBLOG(logINFO, "2. PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()/1073741824.0<< " GB");
      UBLOG(logINFO, "2. PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed()/1073741824.0<< " GB");
      UBLOG(logINFO, "2. PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()/1073741824.0<< " GB");
      UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");

      UBLOG(logINFO, "MetisPartitioningGridVisitor:start");
      MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW);
      grid->accept(metisVisitor);
      UBLOG(logINFO, "MetisPartitioningGridVisitor:end");

      //set connectors
      UBLOG(logINFO, "D3Q27SetConnectorsBlockVisitor:start");
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
      grid->accept(setConnsVisitor);
      UBLOG(logINFO, "D3Q27SetConnectorsBlockVisitor:end");

      //domain decomposition for threads
      UBLOG(logINFO, "PQueuePartitioningGridVisitor:start");
      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      grid->accept(pqPartVisitor);
      UBLOG(logINFO, "PQueuePartitioningGridVisitor:end");

      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), ".", WbWriterVtkXmlBinary::getInstance(), comm));
      ppblocks->update(0);
      ppblocks.reset();

      unsigned long nob = grid->getNumberOfBlocks();
      int gl = 3;
      unsigned long nodb = (blocknx1)* (blocknx2)* (blocknx3);
      unsigned long nod = nob * (blocknx1)* (blocknx2)* (blocknx3);
      unsigned long nodg = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);
      double needMemAll = double(nodg*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
      double needMem = needMemAll / double(comm->getNumberOfProcesses());


      UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
      UBLOG(logINFO, "Setup information:");
      UBLOG(logINFO, "Number of blocks = " << nob);
      UBLOG(logINFO, "Number of nodes  = " << nod);
      int minInitLevel = grid->getCoarsestInitializedLevel();
      int maxInitLevel = grid->getFinestInitializedLevel();
      for (int level = minInitLevel; level<=maxInitLevel; level++)
      {
         int nobl = grid->getNumberOfBlocks(level);
         UBLOG(logINFO, "Number of blocks for level " << level <<" = " << nob);
         UBLOG(logINFO, "Number of nodes for level " << level <<" = " << nob*nodb);
      }
      UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
      UBLOG(logINFO, "Necessary memory  = " << needMemAll/1073741824.0  << " GB");
      UBLOG(logINFO, "Necessary memory per process = " << needMem/1073741824.0  << " GB");
      UBLOG(logINFO, "Available memory per process = " << availMem/1073741824.0 << " GB");
      UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");


      LBMKernel3DPtr kernel;
      kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));
      //rhoLB = 1.0;
      //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));

      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
      grid->accept(kernelVisitor);


      //initialization of distributions
      D3Q27ETInitDistributionsBlockVisitor initVisitor(nueLB, rhoLB);
      initVisitor.setVx1(0.0);
      grid->accept(initVisitor);


      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      UbSchedulerPtr visSch(new UbScheduler(500, 500));

      double endTime = 100;

      UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
      UBLOG(logINFO, "System information:");
      UBLOG(logINFO, "Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()/1073741824.0<< " GB");
      UBLOG(logINFO, "Physical Memory currently used: " << Utilities::getPhysMemUsed()/1073741824.0<<" GB");
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
      UBLOG(logINFO, "Simulation-end");
   }
   catch (std::exception& e)
   {
      cerr << "PID = " << myid << " host name: " << machinename << endl << flush;
      cerr << e.what() << endl << flush<<
         boost::current_exception_diagnostic_information();
   }
   catch (std::string& s)
   {
      cerr << s << endl<<boost::current_exception_diagnostic_information();
   }
   catch (...)
   {
      cerr << "unknown exception" << endl<<
         boost::current_exception_diagnostic_information();
   }

   return 0;
}




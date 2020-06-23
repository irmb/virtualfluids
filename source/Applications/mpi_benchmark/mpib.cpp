#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;


void run(string configname)
{
   SPtr<Communicator> comm = MPICommunicator::getInstance();
   int myid = comm->getProcessID();

   // Get the name of the processor
   char machinename[MPI_MAX_PROCESSOR_NAME];
   int name_len;
   MPI_Get_processor_name(machinename, &name_len);

   try
   {
      //UbLog::reportingLevel() = UbLog::logLevelFromString("DEBUG5");

      ConfigurationFile   config;
      config.load(configname);

      string          pathOut = config.getString("pathOut");
      double          endTime = config.getDouble("endTime");
      int             numOfThreads = config.getInt("numOfThreads");
      vector<int>     blockNx = config.getVector<int>("blockNx");
      double          availMem = config.getDouble("availMem");
      bool            logToFile = config.getBool("logToFile");
      bool            oneD = config.getBool("oneD");
      bool            output = config.getBool("output");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");
      bool            priorityQueue = config.getBool("priorityQueue");
      bool            restart = config.getBool("restart");
      double          restartStep = config.getDouble("restartStep");
      double          cpStep = config.getDouble("cpStep");

      if (logToFile)
      {
#if defined(__unix__)
         if (myid==0)
         {
            const char* str = pathOut.c_str();
            mkdir(str, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
         }
#endif 

         if (myid==0)
         {
            stringstream logFilename;
            logFilename<<pathOut+"/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }

      if (myid==0)
      {
         UBLOG(logINFO, "MPI benchmark");
         UBLOG(logINFO, "1. PID = "<<myid<<" host name: "<<machinename);
         UBLOG(logINFO, "1. PID = "<<myid<<" Number of processes = "<<comm->getNumberOfProcesses());
         UBLOG(logINFO, "1. PID = "<<myid<<" Number of threads = "<<numOfThreads);
         UBLOG(logINFO, "1. PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem()/1073741824.0<<" GB");
         UBLOG(logINFO, "1. PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed()/1073741824.0<<" GB");
         UBLOG(logINFO, "1. PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");

      }

      LBMReal uLB = 0.05;
      LBMReal Re = 20.0;
      LBMReal rhoLB = 0.0;
      LBMReal nueLB = 0.05842;
      
      SPtr<Grid3D> grid(new Grid3D(comm));
      
      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> rSch(new UbScheduler(cpStep,cpStep));
      MPIIORestartCoProcessor rcp(grid, rSch, pathOut, comm);

      if (restart)
      {
         rcp.restart((int)restartStep);
      }
      else
      {
         double dx = 1;
         double g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3;
         double factor = 1.0;

         if (oneD)
         {
            factor = comm->getNumberOfProcesses() * numOfThreads;
            g_minX1 = 0;
            g_minX2 = 0;
            g_minX3 = 0;

            g_maxX1 = blockNx[0]*2.0 * factor;
            g_maxX2 = blockNx[1]*2.0;
            g_maxX3 = blockNx[2]*2.0;
         }
         else
         {
            factor = pow(comm->getNumberOfProcesses() * numOfThreads, 1.0/3.0);
            g_minX1 = 0;
            g_minX2 = 0;
            g_minX3 = 0;

            g_maxX1 = blockNx[0]*2.0 * factor;
            g_maxX2 = blockNx[1]*2.0 * factor;
            g_maxX3 = blockNx[2]*2.0 * factor;
         }

         SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

         grid->setDeltaX(dx);
         grid->setBlockNX(blockNx[0], blockNx[1], blockNx[2]);

         SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid==0&&output) GbSystem3D::writeGeoObject(gridCube.get(), pathOut+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         //grid->setPeriodicX1(true);
         //grid->setPeriodicX2(true);
         //grid->setPeriodicX3(true);

         if (myid==0)
         {
            UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
            UBLOG(logINFO, "2. PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem()/1073741824.0<<" GB");
            UBLOG(logINFO, "2. PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed()/1073741824.0<<" GB");
            UBLOG(logINFO, "2. PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
            UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
         }

         if (priorityQueue)
         {
            if (myid==0) UBLOG(logINFO, "MetisPartitioningGridVisitor:start");
            MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW);
            grid->accept(metisVisitor);
            if (myid==0) UBLOG(logINFO, "MetisPartitioningGridVisitor:end");

            //domain decomposition for threads
            if (myid==0) UBLOG(logINFO, "PQueuePartitioningGridVisitor:start");
            PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
            grid->accept(pqPartVisitor);
            if (myid==0) UBLOG(logINFO, "PQueuePartitioningGridVisitor:end");
         }
         else
         {
            if (myid==0) UBLOG(logINFO, "MetisPartitioningGridVisitor:start");
            MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY, true, numOfThreads);
            grid->accept(metisVisitor);
            if (myid==0) UBLOG(logINFO, "MetisPartitioningGridVisitor:end");
         }


         if (output)
         {
            WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
            ppblocks.process(0);
         }

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blockNx[0])* (blockNx[1])* (blockNx[2]);
         unsigned long nod = nob * (blockNx[0])* (blockNx[1])* (blockNx[2]);
         unsigned long nodg = nob * (blockNx[0]+gl) * (blockNx[1]+gl) * (blockNx[2]+gl);
         double needMemAll = double(nodg*(27*sizeof(double)+sizeof(int)+sizeof(float)*4));
         double needMem = needMemAll/double(comm->getNumberOfProcesses());

         if (myid==0)
         {
            UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
            UBLOG(logINFO, "Setup information:");
            UBLOG(logINFO, "Size of block = "<<blockNx[0]<<" x "<<blockNx[1]<<" x "<<blockNx[2]<<" nodes");
            UBLOG(logINFO, "Size of domain = "<<g_maxX1<<" x "<<g_maxX2<<" x "<<g_maxX3<<" dx ");
            UBLOG(logINFO, "Number of blocks = "<<nob);
            UBLOG(logINFO, "Number of nodes  = "<<nod);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for (int level = minInitLevel; level<=maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO, "Number of blocks for level "<<level<<" = "<<nob);
               UBLOG(logINFO, "Number of nodes for level "<<level<<" = "<<nob*nodb);
            }
            UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
            UBLOG(logINFO, "Necessary memory  = "<<needMemAll/1073741824.0<<" GB");
            UBLOG(logINFO, "Necessary memory per process = "<<needMem/1073741824.0<<" GB");
            UBLOG(logINFO, "Available memory per process = "<<availMem/1073741824.0<<" GB");
            UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
         }

         SPtr<LBMKernel> kernel;
         kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel(blockNx[0], blockNx[1], blockNx[2], IncompressibleCumulantLBMKernel::NORMAL));

         SPtr<BCProcessor> bcProc(new BCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
         grid->accept(kernelVisitor);

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor(nueLB, rhoLB);
         initVisitor.setVx1(uLB);
         grid->accept(initVisitor);
      }

      //set connectors
      if (myid==0) UBLOG(logINFO, "SetConnectorsBlockVisitor:start");
      InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
      SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
      grid->accept(setConnsVisitor);
      if (myid==0) UBLOG(logINFO, "SetConnectorsBlockVisitor:end");

      SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      SPtr<UbScheduler> visSch(new UbScheduler(500, 500));

      if (myid==0)
      {
         UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
         UBLOG(logINFO, "System information:");
         UBLOG(logINFO, "Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem()/1073741824.0<<" GB");
         UBLOG(logINFO, "Physical Memory currently used: "<<Utilities::getPhysMemUsed()/1073741824.0<<" GB");
         UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         UBLOG(logINFO, "//////////////////////////////////////////////////////////////////////////");
      }

      const SPtr<ConcreteCalculatorFactory> calculatorFactory = std::make_shared<ConcreteCalculatorFactory>(visSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, calculatorFactory, CalculatorType::MPI));
      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
      if (myid==0) UBLOG(logINFO, "Simulation-end");
   }
   catch (std::exception& e)
   {
      cerr<<"PID = "<<myid<<" host name: "<<machinename<<endl<<flush;
      cerr<<e.what()<<endl<<flush<<
         boost::current_exception_diagnostic_information();
   }
   catch (std::string& s)
   {
      cerr<<s<<endl<<boost::current_exception_diagnostic_information();
   }
   catch (...)
   {
      cerr<<"unknown exception"<<endl<<
         boost::current_exception_diagnostic_information();
   }
}

int main(int argc, char* argv[])
{

   if (argv!=NULL)
   {
      if (argv[1]!=NULL)
      {
         run(string(argv[1]));
         UBLOG(logINFO, "run end");
      }
      else
      {
         cout<<"Configuration file must be set!: "<<argv[0]<<" <config file>"<<endl<<std::flush;
      }
   }

   return 0;
}




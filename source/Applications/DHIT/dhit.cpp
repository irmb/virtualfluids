#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;


void run(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathname = config.getString("pathname");
      int             numOfThreads = config.getInt("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
      double          endTime = config.getDouble("endTime");
      double          outTime = config.getDouble("outTime");
      double          availMem = config.getDouble("availMem");
      vector<double>  length = config.getVector<double>("length");
      bool            logToFile = config.getBool("logToFile");
      string          initFile = config.getString("initFile");
      double          nuLB = config.getDouble("nuLB");
      double          uRMS = config.getDouble("uRMS");
      double          lambda = config.getDouble("lambda");
      double          initTime = config.getDouble("initTime");

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (logToFile)
      {
#if defined(__unix__)
         if (myid == 0)
         {
            const char* str = pathname.c_str();
            mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         }
#endif 

         if (myid == 0)
         {
            stringstream logFilename;
            logFilename << pathname + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }

      //LBMReal uLB = 0.032;
      LBMReal dx = 1.0;
      LBMReal rhoLB = 0.0;



      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      //bounding box
      double g_minX1 = 0.0;
      double g_minX2 = 0.0;
      double g_minX3 = 0.0;

      double g_maxX1 = length[0];//-1.0;
      double g_maxX2 = length[1];//-1.0;
      double g_maxX3 = length[2];//-1.0;

      //geometry
      GbObject3DPtr box(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(box.get(), pathname + "/geo/box", WbWriterVtkXmlBinary::getInstance());

      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


      double blockLength = blocknx[0] * dx;

      Grid3DPtr grid(new Grid3D(comm));

      if (myid == 0)
      {
         //UBLOG(logINFO, "uLb = " << uLB);
         UBLOG(logINFO, "rho = " << rhoLB);
         UBLOG(logINFO, "nuLb = " << nuLB);
         UBLOG(logINFO, "uRMS = " << uRMS);
         UBLOG(logINFO, "lambda = " << lambda);
         UBLOG(logINFO, "Re = " << (uRMS*lambda)/nuLB);
         UBLOG(logINFO, "dx = " << dx);
         UBLOG(logINFO, "length = " << length[0] << " " << length[1] << " " << length[2]);
         UBLOG(logINFO, "blocknx = " << blocknx[0] << " " << blocknx[1] << " " << blocknx[2]);
         UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
         UBLOG(logINFO, "number of threads = " << numOfThreads);
         UBLOG(logINFO, "Preprocess - start");
      }

      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(true);

      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

      Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW));
      InteractorsHelper intHelper(grid, metisVisitor);
      //intHelper.addInteractor(boxInt);
      intHelper.selectBlocks();

      ppblocks->process(0);
      ppblocks.reset();

      //set connectors
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
      grid->accept(setConnsVisitor);

      //domain decomposition for threads
      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      grid->accept(pqPartVisitor);

      unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
      int ghostLayer = 3;
      unsigned long long numberOfNodesPerBlock = (unsigned long long)(blocknx[0])* (unsigned long long)(blocknx[1])* (unsigned long long)(blocknx[2]);
      unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
      unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
      double needMemAll = double(numberOfNodesPerBlockWithGhostLayer*(27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
      double needMem = needMemAll / double(comm->getNumberOfProcesses());

      if (myid == 0)
      {
         UBLOG(logINFO, "Number of blocks = " << numberOfBlocks);
         UBLOG(logINFO, "Number of nodes  = " << numberOfNodes);
         int minInitLevel = grid->getCoarsestInitializedLevel();
         int maxInitLevel = grid->getFinestInitializedLevel();
         for (int level = minInitLevel; level <= maxInitLevel; level++)
         {
            int nobl = grid->getNumberOfBlocks(level);
            UBLOG(logINFO, "Number of blocks for level " << level << " = " << nobl);
            UBLOG(logINFO, "Number of nodes for level " << level << " = " << nobl*numberOfNodesPerBlock);
         }
         UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
         UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
         UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
      }

      LBMKernel3DPtr kernel;

      //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));
      kernel = LBMKernel3DPtr(new InitDensityLBMKernel(blocknx[0], blocknx[1], blocknx[2]));

      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
      grid->accept(kernelVisitor);

      intHelper.setBC();

      //initialization of distributions
      //D3Q27ETInitDistributionsBlockVisitor initVisitor(nuLB, rhoLB, uLB, uLB, uLB);
      InitDistributionsFromFileBlockVisitor initVisitor(nuLB, rhoLB, initFile);
      grid->accept(initVisitor);

      //boundary conditions grid
      {
         UbSchedulerPtr geoSch(new UbScheduler(1));
         MacroscopicQuantitiesCoProcessorPtr ppgeo(
            new MacroscopicQuantitiesCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
         grid->coProcess(0);
      }

      if (myid == 0) UBLOG(logINFO, "Preprocess - end");

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      UbSchedulerPtr outputSch(new UbScheduler(outTime));
      MacroscopicQuantitiesCoProcessor pp(grid, outputSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      CalculationManagerPtr initialisation(new CalculationManager(grid, numOfThreads, initTime, outputSch));
      if (myid == 0) UBLOG(logINFO, "Initialisation-start");
      initialisation->calculate();
      if (myid == 0) UBLOG(logINFO, "Initialisation-end");


      kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));
      kernel->setBCProcessor(bcProc);
      SetKernelBlockVisitor kernelVisitor2(kernel, nuLB, availMem, needMem, SetKernelBlockVisitor::Change);
      grid->accept(kernelVisitor2);

      UbSchedulerPtr visSch(new UbScheduler(outTime));

      if (myid==0) UBLOG(logINFO, "Simulation-start");
      grid->setTimeStep(initTime+1.0);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      calculation->calculate();
      if (myid==0) UBLOG(logINFO, "Simulation-end");

   }
   catch (std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch (std::string& s)
   {
      cerr << s << endl;
   }
   catch (...)
   {
      cerr << "unknown exception" << endl;
   }

}
int main(int argc, char* argv[])
{
   if (argv != NULL)
   {
      if (argv[1] != NULL)
      {
         run(string(argv[1]));
      }
      else
      {
         cout << "Configuration file is missing!" << endl;
      }
   }

}


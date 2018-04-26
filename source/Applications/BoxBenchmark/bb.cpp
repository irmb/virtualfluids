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
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
      double          uLB = config.getValue<double>("uLB");
      double          endTime = config.getValue<double>("endTime");
      double          outTime = config.getValue<double>("outTime");
      double          availMem = config.getValue<double>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      double          Re = config.getValue<double>("Re");
      double          dx = config.getValue<double>("dx");
      vector<double>  length = config.getVector<double>("length");
      bool            logToFile = config.getValue<bool>("logToFile");

      double          cpStep      = config.getValue<double>("cpStep");
      double          cpStepStart = config.getValue<double>("cpStepStart");
      bool            restart     = config.getValue<bool>("restart");
      double          restartStep = config.getValue<double>("restartStep");

      //UbLog::reportingLevel() = UbLog::logLevelFromString("DEBUG3");

      SPtr<Communicator> comm = MPICommunicator::getInstance();
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

      LBMReal dLB = length[1] / dx;
      LBMReal rhoLB = 0.0;
      LBMReal nuLB = (uLB*dLB) / Re;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      //bounding box
      double g_minX1 = 0.0;
      double g_minX2 = 0.0;
      double g_minX3 = 0.0;

      double g_maxX1 = length[0];
      double g_maxX2 = length[1];
      double g_maxX3 = length[2];

      //geometry
      SPtr<GbObject3D> box(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      //if (myid == 0) GbSystem3D::writeGeoObject(box.get(), pathname + "/geo/box", WbWriterVtkXmlBinary::getInstance());

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      //if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


      double blockLength = blocknx[0] * dx;

      SPtr<Grid3D> grid(new Grid3D(comm));

      SPtr<UbScheduler> rSch2(new UbScheduler(cpStep, cpStepStart));
      MPIIORestartCoProcessor rcp(grid, rSch2, pathname, comm);

      if (!restart)
      {
         if (myid == 0)
         {
            UBLOG(logINFO, "uLb = " << uLB);
            UBLOG(logINFO, "rho = " << rhoLB);
            UBLOG(logINFO, "nuLb = " << nuLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << dx);
            UBLOG(logINFO, "length = " << length[0] << " " << length[1] << " " << length[2]);
            UBLOG(logINFO, "blocknx = " << blocknx[0] << " " << blocknx[1] << " " << blocknx[2]);
            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
            UBLOG(logINFO, "number of threads = " << numOfThreads);
            UBLOG(logINFO, "Preprocess - start");
         }

         grid->setDeltaX(dx);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

         //if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         SPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         //int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         //D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
         //SPtr<D3Q27Interactor> boxInt(new D3Q27Interactor(box, grid, bcObst, Interactor3D::INVERSESOLID));

         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         InteractorsHelper intHelper(grid, metisVisitor);
         //intHelper.addInteractor(boxInt);
         intHelper.selectBlocks();

         ppblocks->process(0);
         ppblocks.reset();

         //set connectors
         //InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         //SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         //grid->accept(setConnsVisitor);

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

         SPtr<LBMKernel> kernel;
         kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel());
         SPtr<BCProcessor> bcProc(new BCProcessor());
         kernel->setBCProcessor(bcProc);
         kernel->setForcingX1(0.1);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         intHelper.setBC();

         //BoundaryConditionBlockVisitor bcVisitor;
         //grid->accept(bcVisitor);

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor;
         initVisitor.setVx1(0.5);
         grid->accept(initVisitor);

         if (myid == 0) UBLOG(logINFO, "Preprocess - end");
      }
      else
      {
         rcp.restart((int)restartStep);
         grid->setTimeStep(restartStep);
         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         if (myid==0) UBLOG(logINFO, "Restart - end");
      }


      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      SPtr<CoProcessor> mqCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, pathname, WbWriterVtkXmlASCII::getInstance(), conv, comm));

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      SPtr<CoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      omp_set_num_threads(numOfThreads);
      SPtr<Calculator> calculator(new BasicCalculator(grid, visSch, (int)endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(mqCoProcessor);
      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      if (myid == 0) UBLOG(logINFO, "Simulation-end");
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


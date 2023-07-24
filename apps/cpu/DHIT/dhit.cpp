#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;


void run(string configname)
{
   try
   {
      //Sleep(30000);

      vf::basics::ConfigurationFile   config;
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

      SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
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



      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      //bounding box
      double g_minX1 = 0.0;
      double g_minX2 = 0.0;
      double g_minX3 = 0.0;

      double g_maxX1 = length[0];//-1.0;
      double g_maxX2 = length[1];//-1.0;
      double g_maxX3 = length[2];//-1.0;

      //geometry
      SPtr<GbObject3D> box(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(box.get(), pathname + "/geo/box", WbWriterVtkXmlBinary::getInstance());

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


      double blockLength = blocknx[0] * dx;

      SPtr<Grid3D> grid(new Grid3D(comm));

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

      WriteBlocksSPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

      SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW));
      InteractorsHelper intHelper(grid, metisVisitor);
      //intHelper.addInteractor(boxInt);
      intHelper.selectBlocks();

      ppblocks->process(0);
      ppblocks.reset();

      //set connectors
      InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolator());
      SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
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

      SPtr<LBMKernel> kernel;

      //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));
      kernel = SPtr<LBMKernel>(new InitDensityLBMKernel(blocknx[0], blocknx[1], blocknx[2]));

      SPtr<BCProcessor> bcProc(new BCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
      grid->accept(kernelVisitor);

      intHelper.setBC();

      //initialization of distributions
      InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
      double u_LB = 0.01;
      mu::Parser inflowProfileVx1, inflowProfileVx2, inflowProfileVx3;
      inflowProfileVx1.DefineConst("U", u_LB);
      inflowProfileVx1.DefineConst("PI", PI);
      inflowProfileVx1.DefineConst("L1", g_maxX1-g_minX1);
      inflowProfileVx1.DefineConst("L2", g_maxX2-g_minX2);
      inflowProfileVx1.DefineConst("L3", g_maxX3-g_minX3);
      inflowProfileVx1.SetExpr("U*cos(2.0*PI*x1/L1)*sin(2.0*PI*x2/L2)*sin(2.0*PI*x3/L3)");
      inflowProfileVx2.DefineConst("U", u_LB);
      inflowProfileVx2.DefineConst("PI", PI);
      inflowProfileVx2.DefineConst("L1", g_maxX1-g_minX1);
      inflowProfileVx2.DefineConst("L2", g_maxX2-g_minX2);
      inflowProfileVx2.DefineConst("L3", g_maxX3-g_minX3);
      inflowProfileVx2.SetExpr("-U*cos(2.0*PI*x1/L1)*sin(2.0*PI*x2/L2)*cos(2.0*PI*x3/L3)");
      inflowProfileVx3.DefineConst("U", u_LB);
      inflowProfileVx3.DefineConst("PI", PI);
      inflowProfileVx3.DefineConst("L1", g_maxX1-g_minX1);
      inflowProfileVx3.DefineConst("L2", g_maxX2-g_minX2);
      inflowProfileVx3.DefineConst("L3", g_maxX3-g_minX3);
      inflowProfileVx3.SetExpr("-U/2.0*sin(8.0*PI*(x1)/(L1))*cos(8.0*PI*(x3)/L3)");
      initVisitor.setVx1(inflowProfileVx1);
      initVisitor.setVx2(inflowProfileVx2);
      initVisitor.setVx3(inflowProfileVx3);
      //InitDistributionsFromFileBlockVisitor initVisitor(nuLB, rhoLB, initFile);
      grid->accept(initVisitor);

      //boundary conditions grid
      {
         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         WriteBoundaryConditionsSPtr<CoProcessor> ppgeo(
            new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         grid->coProcess(0);
      }

      if (myid == 0) UBLOG(logINFO, "Preprocess - end");

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      SPtr<UbScheduler> outputSch(new UbScheduler(outTime));
      WriteMacroscopicQuantitiesCoProcessor pp(grid, outputSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      const SPtr<ConcreteCalculatorFactory> calculatorFactory = std::make_shared<ConcreteCalculatorFactory>(outputSch);
      CalculationManagerPtr initialisation(new CalculationManager(grid, numOfThreads, endTime, calculatorFactory, CalculatorType::HYBRID));
      if (myid == 0) UBLOG(logINFO, "Initialisation-start");
      initialisation->calculate();
      if (myid == 0) UBLOG(logINFO, "Initialisation-end");


      kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], IncompressibleCumulantLBMKernel::NORMAL));
      kernel->setBCProcessor(bcProc);
      SetKernelBlockVisitor kernelVisitor2(kernel, nuLB, availMem, needMem, SetKernelBlockVisitor::ChangeKernel);
      grid->accept(kernelVisitor2);

      SPtr<UbScheduler> visSch(new UbScheduler(outTime));

      if (myid==0) UBLOG(logINFO, "Simulation-start");
      grid->setTimeStep(initTime+1.0);

      const SPtr<ConcreteCalculatorFactory> calculatorFactory2 = std::make_shared<ConcreteCalculatorFactory>(visSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, calculatorFactory2, CalculatorType::HYBRID));
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


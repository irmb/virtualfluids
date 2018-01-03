#include "pf.h"
#include "VirtualFluids.h"

using namespace std;

//pipe flow with forcing
void pf1()
{
   CommunicatorPtr comm = MPICommunicator::getInstance();
   int myid = comm->getProcessID();

   //parameters
   string          pathname = "d:/temp/pflow_pipe_forcing";
   int             numOfThreads = 1;
   int             blocknx[3] ={ 10,10,10 };
   double          endTime = 100000;
   double          outTime = 100;
   double          availMem = 8e9;
   double          deltax = 1;
   double          rhoLB = 0.0;
   double          nuLB = 0.005;

   //geometry definition

   //simulation bounding box
   double g_minX1 = 0.0;
   double g_minX2 = -5.0;
   double g_minX3 = -5.0;

   double g_maxX1 = 30;
   double g_maxX2 = 5;
   double g_maxX3 = 5;

   GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
   if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

   //cylinder
   GbObject3DPtr cylinder(new GbCylinder3D(g_minX1 - 2.0*deltax, 0.0, 0.0, g_maxX1 + 2.0*deltax, 0.0, 0.0, g_maxX2));
   GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

   if (myid == 0)
   {
      UBLOG(logINFO, "rhoLB = " << rhoLB);
      UBLOG(logINFO, "nuLB = " << nuLB);
      UBLOG(logINFO, "deltaX = " << deltax);
      UBLOG(logINFO, "Preprocess - start");
   }

   //Grid definition
   Grid3DPtr grid(new Grid3D(comm));
   grid->setDeltaX(deltax);
   grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
   grid->setPeriodicX1(true);
   grid->setPeriodicX2(false);
   grid->setPeriodicX3(false);

   if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

   //blocks generating
   GenBlocksGridVisitor genBlocks(gridCube);
   grid->accept(genBlocks);

   //boundary conditions definition 
   //boundary conditions adapters
   //////////////////////////////////////////////////////////////////////////////
   BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
   noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithm()));

   //boundary conditions visitor
   BoundaryConditionsBlockVisitor bcVisitor;
   bcVisitor.addBC(noSlipBCAdapter);
   //////////////////////////////////////////////////////////////////////////////////

   //set boundary conditions for blocks and create process decomposition for MPI
   D3Q27InteractorPtr cylinderInt(new D3Q27Interactor(cylinder, grid, noSlipBCAdapter, Interactor3D::INVERSESOLID));
   Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
   InteractorsHelper intHelper(grid, metisVisitor);
   intHelper.addInteractor(cylinderInt);
   intHelper.selectBlocks();

   //write data for visualization of block grid
   WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
   ppblocks->process(0);
   ppblocks.reset();

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

   //LBM kernel definition
   LBMKernelPtr kernel;
   kernel = LBMKernelPtr(new IncompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], IncompressibleCumulantLBMKernel::NORMAL));
   BCProcessorPtr bcProc(new BCProcessor());
   kernel->setBCProcessor(bcProc);

   //set forcing
   mu::Parser fctForcingX1;
   fctForcingX1.SetExpr("Fx1");
   fctForcingX1.DefineConst("Fx1", 9e-7);
   kernel->setWithForcing(true);
   kernel->setForcingX1(fctForcingX1);

   //create LBM kernel
   SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
   grid->accept(kernelVisitor);

   //set boundary conditions for nodes
   intHelper.setBC();
   grid->accept(bcVisitor);

   //initialization of distributions
   InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
   grid->accept(initVisitor);

   //set connectors
   InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
   SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
   grid->accept(setConnsVisitor);

   //domain decomposition for threads
   PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
   grid->accept(pqPartVisitor);

   //write data for visualization of boundary conditions
   UbSchedulerPtr geoSch(new UbScheduler(1));
   WriteBoundaryConditionsCoProcessorPtr ppgeo(
      new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), LBMUnitConverterPtr(new LBMUnitConverter()), comm));
   ppgeo->process(0);
   ppgeo.reset();
   
   if (myid == 0) UBLOG(logINFO, "Preprocess - end");

   //write data for visualization of macroscopic quantities
   UbSchedulerPtr visSch(new UbScheduler(outTime));
   WriteMacroscopicQuantitiesCoProcessor pp(grid, visSch, pathname, WbWriterVtkXmlASCII::getInstance(), LBMUnitConverterPtr(new LBMUnitConverter()), comm);

   //performance control
   UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
   NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

   //start simulation 
   const std::shared_ptr<ConcreteCalculatorFactory> calculatorFactory = std::make_shared<ConcreteCalculatorFactory>(visSch);
   CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, calculatorFactory, CalculatorType::HYBRID));
   if (myid == 0) UBLOG(logINFO, "Simulation-start");
   calculation->calculate();
   if (myid == 0) UBLOG(logINFO, "Simulation-end");
}



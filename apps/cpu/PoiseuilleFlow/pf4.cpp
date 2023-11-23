//#include "pf.h"
//#include "VirtualFluids.h"
//
//using namespace std;
//
////two plates flow with pressure drop
//void pf4()
//{
//   SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
//   int myid = comm->getProcessID();
//
//   //parameters
//   string          pathname = "d:/temp/pflow_plates_dp";
//   int             numOfThreads = 4;
//   int             blocknx[3] ={ 10,10,10 };
//   double          endTime = 100000;
//   double          outTime = 100;
//   double          availMem = 8e9;
//   double          deltax = 1;
//   double          rhoLBInflow = 0.001;
//   double          rhoLB = 0.0;
//   double          nuLB = 0.005;
//
//   //geometry definition
//
//   //simulation bounding box
//   double g_minX1 = 0.0;
//   double g_minX2 = -5.0;
//   double g_minX3 = -5.0;
//
//   double g_maxX1 = 30;
//   double g_maxX2 = 5;
//   double g_maxX3 = 5;
//
//   SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
//   if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
//
//   //walls
//   GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1-2.0*deltax, g_minX2-2.0*deltax, g_minX3-2.0*deltax, g_maxX1+2.0*deltax, g_minX2, g_maxX3+2.0*deltax));
//   if (myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());
//
//   GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1-2.0*deltax, g_maxX2, g_minX3-2.0*deltax, g_maxX1+2.0*deltax, g_maxX2+2.0*deltax, g_maxX3+2.0*deltax));
//   if (myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());
//
//   //inflow
//   GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-2.0*deltax, g_minX2-2.0*deltax, g_minX3-2.0*deltax, g_minX1, g_maxX2+2.0*deltax, g_maxX3+2.0*deltax));
//   if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());
//
//   //outflow
//   GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-2.0*deltax, g_minX3-2.0*deltax, g_maxX1+2.0*deltax, g_maxX2+2.0*deltax, g_maxX3+2.0*deltax));
//   if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());
//
//   if (myid == 0)
//   {
//      UBLOG(logINFO, "rhoLB = " << rhoLB);
//      UBLOG(logINFO, "nuLB = " << nuLB);
//      UBLOG(logINFO, "deltaX = " << deltax);
//      UBLOG(logINFO, "Preprocess - start");
//   }
//
//   //Grid definition
//   SPtr<Grid3D> grid(new Grid3D(comm));
//   grid->setDeltaX(deltax);
//   grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
//   grid->setPeriodicX1(false);
//   grid->setPeriodicX2(false);
//   grid->setPeriodicX3(true);
//
//   if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
//
//   //blocks generating
//   GenBlocksGridVisitor genBlocks(gridCube);
//   grid->accept(genBlocks);
//
//   //boundary conditions definition 
//   //boundary conditions adapters
//   //////////////////////////////////////////////////////////////////////////////
//   SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
//   noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
//   SPtr<BCAdapter> denInflowBCAdapter(new DensityBCAdapter(rhoLBInflow));
//   denInflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));
//   SPtr<BCAdapter> denOutflowBCAdapter(new DensityBCAdapter(rhoLB));
//   denOutflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));
//
//   //BC visitor
//   BoundaryConditionsBlockVisitor bcVisitor;
//   bcVisitor.addBC(noSlipBCAdapter);
//   bcVisitor.addBC(denInflowBCAdapter);
//   //////////////////////////////////////////////////////////////////////////////////
//
//   //set boundary conditions for blocks and create process decomposition for MPI
//   //walls
//   SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
//   SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
//   //inflow
//   SPtr<D3Q27Interactor> inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, denInflowBCAdapter, Interactor3D::SOLID));
//   //outflow
//   SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, denOutflowBCAdapter, Interactor3D::SOLID));
//
//   SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
//   InteractorsHelper intHelper(grid, metisVisitor);
//   intHelper.addInteractor(addWallYminInt);
//   intHelper.addInteractor(addWallYmaxInt);
//   intHelper.addInteractor(inflowInt);
//   intHelper.addInteractor(outflowInt);
//   intHelper.selectBlocks();
//
//   //write data for visualization of block grid
//   WriteBlocksSPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
//   ppblocks->process(0);
//   ppblocks.reset();
//
//   unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
//   int ghostLayer = 3;
//   unsigned long long numberOfNodesPerBlock = (unsigned long long)(blocknx[0])* (unsigned long long)(blocknx[1])* (unsigned long long)(blocknx[2]);
//   unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
//   unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
//   double needMemAll = double(numberOfNodesPerBlockWithGhostLayer*(27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
//   double needMem = needMemAll / double(comm->getNumberOfProcesses());
//
//   if (myid == 0)
//   {
//      UBLOG(logINFO, "Number of blocks = " << numberOfBlocks);
//      UBLOG(logINFO, "Number of nodes  = " << numberOfNodes);
//      int minInitLevel = grid->getCoarsestInitializedLevel();
//      int maxInitLevel = grid->getFinestInitializedLevel();
//      for (int level = minInitLevel; level <= maxInitLevel; level++)
//      {
//         int nobl = grid->getNumberOfBlocks(level);
//         UBLOG(logINFO, "Number of blocks for level " << level << " = " << nobl);
//         UBLOG(logINFO, "Number of nodes for level " << level << " = " << nobl*numberOfNodesPerBlock);
//      }
//      UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
//      UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
//      UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
//   }
//
//   //LBM kernel definition
//   SPtr<LBMKernel> kernel;
//   kernel = SPtr<LBMKernel>(new InK15CompressibleNavierStokes(blocknx[0], blocknx[1], blocknx[2], InK15CompressibleNavierStokes::NORMAL));
//   SPtr<BCProcessor> bcProc(new BCProcessor());
//   kernel->setBCProcessor(bcProc);
//
//   //create LBM kernel
//   SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
//   grid->accept(kernelVisitor);
//
//   //set boundary conditions for nodes
//   intHelper.setBC();
//   grid->accept(bcVisitor);
//
//   //initialization of distributions
//   InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
//   grid->accept(initVisitor);
//
//   //set connectors
//   InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolator());
//   SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
//   grid->accept(setConnsVisitor);
//
//   //domain decomposition for threads
//   PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
//   grid->accept(pqPartVisitor);
//
//   //write data for visualization of boundary conditions
//   SPtr<UbScheduler> geoSch(new UbScheduler(1));
//   WriteBoundaryConditionsSPtr<CoProcessor> ppgeo(
//      new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
//   ppgeo->process(0);
//   ppgeo.reset();
//   
//   if (myid == 0) UBLOG(logINFO, "Preprocess - end");
//
//   //write data for visualization of macroscopic quantities
//   SPtr<UbScheduler> visSch(new UbScheduler(outTime));
//   WriteMacroscopicQuantitiesCoProcessor pp(grid, visSch, pathname, WbWriterVtkXmlASCII::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm);
//
//   //performance control
//   SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
//   NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);
//
//   //start solver
//   CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
//   if (myid == 0) UBLOG(logINFO, "Simulation-start");
//   calculation->calculate();
//   if (myid == 0) UBLOG(logINFO, "Simulation-end");
//}
//
//

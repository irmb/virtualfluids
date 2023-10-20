#include "pf.h"
#include "VirtualFluids.h"
#include "CheckpointConverter.h"

using namespace std;

//pipe flow with forcing
void pf1()
{
    using namespace vf::lbm::dir;

   SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
   int myid = comm->getProcessID();

   //parameters
   string          pathOut = "d:/temp/test_dir_naming";  //"/gfs1/work/niikonst/pflow_pipe_forcing";
   int             numOfThreads = 1;
   int             blocknx[3] ={ 10,10,10 };
   real          endTime = 10;
   real          cpStart = 10;
   real          cpStep = 10;
   real          outTime = 10;
   real          availMem = 8e9;
   real          deltax = 1;
   real          rhoLB = 0.0;
   real          nuLB = 0.005;

   //geometry definition

   //simulation bounding box
   real g_minX1 = 0.0;
   real g_minX2 = -10.0;
   real g_minX3 = -10.0;

   real g_maxX1 = 50;
   real g_maxX2 = 10;
   real g_maxX3 = 10;

   //Sleep(15000);

   SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
   if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

   //cylinder
   SPtr<GbObject3D> cylinder(new GbCylinder3D(g_minX1 - 2.0*deltax, 0.0, 0.0, g_maxX1 + 2.0*deltax, 0.0, 0.0, g_maxX2));
   GbSystem3D::writeGeoObject(cylinder.get(), pathOut + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

   if (myid == 0)
   {
      UBLOG(logINFO, "rhoLB = " << rhoLB);
      UBLOG(logINFO, "nuLB = " << nuLB);
      UBLOG(logINFO, "deltaX = " << deltax);
      UBLOG(logINFO, "Preprocess - start");
   }

   //Grid definition
   SPtr<Grid3D> grid(new Grid3D(comm));
   grid->setDeltaX(deltax);
   grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
   grid->setPeriodicX1(true);
   grid->setPeriodicX2(false);
   grid->setPeriodicX3(false);

   //blocks generating
   GenBlocksGridVisitor genBlocks(gridCube);
   grid->accept(genBlocks);

   //boundary conditions definition 
   //boundary conditions adapters
   //////////////////////////////////////////////////////////////////////////////
   SPtr<BC> noSlipBC(new NoSlipBC());
   noSlipBC->setBCStrategy(SPtr<BCStrategy>(new NoSlipBCStrategy()));

   //boundary conditions visitor
   BoundaryConditionsBlockVisitor bcVisitor;
   //////////////////////////////////////////////////////////////////////////////////

   //set boundary conditions for blocks and create process decomposition for MPI
   SPtr<D3Q27Interactor> cylinderInt(new D3Q27Interactor(cylinder, grid, noSlipBC, Interactor3D::INVERSESOLID));
   SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, DIR_00M));
   InteractorsHelper intHelper(grid, metisVisitor);
   intHelper.addInteractor(cylinderInt);
   intHelper.selectBlocks();

   //write data for visualization of block grid
   SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm));
   ppblocks->update(0);
   //ppblocks.reset();

   unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
   int ghostLayer = 3;
   unsigned long long numberOfNodesPerBlock = (unsigned long long)(blocknx[0])* (unsigned long long)(blocknx[1])* (unsigned long long)(blocknx[2]);
   unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
   unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
   real needMemAll = real(numberOfNodesPerBlockWithGhostLayer*(27 * sizeof(real) + sizeof(int) + sizeof(float) * 4));
   real needMem = needMemAll / real(comm->getNumberOfProcesses());

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
   SPtr<LBMKernel> kernel;
   kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel());
   SPtr<BCSet> bcProc(new BCSet());
   kernel->setBCSet(bcProc);

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
   InitDistributionsBlockVisitor initVisitor;
   grid->accept(initVisitor);

   //set connectors
   OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
   grid->accept(setConnsVisitor);

   SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
   //SPtr<MPIIORestartSimulationObserver> restartSimulationObserver(new MPIIORestartSimulationObserver(grid, mSch, pathOut, comm));
   //restartSimulationObserver->setLBMKernel(kernel);
   //restartSimulationObserver->setBCSet(bcProc);

   /*SPtr<MPIIOMigrationSimulationObserver> migSimulationObserver(new MPIIOMigrationSimulationObserver(grid, mSch, pathOut + "/mig", comm));
   migSimulationObserver->setLBMKernel(kernel);
   migSimulationObserver->setBCSet(bcProc);*/

   //SPtr<MPIIOMigrationBESimulationObserver> migSimulationObserver(new MPIIOMigrationBESimulationObserver(grid, mSch, pathOut + "/mig", comm));
   //migSimulationObserver->setLBMKernel(kernel);
   //migSimulationObserver->setBCSet(bcProc);
   //migSimulationObserver->setNu(nuLB);

   //SPtr<UtilConvertor> convertProcessor(new UtilConvertor(grid, pathOut, comm));
   //convertProcessor->convert(300, 4);
   //return;

   //write data for visualization of boundary conditions
   {
      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      WriteBoundaryConditionsSimulationObserver ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), /*SPtr<LBMUnitConverter>(new LBMUnitConverter()),*/ comm);
      ppgeo.update(0);
   }
   
   if (myid == 0) UBLOG(logINFO, "Preprocess - end");

   //grid=SPtr<Grid3D>(new Grid3D(comm));
   //restartSimulationObserver->restart(200);
   //SPtr<MPIIOMigrationBESimulationObserver> migSimulationObserver(new MPIIOMigrationBESimulationObserver(grid, mSch, metisVisitor, pathOut + "/mig", comm));
   //migSimulationObserver->setLBMKernel(kernel);
   //migSimulationObserver->setBCSet(bcProc);
   //migSimulationObserver->setNu(nuLB);
   //migSimulationObserver->restart(10);

   ppblocks->update(1);

   //write data for visualization of macroscopic quantities
   SPtr<UbScheduler> visSch(new UbScheduler(outTime));
   SPtr<WriteMacroscopicQuantitiesSimulationObserver> writeMQSimulationObserver(new WriteMacroscopicQuantitiesSimulationObserver(grid, visSch, pathOut, 
      WbWriterVtkXmlASCII::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

   //performance control
   SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
   SPtr<NUPSCounterSimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

   //start simulation 
   //omp_set_num_threads(numOfThreads);
   SPtr<UbScheduler> stepGhostLayer(new UbScheduler(outTime));
   SPtr<Simulation> simulation(new Simulation(grid, stepGhostLayer, endTime));
   simulation->addSimulationObserver(npr);
   simulation->addSimulationObserver(writeMQSimulationObserver);
   //simulation->addSimulationObserver(migSimulationObserver);
   //simulation->addSimulationObserver(restartSimulationObserver);

   if (myid == 0) UBLOG(logINFO, "Simulation-start");
   simulation->run();
   if (myid == 0) UBLOG(logINFO, "Simulation-end");
   
   ppblocks->update(10);
}



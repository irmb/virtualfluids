#include <iostream>
#include <string>

#include "PointerDefinitions.h"

#include <iostream>
#include <string>
#include <memory>
#include <array>

#include "VirtualFluids.h"
#include <MuParser/include/muParser.h>
#include "ForceCalculator.h"


#include <MovableObjectInteractor.h>
#include <DemCalculator.h>
#include <DemCoProcessor.h>
#include <PePartitioningGridVisitor.h>

#include <PePhysicsEngineMaterialAdapter.h>
#include <PePhysicsEngineGeometryAdapter.h>
#include <PePhysicsEngineSolverAdapter.h>

#include <VelocityBcReconstructor.h>
#include <EquilibriumReconstructor.h>
#include <ExtrapolationReconstructor.h>

#include <DummyPhysicsEngineSolverAdapter.h>
#include <DummyPhysicsEngineMaterialAdapter.h>
#include <DummyPhysicsEngineGeometryAdapter.h>

#include "CreateGeoObjectsCoProcessor.h"

using namespace std;

std::shared_ptr<DemCoProcessor> makePeCoProcessor(SPtr<Grid3D> grid, SPtr<Communicator> comm, const SPtr<UbScheduler> peScheduler, const std::shared_ptr<LBMUnitConverter> lbmUnitConverter)
{
   double peRelaxtion = 0.7;
   int maxpeIterations = 10000;
   //Beschleunigung g
   //double g = 2.0 * 9.81 * lbmUnitConverter->getFactorAccWToLb();
   Vector3D globalLinearAcc(0.0, 0.0, 0.0);

   std::shared_ptr<PePhysicsEngineMaterialAdapter> planeMaterial = std::make_shared<PePhysicsEngineMaterialAdapter>("granular", 1.0, 0, 0.1 / 2, 0.1 / 2, 0.5, 1, 1, 0, 0);

   const int gridNx = val<1>(grid->getBlockNX()) * grid->getNX1();
   const int gridNy = val<2>(grid->getBlockNX()) * grid->getNX2();
   const int gridNz = val<3>(grid->getBlockNX()) * grid->getNX3();

   UbTupleInt3 simulationDomain(gridNx, gridNy, gridNz);
   UbTupleInt3 numberOfBlocks(grid->getNX1(), grid->getNX2(), grid->getNX3());
   UbTupleBool3 isPeriodic(grid->isPeriodicX1(), grid->isPeriodicX2(), grid->isPeriodicX3());

   std::shared_ptr<PeParameter> peParamter = std::make_shared<PeParameter>(peRelaxtion, maxpeIterations, globalLinearAcc, planeMaterial, simulationDomain, numberOfBlocks, isPeriodic);
   std::shared_ptr<PhysicsEngineSolverAdapter> peSolver = std::make_shared<PePhysicsEngineSolverAdapter>(peParamter);

   const std::shared_ptr<ForceCalculator> forceCalculator = std::make_shared<ForceCalculator>(comm);

   return std::make_shared<DemCoProcessor>(grid, peScheduler, comm, forceCalculator, peSolver);
}

//pipe flow with forcing
void pf1()
{
   SPtr<Communicator> comm = MPICommunicator::getInstance();
   int myid = comm->getProcessID();

   //parameters
   string          pathname = "d:/temp/thermoplast";
   int             numOfThreads = 1;
   int             blocknx[3] ={ 16,16,16 };
   double          endTime = 1000000;
   double          outTime = 10;
   double          availMem = 8e9;
   double          deltax = 1;
   double          rhoLB = 0.0;
   double          nuLB = 0.005;
   double          uLB =  0.1;

   //geometry definition

   //simulation bounding box
   double g_minX1 = 0.0;
   double g_minX2 = 0.0;
   double g_minX3 = 0.0;

   double g_maxX1 = 32;
   double g_maxX2 = 32;
   double g_maxX3 = 32;

   double blockLength = 10;

   SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
   if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

   //box
   SPtr<GbObject3D> box(new GbCuboid3D(g_minX1-blockLength, g_minX2, g_minX3, g_maxX1+blockLength, g_maxX2, g_maxX3));
   GbSystem3D::writeGeoObject(box.get(), pathname + "/geo/box", WbWriterVtkXmlBinary::getInstance());

   //sphere
   double radius = 5;
   SPtr<GbObject3D> sphere1(new GbSphere3D(10, g_maxX2/2.0, g_maxX3/2.0, radius));
   if (myid == 0) GbSystem3D::writeGeoObject(sphere1.get(), pathname + "/geo/sphere", WbWriterVtkXmlASCII::getInstance());

   if (myid == 0)
   {
      UBLOG(logINFO, "uLB = " << uLB);
      UBLOG(logINFO, "rhoLB = " << rhoLB);
      UBLOG(logINFO, "nuLB = " << nuLB);
      UBLOG(logINFO, "deltaX = " << deltax);
      UBLOG(logINFO, "Preprocess - start");
   }

   //Grid definition
   SPtr<Grid3D> grid(new Grid3D(comm));
   grid->setDeltaX(deltax);
   grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
   grid->setPeriodicX1(false);
   grid->setPeriodicX2(false);
   grid->setPeriodicX3(false);

   //blocks generating
   GenBlocksGridVisitor genBlocks(gridCube);
   grid->accept(genBlocks);

   //boundary conditions definition 
   //boundary conditions adapters
   //////////////////////////////////////////////////////////////////////////////
   SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
   noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

   mu::Parser fct;
   fct.SetExpr("U");
   fct.DefineConst("U", uLB);
   SPtr<BCAdapter> inflowAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
   inflowAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));

   SPtr<BCAdapter> outflowAdapter(new DensityBCAdapter(rhoLB));
   outflowAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

   //sphere
   mu::Parser fct2;
   fct2.SetExpr("U");
   fct2.DefineConst("U", 0.0);
   SPtr<BCAdapter> velocityBcParticleAdapter(new VelocityBCAdapter(true, false, false, fct2, 0, BCFunction::INFCONST));
   velocityBcParticleAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

   //inflow
   GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_minX1, g_maxX2 + blockLength, g_maxX3 + blockLength));
   if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

   //outflow
   GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
   if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

   //boundary conditions visitor
   SPtr<BoundaryConditionsBlockVisitor> bcVisitor(new BoundaryConditionsBlockVisitor());
   bcVisitor->addBC(noSlipBCAdapter);
   bcVisitor->addBC(inflowAdapter);
   bcVisitor->addBC(outflowAdapter);
   bcVisitor->addBC(velocityBcParticleAdapter);
   //////////////////////////////////////////////////////////////////////////////////

   //set boundary conditions for blocks and create process decomposition for MPI
   SPtr<D3Q27Interactor> boxInt(new D3Q27Interactor(box, grid, noSlipBCAdapter, Interactor3D::INVERSESOLID));

   //sphere bc object
   SPtr<MovableObjectInteractor> sphereInt1;
   const std::shared_ptr<Reconstructor> velocityReconstructor = std::make_shared<VelocityBcReconstructor>();
   const std::shared_ptr<Reconstructor> extrapolationReconstructor = std::make_shared<ExtrapolationReconstructor>(velocityReconstructor);

   sphereInt1 = SPtr<MovableObjectInteractor>(new MovableObjectInteractor(sphere1, grid, velocityBcParticleAdapter, Interactor3D::SOLID, extrapolationReconstructor, State::UNPIN));

   //inflow
   SPtr<D3Q27Interactor> inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, inflowAdapter, Interactor3D::SOLID));

   //outflow
   SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowAdapter, Interactor3D::SOLID));

   //PE
   double refLengthLb = radius*2.0;
   double refLengthWorld = refLengthLb * deltax;
   double refVelWorld = 0.1;
   //const std::shared_ptr<LBMUnitConverter> lbmUnitConverter = std::make_shared<LBMUnitConverter>(refLengthWorld, LBMUnitConverter::WORLD_MATERIAL::AIR_20C, refLengthLb);
   const SPtr<LBMUnitConverter> lbmUnitConverter(new LBMUnitConverter(refLengthWorld, LBMUnitConverter::WORLD_MATERIAL::OIL, refLengthLb));
   if (myid == 0) std::cout << lbmUnitConverter->toString() << std::endl;
   double rhoSphere = 915 * lbmUnitConverter->getFactorDensityWToLb(); // 2 // kg/m^3
   if (myid == 0) UBLOG(logINFO, "rhoSphere = "<<rhoSphere);
   SPtr<PhysicsEngineMaterialAdapter> sphereMaterial(new PePhysicsEngineMaterialAdapter("Polyethylen", rhoSphere, 0, 0.15, 0.1, 0.45, 0.5, 1, 0, 0));
   const int timestep = 2;
   const SPtr<UbScheduler> peScheduler(new UbScheduler(timestep));
   SPtr<DemCoProcessor> demCoProcessor = makePeCoProcessor(grid, comm, peScheduler, lbmUnitConverter);
   //demCoProcessor->addInteractor(sphereInt1, sphereMaterial, Vector3D(0.0, 0.0, 0.0));
   //demCoProcessor->distributeIDs();
   demCoProcessor->setBlockVisitor(bcVisitor);
   double g = 9.81 * lbmUnitConverter->getFactorAccWToLb();
   double rhoFluid = lbmUnitConverter->getFactorDensityWToLb() * 830.0; // 1 // kg/m^3
   //////////////////////////////////////////////////////////////////////////

   SPtr<Grid3DVisitor> peVisitor(new PePartitioningGridVisitor(comm, demCoProcessor));
   InteractorsHelper intHelper(grid, peVisitor);
   intHelper.addInteractor(boxInt);
   //intHelper.addInteractor(sphereInt1);
   intHelper.addInteractor(inflowInt);
   intHelper.addInteractor(outflowInt);

   intHelper.selectBlocks();

   //write data for visualization of block grid
   SPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
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
   SPtr<LBMKernel> kernel;
   kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel());
   SPtr<BCProcessor> bcProc(new BCProcessor());
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
   grid->accept(*bcVisitor.get());

   //initialization of distributions
   InitDistributionsBlockVisitor initVisitor;
   grid->accept(initVisitor);

   //set connectors
   InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
   SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
   grid->accept(setConnsVisitor);

   //write data for visualization of boundary conditions
   {
      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm);
      ppgeo.process(0);
   }
   
   if (myid == 0) UBLOG(logINFO, "Preprocess - end");

   //write data for visualization of macroscopic quantities
   SPtr<UbScheduler> visSch(new UbScheduler(outTime,1));
   SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, pathname, 
      WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

   SPtr<WriteBoundaryConditionsCoProcessor> writeBCCoProcessor(new WriteBoundaryConditionsCoProcessor(grid, visSch, pathname,
      WbWriterVtkXmlBinary::getInstance(), comm));


   //performance control
   SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
   SPtr<NUPSCounterCoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

   //////////////////////////////////////////////////////////////////////////

   SPtr<UbScheduler> sphereScheduler(new UbScheduler(1,1,2));
   SPtr<CreateGeoObjectsCoProcessor> createSphereCoProcessor(new CreateGeoObjectsCoProcessor(grid,sphereScheduler,demCoProcessor,sphere1,sphereMaterial));


   //start simulation 
   //omp_set_num_threads(numOfThreads);
   SPtr<UbScheduler> stepGhostLayer(new UbScheduler(outTime));
   SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
   //SPtr<Calculator> calculator(new DemCalculator(grid, stepGhostLayer, endTime));
   calculator->addCoProcessor(npr);
   
   calculator->addCoProcessor(createSphereCoProcessor);
   calculator->addCoProcessor(demCoProcessor);
   
   
   calculator->addCoProcessor(writeBCCoProcessor);
   calculator->addCoProcessor(writeMQCoProcessor);

   if (myid == 0) UBLOG(logINFO, "Simulation-start");
   calculator->calculate();
   if (myid == 0) UBLOG(logINFO, "Simulation-end");
}


//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   try
   {
      //Sleep(30000);
      walberla::Environment env(argc, argv);

      pf1();
      //peFlow(std::string(argv[1]));
      return 0;
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

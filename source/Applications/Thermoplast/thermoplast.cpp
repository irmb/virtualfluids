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
#include <WriteDemObjectsCoProcessor.h>

#include "CreateDemObjectsCoProcessor.h"

using namespace std;

//simulation bounding box
double g_minX1 = 0;
double g_minX2 = 0;
double g_minX3 = 0;

double g_maxX1 = 0;
double g_maxX2 = 0;
double g_maxX3 = 0;

vector<double> peMinOffset;
vector<double> peMaxOffset;

string          pathOut = "d:/temp/thermoplastCluster";
string          pathGeo = "d:/Projects/ThermoPlast/Geometrie";

std::shared_ptr<DemCoProcessor> makePeCoProcessor(SPtr<Grid3D> grid, SPtr<Communicator> comm, const SPtr<UbScheduler> peScheduler, const std::shared_ptr<LBMUnitConverter> lbmUnitConverter,  int maxpeIterations)
{
   double peRelaxtion = 0.7;
   //int maxpeIterations = 10000;
   //Beschleunigung g
   double g = 9.81 * lbmUnitConverter->getFactorAccWToLb();
   Vector3D globalLinearAcc(0.0, -g, 0.0);

   std::shared_ptr<PePhysicsEngineMaterialAdapter> planeMaterial = std::make_shared<PePhysicsEngineMaterialAdapter>("granular", 1.0, 0, 0.1 / 2, 0.1 / 2, 0.5, 1, 1, 0, 0);

   const int gridNx = val<1>(grid->getBlockNX()) * grid->getNX1();
   const int gridNy = val<2>(grid->getBlockNX()) * grid->getNX2();
   const int gridNz = val<3>(grid->getBlockNX()) * grid->getNX3();

   //UbTupleInt3 simulationDomain(gridNx, gridNy, gridNz);
   //std::array<double, 6> simulationDomain = {1, 1, 1, 30, 30, 30};
   std::array<double, 6> simulationDomain = {g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3};
   UbTupleInt3 numberOfBlocks(grid->getNX1(), grid->getNX2(), grid->getNX3());
   //UbTupleInt3 numberOfBlocks((simulationDomain[3]-simulationDomain[0])/val<1>(grid->getBlockNX()), (simulationDomain[4]-simulationDomain[1])/val<2>(grid->getBlockNX()), (simulationDomain[5]-simulationDomain[2])/val<3>(grid->getBlockNX()));
   UbTupleBool3 isPeriodic(grid->isPeriodicX1(), grid->isPeriodicX2(), grid->isPeriodicX3());
   Vector3D minOffset(peMinOffset[0], peMinOffset[1], peMinOffset[2]);
   Vector3D maxOffset(peMaxOffset[0], peMaxOffset[1], peMaxOffset[2]);
   
   SPtr<GbObject3D> boxPE(new GbCuboid3D(simulationDomain[0]+minOffset[0], simulationDomain[1]+minOffset[1], simulationDomain[2]+minOffset[2], simulationDomain[3]+maxOffset[0], simulationDomain[4]+maxOffset[1], simulationDomain[5]+maxOffset[2]));
   GbSystem3D::writeGeoObject(boxPE.get(), pathOut + "/geo/boxPE", WbWriterVtkXmlBinary::getInstance());

   std::shared_ptr<PeParameter> peParamter = std::make_shared<PeParameter>(peRelaxtion, maxpeIterations, globalLinearAcc, 
      planeMaterial, simulationDomain, numberOfBlocks, isPeriodic, minOffset, maxOffset);
   std::shared_ptr<PhysicsEngineSolverAdapter> peSolver = std::make_shared<PePhysicsEngineSolverAdapter>(peParamter);

   const std::shared_ptr<ForceCalculator> forceCalculator = std::make_shared<ForceCalculator>(comm);

   return std::make_shared<DemCoProcessor>(grid, peScheduler, comm, forceCalculator, peSolver);
}

void createGridObjects()
{

}

//pipe flow with forcing
void thermoplast(string configname)
{
   try
   {
   SPtr<Communicator> comm = MPICommunicator::getInstance();
   int myid = comm->getProcessID();

   ConfigurationFile   config;
   config.load(configname);

   vector<int>     blocknx = config.getVector<int>("blocknx");
   vector<double>  boundingBox = config.getVector<double>("boundingBox");
                   peMinOffset = config.getVector<double>("peMinOffset");
                   peMaxOffset = config.getVector<double>("peMaxOffset");
   int             endTime = config.getValue<int>("endTime");
   double          outTime = config.getValue<double>("outTime");
   double          availMem = config.getValue<double>("availMem");
   double          uLB = config.getValue<double>("uLB");
   string          pathOut = config.getValue<string>("pathOut");
   string          pathGeo = config.getValue<string>("pathGeo");
   string          michel = config.getValue<string>("michel");
   string          plexiglas = config.getValue<string>("plexiglas");
   double          sphereTime = config.getValue<double>("sphereTime");

   //parameters
   //string          pathOut = "d:/temp/thermoplast3";
   //string          pathGeo = "d:/Projects/ThermoPlast/Geometrie";
   int             numOfThreads = 1;
   //int             blocknx[3] ={ 10,10,10 };
   //double          endTime = 1000000;
   //double          outTime = 300;
   //double          availMem = 8e9;
   double          deltax = 1;
   double          rhoLB = 0.0;
   //double          uLB =  0.1;
   double          radius = 5;
   double          Re = 300;
   double          nuLB = (uLB*2.0*radius)/Re;

             

   //geometry definition

   //simulation bounding box
   g_minX1 = boundingBox[0];
   g_minX2 = boundingBox[1];
   g_minX3 = boundingBox[2];

   g_maxX1 = boundingBox[3];
   g_maxX2 = boundingBox[4];
   g_maxX3 = boundingBox[5];

   double blockLength = blocknx[0]*deltax;

   SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
   if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

   //box
   SPtr<GbObject3D> box(new GbCuboid3D(g_minX1-blockLength, g_minX2, g_minX3, g_maxX1+blockLength, g_maxX2, g_maxX3));
   GbSystem3D::writeGeoObject(box.get(), pathOut + "/geo/box", WbWriterVtkXmlBinary::getInstance());

   //michel
   if (myid==0) UBLOG(logINFO, "Read michelGeo:start");
   SPtr<GbTriFaceMesh3D> michelGeo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+michel, "michelGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
   if (myid==0) UBLOG(logINFO, "Read michelGeo:end");
   if (myid==0) GbSystem3D::writeGeoObject(michelGeo.get(), pathOut+"/geo/michelGeo", WbWriterVtkXmlBinary::getInstance());

   //plexiglas
   if (myid==0) UBLOG(logINFO, "Read plexiglasGeo:start");
   SPtr<GbTriFaceMesh3D> plexiglasGeo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+plexiglas, "plexiglasGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
   if (myid==0) UBLOG(logINFO, "Read plexiglasGeo:end");
   if (myid==0) GbSystem3D::writeGeoObject(plexiglasGeo.get(), pathOut+"/geo/plexiglasGeo", WbWriterVtkXmlBinary::getInstance());

   if (myid == 0)
   {
      UBLOG(logINFO, "Parameters:");
      UBLOG(logINFO, "* uLB    = " << uLB);
      UBLOG(logINFO, "* rhoLB  = " << rhoLB);
      UBLOG(logINFO, "* nuLB   = " << nuLB);
      UBLOG(logINFO, "* deltaX = " << deltax);
      UBLOG(logINFO, "* radius = " << radius);
      UBLOG(logINFO, "* Re     = " << Re);
      UBLOG(logINFO, "* number of threads   = "<<numOfThreads);
      UBLOG(logINFO, "* number of processes = "<<comm->getNumberOfProcesses());
      UBLOG(logINFO, "* path = "<<pathOut);
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
   //outflowAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

   //sphere
   mu::Parser fct2;
   fct2.SetExpr("U");
   fct2.DefineConst("U", 0.0);
   SPtr<BCAdapter> velocityBcParticleAdapter(new VelocityBCAdapter(true, false, false, fct2, 0, BCFunction::INFCONST));
   velocityBcParticleAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

   //inflow
   GbCuboid3DPtr geoInflow1(new GbCuboid3D(g_minX1-blockLength, g_maxX2-120.0, g_minX3+190.0, g_minX1+1, g_maxX2+20.0, g_minX3+320.0));   if (myid == 0) GbSystem3D::writeGeoObject(geoInflow1.get(), pathOut + "/geo/geoInflow1", WbWriterVtkXmlASCII::getInstance());

   GbCuboid3DPtr geoInflow2(new GbCuboid3D(g_minX1-blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_minX1, g_maxX2 + blockLength, g_maxX3 + blockLength));
   if (myid == 0) GbSystem3D::writeGeoObject(geoInflow2.get(), pathOut + "/geo/geoInflow2", WbWriterVtkXmlASCII::getInstance());

   GbCuboid3DPtr geoInflow3(new GbCuboid3D(g_minX1-blockLength, g_minX2-radius, g_maxX3-4.0*radius-1, g_minX1+1, g_maxX2+radius, g_maxX3-4.0*radius));
   if (myid == 0) GbSystem3D::writeGeoObject(geoInflow3.get(), pathOut + "/geo/geoInflow3", WbWriterVtkXmlASCII::getInstance());

   GbCuboid3DPtr geoInflow4(new GbCuboid3D(g_minX1-blockLength, g_minX2+4.0*radius, g_maxX3-4.0*radius-1.0, g_minX1+1, g_minX2+4.0*radius+1.0, g_maxX3+radius));
   if (myid == 0) GbSystem3D::writeGeoObject(geoInflow4.get(), pathOut + "/geo/geoInflow4", WbWriterVtkXmlASCII::getInstance());

   //outflow
   GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
   if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

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
   //SPtr<MovableObjectInteractor> sphereInt1;
   //const std::shared_ptr<Reconstructor> velocityReconstructor = std::make_shared<VelocityBcReconstructor>();
   //const std::shared_ptr<Reconstructor> extrapolationReconstructor = std::make_shared<ExtrapolationReconstructor>(velocityReconstructor);
   //sphereInt1 = SPtr<MovableObjectInteractor>(new MovableObjectInteractor(sphere1, grid, velocityBcParticleAdapter, Interactor3D::SOLID, extrapolationReconstructor, State::UNPIN));

   //inflow
   SPtr<D3Q27Interactor> inflowInt1 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow1, grid, inflowAdapter, Interactor3D::SOLID));
   SPtr<D3Q27Interactor> inflowInt2 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow2, grid, outflowAdapter, Interactor3D::SOLID));

   SPtr<D3Q27Interactor> inflowInt3 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow3, grid, noSlipBCAdapter, Interactor3D::SOLID));
   SPtr<D3Q27Interactor> inflowInt4 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow4, grid, noSlipBCAdapter, Interactor3D::SOLID));

   //outflow
   SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowAdapter, Interactor3D::SOLID));

   //michel
   SPtr<Interactor3D> michelInt = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(michelGeo, grid, noSlipBCAdapter, Interactor3D::SOLID));

   //plexiglas
   SPtr<Interactor3D> plexiglasInt = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(plexiglasGeo, grid, noSlipBCAdapter, Interactor3D::SOLID));

   //PE
   double refLengthLb = radius*2.0;
   double refLengthWorld = refLengthLb * deltax;
   const std::shared_ptr<LBMUnitConverter> lbmUnitConverter = std::make_shared<LBMUnitConverter>(refLengthWorld, LBMUnitConverter::WORLD_MATERIAL::AIR_20C, refLengthLb);
   if (myid == 0) std::cout << lbmUnitConverter->toString() << std::endl;
   double rhoSphere = 915 * lbmUnitConverter->getFactorDensityWToLb();  // kg/m^3
   if (myid == 0) UBLOG(logINFO, "rhoSphere = "<<rhoSphere);
   SPtr<PhysicsEngineMaterialAdapter> sphereMaterial(new PePhysicsEngineMaterialAdapter("Polypropylen", rhoSphere, 0, 0.15, 0.1, 0.45, 0.5, 1, 0, 0));
   const int timestep = 2;
   const SPtr<UbScheduler> peScheduler(new UbScheduler(timestep));
   int maxpeIterations = 10;//endTime/2;
   SPtr<DemCoProcessor> demCoProcessor = makePeCoProcessor(grid, comm, peScheduler, lbmUnitConverter, maxpeIterations);
   //demCoProcessor->addInteractor(sphereInt1, sphereMaterial, Vector3D(0.0, 0.0, 0.0));
   //demCoProcessor->distributeIDs();
   demCoProcessor->setBlockVisitor(bcVisitor);
   //double g = 9.81 * lbmUnitConverter->getFactorAccWToLb();
   //double rhoFluid = lbmUnitConverter->getFactorDensityWToLb() * 830.0; // 1 // kg/m^3
   //////////////////////////////////////////////////////////////////////////
   SPtr<Grid3DVisitor> peVisitor(new PePartitioningGridVisitor(comm, demCoProcessor));
   //grid->accept(peVisitor);

   //SetSolidBlocksBlockVisitor solidBlocksVisitor1(michelInt);
   //grid->accept(solidBlocksVisitor1);
   //SetSolidBlocksBlockVisitor solidBlocksVisitor2(plexiglasInt);
   //grid->accept(solidBlocksVisitor2);

   //SetBcBlocksBlockVisitor bcBlocksVisitor1(michelInt);
   //grid->accept(bcBlocksVisitor1);
   //SetBcBlocksBlockVisitor bcBlocksVisitor2(plexiglasInt);
   //grid->accept(bcBlocksVisitor2);
   
   InteractorsHelper intHelper(grid, peVisitor);
   intHelper.addInteractor(boxInt);
   //intHelper.addInteractor(sphereInt1);
   //intHelper.addInteractor(inflowInt2);
   intHelper.addInteractor(michelInt);
   intHelper.addInteractor(plexiglasInt);
   intHelper.addInteractor(inflowInt1);
   intHelper.addInteractor(outflowInt);
   intHelper.addInteractor(inflowInt2);
   //intHelper.addInteractor(inflowInt3);
   //intHelper.addInteractor(inflowInt4);
   
   //intHelper.addInteractor(plexiglasInt);

   intHelper.selectBlocks();

   //write data for visualization of block grid
   SPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm));
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
   //mu::Parser fctForcingX1;
   //fctForcingX1.SetExpr("Fx1");
   //fctForcingX1.DefineConst("Fx1", 9e-7);
   //kernel->setWithForcing(true);
   //kernel->setForcingX1(fctForcingX1);

   //create LBM kernel
   SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
   grid->accept(kernelVisitor);

   //set boundary conditions for nodes
   //michelInt->initInteractor();
   //plexiglasInt->initInteractor();

   intHelper.setBC();
   grid->accept(*bcVisitor.get());

   //initialization of distributions
   InitDistributionsBlockVisitor initVisitor;
   initVisitor.setVx1(uLB);
   grid->accept(initVisitor);

   //set connectors
   InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
   SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
   grid->accept(setConnsVisitor);

   //write data for visualization of boundary conditions
   {
      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
      ppgeo.process(0);

      WriteMacroscopicQuantitiesCoProcessor ppInit(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm);
      ppInit.process(0);
   }
   
   if (myid == 0) UBLOG(logINFO, "Preprocess - end");

   //write data for visualization of macroscopic quantities
   SPtr<UbScheduler> visSch(new UbScheduler(outTime));
   SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, pathOut, 
      WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

   SPtr<WriteBoundaryConditionsCoProcessor> writeBCCoProcessor(new WriteBoundaryConditionsCoProcessor(grid, visSch, pathOut,
      WbWriterVtkXmlBinary::getInstance(), comm));

   SPtr<WriteDemObjectsCoProcessor> writeDemObjectsCoProcessor(new WriteDemObjectsCoProcessor(grid, visSch, pathOut, SPtr<WbWriter>(WbWriterVtkXmlBinary::getInstance()), demCoProcessor, comm));
   writeDemObjectsCoProcessor->process(0);

   //performance control
   SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
   SPtr<NUPSCounterCoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

   //////////////////////////////////////////////////////////////////////////
   //generating spheres 
   double d = 2.0*radius;   SPtr<UbScheduler> sphereScheduler(new UbScheduler(sphereTime));   Vector3D origin(g_minX1+60.0+4.0+radius, geoInflow1->getX2Minimum()+4.0*d,  geoInflow1->getX3Minimum()+2.0*d);
   //std::array<double, 6> AABB={origin[0]-radius-1,origin[1]-radius,origin[2]-radius-1,origin[0]+radius+1, origin[1]+2.0*d+radius+1, origin[2]+2.0*d+radius+1};
   //SPtr<GbObject3D> boxAABB(new GbCuboid3D(AABB[0],AABB[1],AABB[2],AABB[3],AABB[4],AABB[5]));
   //GbSystem3D::writeGeoObject(boxAABB.get(), pathOut + "/geo/boxAABB", WbWriterVtkXmlBinary::getInstance());
   SPtr<CreateDemObjectsCoProcessor> createSphereCoProcessor(new CreateDemObjectsCoProcessor(grid,sphereScheduler,demCoProcessor,sphereMaterial,Vector3D(uLB, 0.0, 0.0)));
   //spheres
   for (int x3 = 0; x3 < 6; x3++)      for (int x2 = 0; x2 < 4; x2++)         for (int x1 = 0; x1 < 1; x1++)
         {
            //SPtr<GbObject3D> sphere(new GbSphere3D(origin[0]+x1*d, origin[1]+x2*2.0*d, origin[2]+x3*2.0*d, radius));
            SPtr<GbObject3D> sphere(new GbSphere3D(origin[0]+x1*d, origin[1]+x2*1.5*d, origin[2]+x3*1.5*d, radius));
            if (myid == 0) GbSystem3D::writeGeoObject(sphere.get(), pathOut + "/geo/sphere"+UbSystem::toString(x1)+UbSystem::toString(x2)+UbSystem::toString(x3), WbWriterVtkXmlASCII::getInstance());
            createSphereCoProcessor->addGeoObject(sphere);
         }

   createSphereCoProcessor->process(0);

   //start simulation 
   omp_set_num_threads(numOfThreads);
   SPtr<UbScheduler> stepGhostLayer(visSch);
   SPtr<Calculator> calculator(new BasicCalculator(grid, peScheduler, endTime));
   calculator->addCoProcessor(npr);
   
   calculator->addCoProcessor(createSphereCoProcessor);
   calculator->addCoProcessor(demCoProcessor);
      
   calculator->addCoProcessor(writeBCCoProcessor);
   calculator->addCoProcessor(writeDemObjectsCoProcessor);
   calculator->addCoProcessor(writeMQCoProcessor);

   if (myid == 0) UBLOG(logINFO, "Simulation-start");
   calculator->calculate();
   if (myid == 0) UBLOG(logINFO, "Simulation-end");

      }
      catch (std::exception& e)      {         UBLOG(logERROR, e.what());      }      catch (std::string& s)      {         UBLOG(logERROR, s);      }      catch (...)      {         UBLOG(logERROR, "unknown exception");      }
}


//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   //try
   //{
      //Sleep(30000);
      walberla::Environment env(argc, argv);

      if (argv!=NULL)
      {
         if (argv[1]!=NULL)
         {
            thermoplast(string(argv[1]));
         }
         else
         {
            cout<<"Configuration file must be set!: "<<argv[0]<<" <config file>"<<endl<<std::flush;
         }
      }
      return 0;
   //}
   //catch (std::exception& e)   //{   //   UBLOG(logERROR, e.what());   //}   //catch (std::string& s)   //{   //   UBLOG(logERROR, s);   //}   //catch (...)   //{   //   UBLOG(logERROR, "unknown exception");   //}
}

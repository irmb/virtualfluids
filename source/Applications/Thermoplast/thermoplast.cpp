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
#include "PeLoadBalancerAdapter.h"

#include <VelocityBcReconstructor.h>
#include <EquilibriumReconstructor.h>
#include <ExtrapolationReconstructor.h>

#include <DummyPhysicsEngineSolverAdapter.h>
#include <DummyPhysicsEngineMaterialAdapter.h>
#include <DummyPhysicsEngineGeometryAdapter.h>
#include <WriteDemObjectsCoProcessor.h>

#include "CreateDemObjectsCoProcessor.h"
#include "RestartDemObjectsCoProcessor.h"

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

string          pathOut;// = "d:/temp/thermoplastCluster";
string          pathGeo;// = "d:/Projects/ThermoPlast/Geometrie";

std::shared_ptr<DemCoProcessor> makePeCoProcessor(SPtr<Grid3D> grid, SPtr<Communicator> comm, const SPtr<UbScheduler> peScheduler, const std::shared_ptr<LBMUnitConverter> lbmUnitConverter, int maxpeIterations)
{
   double peRelaxtion = 0.7;
   //int maxpeIterations = 10000;
   //Beschleunigung g
   double g = 9.81 * lbmUnitConverter->getFactorAccWToLb();
   //Vector3D globalLinearAcc(0.0, -g, 0.0);
   Vector3D globalLinearAcc(0.0, 0.0, -g);

   std::shared_ptr<PePhysicsEngineMaterialAdapter> planeMaterial = std::make_shared<PePhysicsEngineMaterialAdapter>("granular", 1.0, 0, 0.1 / 2, 0.1 / 2, 0.5, 1, 1, 0, 0);

   const int gridNx = val<1>(grid->getBlockNX()) * grid->getNX1();
   const int gridNy = val<2>(grid->getBlockNX()) * grid->getNX2();
   const int gridNz = val<3>(grid->getBlockNX()) * grid->getNX3();

   //UbTupleInt3 simulationDomain(gridNx, gridNy, gridNz);
   //std::array<double, 6> simulationDomain = {1, 1, 1, 30, 30, 30};
   std::array<double, 6> simulationDomain ={ g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3 };
   UbTupleInt3 numberOfBlocks(grid->getNX1(), grid->getNX2(), grid->getNX3());
   //UbTupleInt3 numberOfBlocks((simulationDomain[3]-simulationDomain[0])/val<1>(grid->getBlockNX()), (simulationDomain[4]-simulationDomain[1])/val<2>(grid->getBlockNX()), (simulationDomain[5]-simulationDomain[2])/val<3>(grid->getBlockNX()));
   UbTupleBool3 isPeriodic(grid->isPeriodicX1(), grid->isPeriodicX2(), grid->isPeriodicX3());
   Vector3D minOffset(peMinOffset[0], peMinOffset[1], peMinOffset[2]);
   Vector3D maxOffset(peMaxOffset[0], peMaxOffset[1], peMaxOffset[2]);

   SPtr<GbObject3D> boxPE(new GbCuboid3D(simulationDomain[0]+minOffset[0], simulationDomain[1]+minOffset[1], simulationDomain[2]+minOffset[2], simulationDomain[3]+maxOffset[0], simulationDomain[4]+maxOffset[1], simulationDomain[5]+maxOffset[2]));
   GbSystem3D::writeGeoObject(boxPE.get(), pathOut + "/geo/boxPE", WbWriterVtkXmlBinary::getInstance());

   std::shared_ptr<PeParameter> peParamter = std::make_shared<PeParameter>(peRelaxtion, maxpeIterations, globalLinearAcc,
      planeMaterial, simulationDomain, numberOfBlocks, isPeriodic, minOffset, maxOffset);
   std::shared_ptr<PeLoadBalancerAdapter> loadBalancer(new PeLoadBalancerAdapter(grid, comm->getNumberOfProcesses(), comm->getProcessID()));
   std::shared_ptr<PhysicsEngineSolverAdapter> peSolver = std::make_shared<PePhysicsEngineSolverAdapter>(peParamter, loadBalancer);

   const std::shared_ptr<ForceCalculator> forceCalculator = std::make_shared<ForceCalculator>(comm);

   return std::make_shared<DemCoProcessor>(grid, peScheduler, comm, forceCalculator, peSolver);
}

//pipe flow with forcing
void thermoplast(string configname)
{
   SPtr<Communicator> comm = MPICommunicator::getInstance();
   int myid = comm->getProcessID();

   ConfigurationFile   config;
   config.load(configname);

   vector<int>     blocknx = config.getVector<int>("blocknx");
   vector<double>  boundingBox = config.getVector<double>("boundingBox");

   int             endTime = config.getValue<int>("endTime");
   double          outTime = config.getValue<double>("outTime");
   double          availMem = config.getValue<double>("availMem");
   double          uLB = config.getValue<double>("uLB");

   string          michel = config.getValue<string>("michel");
   string          plexiglas = config.getValue<string>("plexiglas");
   double          sphereTime = config.getValue<double>("sphereTime");

   double          cpStart = config.getValue<double>("cpStart");
   double          cpStep = config.getValue<double>("cpStep");
   bool            restart = config.getValue<bool>("restart");

   peMinOffset = config.getVector<double>("peMinOffset");
   peMaxOffset = config.getVector<double>("peMaxOffset");

   pathOut = config.getValue<string>("pathOut");
   pathGeo = config.getValue<string>("pathGeo");

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

   //Grid definition
   SPtr<Grid3D> grid(new Grid3D(comm));
   grid->setDeltaX(deltax);
   grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
   grid->setPeriodicX1(false);
   grid->setPeriodicX2(false);
   grid->setPeriodicX3(false);

   //boundary conditions definition 
   //////////////////////////////////////////////////////////////////////////////
   SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
   noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

   mu::Parser fct;
   fct.SetExpr("U");
   fct.DefineConst("U", uLB);
   SPtr<BCAdapter> inflowAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
   inflowAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));
   //inflowAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

   SPtr<BCAdapter> outflowAdapter(new DensityBCAdapter(rhoLB));
   outflowAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new EqDensityBCAlgorithm()));
   //outflowAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));
   //outflowAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

   //sphere BC
   mu::Parser fct2;
   fct2.SetExpr("U");
   fct2.DefineConst("U", 0.0);
   SPtr<BCAdapter> velocityBcParticleAdapter(new VelocityBCAdapter(true, false, false, fct2, 0, BCFunction::INFCONST));
   velocityBcParticleAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

   //boundary conditions visitor
   SPtr<BoundaryConditionsBlockVisitor> bcVisitor(new BoundaryConditionsBlockVisitor());
   bcVisitor->addBC(noSlipBCAdapter);
   bcVisitor->addBC(inflowAdapter);
   bcVisitor->addBC(outflowAdapter);
   bcVisitor->addBC(velocityBcParticleAdapter);
   //////////////////////////////////////////////////////////////////////////////////

   //LBM kernel definition
   SPtr<LBMKernel> kernel;
   kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel());
   SPtr<BCProcessor> bcProc(new BCProcessor());
   kernel->setBCProcessor(bcProc);

   //blocks generating
   SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
   if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
   GenBlocksGridVisitor genBlocks(gridCube);
   grid->accept(genBlocks);


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

   GbCuboid3DPtr geoInflow1(new GbCuboid3D(g_minX1-blockLength, g_maxX2-120.0, g_minX3+190.0, g_minX1+1, g_maxX2+20.0, g_minX3+320.0));
   if (myid == 0) GbSystem3D::writeGeoObject(geoInflow1.get(), pathOut + "/geo/geoInflow1", WbWriterVtkXmlASCII::getInstance());

   if (!restart)
   {
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

      //Duese
      if (myid==0) UBLOG(logINFO, "Read Duese:start");
      SPtr<GbTriFaceMesh3D> s1Geo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/Duese/s1.stl", "plexiglasGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
      SPtr<GbTriFaceMesh3D> b1Geo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/Duese/b1.stl", "plexiglasGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
      SPtr<GbTriFaceMesh3D> p1Geo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/Duese/p1.stl", "plexiglasGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
      SPtr<GbTriFaceMesh3D> p2Geo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/Duese/p2.stl", "plexiglasGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
      if (myid==0) UBLOG(logINFO, "Read Duese:end");



      //inflow
      GbCuboid3DPtr geoInflow2(new GbCuboid3D(g_minX1-blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_minX1, g_maxX2 + blockLength, g_maxX3 + blockLength));
      if (myid == 0) GbSystem3D::writeGeoObject(geoInflow2.get(), pathOut + "/geo/geoInflow2", WbWriterVtkXmlASCII::getInstance());

      //GbCuboid3DPtr geoInflow3(new GbCuboid3D(g_minX1-blockLength, g_minX2-radius, g_maxX3-4.0*radius-1, g_minX1+1, g_maxX2+radius, g_maxX3-4.0*radius));
      //if (myid == 0) GbSystem3D::writeGeoObject(geoInflow3.get(), pathOut + "/geo/geoInflow3", WbWriterVtkXmlASCII::getInstance());

      //GbCuboid3DPtr geoInflow4(new GbCuboid3D(g_minX1-blockLength, g_minX2+4.0*radius, g_maxX3-4.0*radius-1.0, g_minX1+1, g_minX2+4.0*radius+1.0, g_maxX3+radius));
      //if (myid == 0) GbSystem3D::writeGeoObject(geoInflow4.get(), pathOut + "/geo/geoInflow4", WbWriterVtkXmlASCII::getInstance());

      //outflow
      GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
      if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

      //set boundary conditions for blocks and create process decomposition for MPI
      SPtr<D3Q27Interactor> boxInt(new D3Q27Interactor(box, grid, noSlipBCAdapter, Interactor3D::INVERSESOLID));

      //inflow
      SPtr<D3Q27Interactor> inflowInt1 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow1, grid, inflowAdapter, Interactor3D::SOLID));
      SPtr<D3Q27Interactor> inflowInt2 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow2, grid, outflowAdapter, Interactor3D::SOLID));

      //SPtr<D3Q27Interactor> inflowInt3 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow3, grid, noSlipBCAdapter, Interactor3D::SOLID));
      //SPtr<D3Q27Interactor> inflowInt4 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow4, grid, noSlipBCAdapter, Interactor3D::SOLID));

      //outflow
      SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowAdapter, Interactor3D::SOLID));

      //michel
      SPtr<Interactor3D> michelInt = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(michelGeo, grid, noSlipBCAdapter, Interactor3D::SOLID));

      //plexiglas
      SPtr<Interactor3D> plexiglasInt = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(plexiglasGeo, grid, noSlipBCAdapter, Interactor3D::SOLID));

      //Duese
      SPtr<Interactor3D> s1Int = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(s1Geo, grid, noSlipBCAdapter, Interactor3D::SOLID));
      SPtr<Interactor3D> b1Int = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(b1Geo, grid, noSlipBCAdapter, Interactor3D::SOLID));
      SPtr<Interactor3D> p1Int = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(p1Geo, grid, noSlipBCAdapter, Interactor3D::SOLID));
      SPtr<Interactor3D> p2Int = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(p2Geo, grid, noSlipBCAdapter, Interactor3D::SOLID));

      //////////////////////////////////////////////////////////////////////////
      //SPtr<Grid3DVisitor> peVisitor(new PePartitioningGridVisitor(comm, demCoProcessor));
      SPtr<Grid3DVisitor> peVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
      InteractorsHelper intHelper(grid, peVisitor);
      intHelper.addInteractor(boxInt);
      intHelper.addInteractor(michelInt);
      intHelper.addInteractor(plexiglasInt);
      intHelper.addInteractor(s1Int);
      intHelper.addInteractor(b1Int);
      intHelper.addInteractor(p1Int);
      intHelper.addInteractor(p2Int);
      intHelper.addInteractor(inflowInt1);
      intHelper.addInteractor(outflowInt);
      intHelper.addInteractor(inflowInt2);
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

      //create LBM kernel
      SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
      grid->accept(kernelVisitor);

      intHelper.setBC();


      //initialization of distributions
      InitDistributionsBlockVisitor initVisitor;
      //initVisitor.setVx1(uLB);
      grid->accept(initVisitor);

      //write data for visualization of boundary conditions
      {
         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
         ppgeo.process(0);

         WriteMacroscopicQuantitiesCoProcessor ppInit(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm);
         ppInit.process(0);
      }

      if (myid == 0) UBLOG(logINFO, "Preprocess - end");
   }

   //PE initialization
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
   demCoProcessor->setBlockVisitor(bcVisitor);

   ////////////////////////////////////////////////////////////////////////////
   ////generating spheres 
   //UBLOG(logINFO, "generating spheres - start, rank="<<myid);
   SPtr<UbScheduler> sphereScheduler(new UbScheduler(sphereTime));
   SPtr<CreateDemObjectsCoProcessor> createSphereCoProcessor(new CreateDemObjectsCoProcessor(grid, sphereScheduler, comm, demCoProcessor, sphereMaterial));
   //UBLOG(logINFO, "generating spheres - stop, rank="<<myid);

   //restart
   //UBLOG(logINFO, "restart definition - start, rank="<<myid);
   SPtr<UbScheduler> restartSch(new UbScheduler(cpStep, cpStart));
   SPtr<MPIIORestartCoProcessor> restartCoProcessor(new MPIIORestartCoProcessor(grid, restartSch, pathOut, comm));
   restartCoProcessor->setLBMKernel(kernel);
   restartCoProcessor->setBCProcessor(bcProc);
   SPtr<RestartDemObjectsCoProcessor> restartDemObjectsCoProcessor(new RestartDemObjectsCoProcessor(grid, restartSch, pathOut, demCoProcessor, createSphereCoProcessor, radius, comm));
   //UBLOG(logINFO, "restart definition - stop, rank="<<myid);

   if (restart)
   {
      int startStep = restartCoProcessor->readCpTimeStep();
      restartCoProcessor->restart(startStep);
      restartDemObjectsCoProcessor->restart(startStep);
   }

   //set connectors
   //UBLOG(logINFO, "set connectors - start, rank="<<myid);
   InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
   SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
   grid->accept(setConnsVisitor);
   //UBLOG(logINFO, "set connectors - stop, rank="<<myid);

   //BC visitor
   //UBLOG(logINFO, "BC visitor - start, rank="<<myid);
   grid->accept(*bcVisitor.get());
   //UBLOG(logINFO, "BC visitor - stop, rank="<<myid);

   //sphere prototypes
   //UBLOG(logINFO, "sphere prototypes - start, rank="<<myid);
   double d = 2.0*radius;
   Vector3D origin(g_minX1+peMinOffset[0]+radius, geoInflow1->getX2Minimum()+4.0*d, geoInflow1->getX3Minimum()+2.0*d);
   for (int x3 = 0; x3 < 2; x3++)
      for (int x2 = 0; x2 < 2; x2++)
         for (int x1 = 0; x1 < 1; x1++)
         {
            //SPtr<GbObject3D> sphere(new GbSphere3D(origin[0]+x1*d, origin[1]+x2*2.0*d, origin[2]+x3*2.0*d, radius));
            SPtr<GbObject3D> sphere(new GbSphere3D(origin[0]+2.0*d, origin[1]+x2*1.5*d, origin[2]+x3*1.5*d, radius));
            if (myid == 0) GbSystem3D::writeGeoObject(sphere.get(), pathOut + "/geo/sphere"+UbSystem::toString(x1)+UbSystem::toString(x2)+UbSystem::toString(x3), WbWriterVtkXmlASCII::getInstance());
            createSphereCoProcessor->addGeoObject(sphere, Vector3D(uLB, 0.0, 0.0));
         }
   //UBLOG(logINFO, "sphere prototypes - stop, rank="<<myid);

   //Vector3D origin(106+radius, 1372+radius, 12+radius);
   //for (int x3 = 0; x3 < 28; x3++)
   //   for (int x2 = 0; x2 < 12; x2++)
   //      for (int x1 = 0; x1 < 7; x1++)
   //      {
   //         //SPtr<GbObject3D> sphere(new GbSphere3D(origin[0]+x1*d, origin[1]+x2*2.0*d, origin[2]+x3*2.0*d, radius));
   //         SPtr<GbObject3D> sphere(new GbSphere3D(origin[0]+x1*1.1*d, origin[1]+x2*1.1*d, origin[2]+x3*1.1*d, radius));
   //         //if (myid == 0) GbSystem3D::writeGeoObject(sphere.get(), pathOut + "/geo/sphere"+UbSystem::toString(x1)+UbSystem::toString(x2)+UbSystem::toString(x3), WbWriterVtkXmlASCII::getInstance());
   //         createSphereCoProcessor->addGeoObject(sphere, Vector3D(uLB, 0.0, 0.0));
   //      }

   createSphereCoProcessor->process(0);

   //write data for visualization of macroscopic quantities
   SPtr<UbScheduler> visSch(new UbScheduler(outTime));
   SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, pathOut,
      WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
   writeMQCoProcessor->process(0);

   SPtr<WriteBoundaryConditionsCoProcessor> writeBCCoProcessor(new WriteBoundaryConditionsCoProcessor(grid, visSch, pathOut,
      WbWriterVtkXmlBinary::getInstance(), comm));
   writeBCCoProcessor->process(0);

   SPtr<WriteDemObjectsCoProcessor> writeDemObjectsCoProcessor(new WriteDemObjectsCoProcessor(grid, visSch, pathOut, WbWriterVtkXmlBinary::getInstance(), demCoProcessor, comm));
   writeDemObjectsCoProcessor->process(0);

   ////performance control
   SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
   SPtr<NUPSCounterCoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

   //start simulation 
   omp_set_num_threads(numOfThreads);
   SPtr<UbScheduler> stepGhostLayer(peScheduler);
   SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
   calculator->addCoProcessor(npr);

   calculator->addCoProcessor(createSphereCoProcessor);
   calculator->addCoProcessor(demCoProcessor);
   calculator->addCoProcessor(writeBCCoProcessor);
   calculator->addCoProcessor(writeDemObjectsCoProcessor);
   calculator->addCoProcessor(writeMQCoProcessor);
   //calculator->addCoProcessor(restartDemObjectsCoProcessor);
   //calculator->addCoProcessor(restartCoProcessor);

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
   }
   catch (std::exception& e)
   {
      UBLOG(logERROR, e.what());
   }
   catch (std::string& s)
   {
      UBLOG(logERROR, s);
   }
   catch (...)
   {
      UBLOG(logERROR, "unknown exception");
   }
}

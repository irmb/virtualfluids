#include <VirtualFluids.h>
#include <set>
#include <map>
using namespace std;


////////////////////////////////////////////////////////////////////////
void run(string configname)
{
   try
   {

      //Sleep(30000);

      string machine = QUOTEME(CAB_MACHINE);

      SPtr<Communicator> comm = MPICommunicator::getInstance();

      int myid = comm->getProcessID();
      int mybundle = comm->getBundleID();
      int root = comm->getRoot();

      vf::basics::ConfigurationFile   config;
      config.load(configname);

      string pathname = config.getValue<string>("path");
      double availMem = config.getValue<double>("memory");
      double outstep = config.getValue<double>("outstep");
      double endstep = config.getValue<double>("endstep");
      int numOfThreads = config.getValue<int>("threads");
      int refineLevel = config.getValue<int>("level");
      vector<double> dim = config.getVector<double>("dim");
      vector<int> blockNx = config.getVector<int>("blockNx");
      double radius = config.getValue<double>("radius");

      //LBMReal radius = 4;
      LBMReal uLB = 0.1;
      LBMReal Re = 1;
      LBMReal rhoLB = 0.0;
      //LBMReal nuLB = (uLB*2.0*radius)/Re;
      //LBMReal nuLB = (uLB*L2)/Re;
      LBMReal nuLB = 0.168666666667/100;

      double dp_LB = 1e-6;
      double rhoLBinflow = dp_LB*3.0;

      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);
      SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapter(new DensityBCAdapter(rhoLB));
      denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(velBCAdapter);
      bcVisitor.addBC(denBCAdapter);

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      double dx = 1;

      const int blocknx1 = blockNx[0];
      const int blocknx2 = blockNx[1];
      const int blocknx3 = blockNx[2];

      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);

      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> restartSch(new UbScheduler(100000, 100000, 100000));
      RestartCoProcessor rp(grid, restartSch, comm, pathname, RestartCoProcessor::BINARY);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep()==0)
      {

         const int baseLevel = 0;

         //bounding box
         double d_minX1 = 0.0;
         double d_minX2 = 0.0;
         double d_minX3 = 0.0;

         double d_maxX1 = dim[0];
         double d_maxX2 = dim[1];
         double d_maxX3 = dim[2];

         double blockLength = blocknx1*dx;

         if (myid==0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "uLB = "<<uLB);
            UBLOG(logINFO, "rhoLB = "<<rhoLB);
            UBLOG(logINFO, "nueLB = "<<nuLB);
            UBLOG(logINFO, "Re = "<<Re);
            UBLOG(logINFO, "dx = "<<dx);
            UBLOG(logINFO, "number of levels = "<<refineLevel+1);
            UBLOG(logINFO, "numOfThreads = "<<numOfThreads);
            UBLOG(logINFO, "Preprozess - start");
         }

         SPtr<GbObject3D> gridCube(new GbCuboid3D(d_minX1, d_minX2, d_minX3, d_maxX1, d_maxX2, d_maxX3));
         if (myid==0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         //SPtr<CoordinateTransformation3D> trafo = grid->getCoordinateTransformator();
         //trafo->setRotationX2Angle(4);

         //sphere
         //SPtr<GbObject3D> sphereRef(new GbSphere3D(L1/4.0, L2*0.5, L3*0.5, radius+1.0));
         //GbSystem3D::writeGeoObject(sphereRef.get(),pathname + "/geo/sphereRef", WbWriterVtkXmlBinary::getInstance());


         //sphere
         SPtr<GbObject3D> sphere(new GbSphere3D(d_maxX1*0.5, d_maxX2*0.5, d_maxX3*0.5, radius));
         //SPtr<GbObject3D> sphere(new GbSphere3D(L1/2.0-4.0, L2*0.5+4.0, L3*0.5+4.0, radius));
         //SPtr<GbObject3D> sphere(new GbCuboid3D(L1/4.0-radius, L2/2.0-radius, L3/2.0-radius, L1/4.0+radius, L2/2.0+radius, L3/2.0+radius));
         GbSystem3D::writeGeoObject(sphere.get(), pathname+"/geo/sphere", WbWriterVtkXmlBinary::getInstance());

         double off = 0.0;
         SPtr<GbObject3D> refCube(new GbCuboid3D(sphere->getX1Minimum()-off, sphere->getX2Minimum()-off, sphere->getX3Minimum(),
            sphere->getX1Maximum()+off, sphere->getX2Maximum()+off, sphere->getX3Maximum()));
         if (myid==0) GbSystem3D::writeGeoObject(refCube.get(), pathname+"/geo/refCube", WbWriterVtkXmlBinary::getInstance());

         if (refineLevel>0)
         {
            if (myid==0) UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
            //refineHelper.addGbObject(sphere, refineLevel);
            refineHelper.addGbObject(refCube, refineLevel);
            refineHelper.refine();
            if (myid==0) UBLOG(logINFO, "Refinement - end");
         }

         //walls
         GbCuboid3DPtr addWallYmin(new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_minX2, d_maxX3+4.0*blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin(new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_minX3));
         if (myid==0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax(new GbCuboid3D(d_minX1-4.0*blockLength, d_maxX2, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_maxX3, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_minX1, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(d_maxX1, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         WriteBlocksSPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));



         //sphere
         SPtr<D3Q27Interactor> sphereInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(sphere, grid, noSlipBCAdapter, Interactor3D::SOLID));

         //walls
         SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));

         mu::Parser fct;
         fct.SetExpr("U");
         fct.DefineConst("U", uLB);

         //inflow
         SPtr<D3Q27Interactor> inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //D3Q27BoundaryConditionAdapterPtr denBCAdapterInflow(new D3Q27DensityBCAdapter(rhoLBinflow));
         //denBCAdapterInflow->setSecondaryBcOption(0);
         //SPtr<D3Q27Interactor> inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, denBCAdapterInflow, Interactor3D::SOLID));

         //outflow
         SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));

         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         InteractorsHelper intHelper(grid, metisVisitor);
         //intHelper.addInteractor(sphereInt);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);


         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //Block3DSPtr<ConnectorFactory> factory(new Block3DConnectorFactory());
         //ConnectorBlockVisitor setConnsVisitor(comm, nuLB, iProcessor, factory);
         //grid->accept(setConnsVisitor);

         ppblocks->process(0);
         ppblocks.reset();

         unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
         int ghostLayer = 3;
         unsigned long long numberOfNodesPerBlock = (unsigned long long)(blockNx[0])* (unsigned long long)(blockNx[1])* (unsigned long long)(blockNx[2]);
         unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
         unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blockNx[0]+ghostLayer) * (blockNx[1]+ghostLayer) * (blockNx[2]+ghostLayer);
         double needMemAll = double(numberOfNodesPerBlockWithGhostLayer*(27*sizeof(double)+sizeof(int)+sizeof(float)*4));
         double needMem = needMemAll/double(comm->getNumberOfProcesses());

         if (myid==0)
         {
            UBLOG(logINFO, "Number of blocks = "<<numberOfBlocks);
            UBLOG(logINFO, "Number of nodes  = "<<numberOfNodes);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for (int level = minInitLevel; level<=maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO, "Number of blocks for level "<<level<<" = "<<nobl);
               UBLOG(logINFO, "Number of nodes for level "<<level<<" = "<<nobl*numberOfNodesPerBlock);
            }
            UBLOG(logINFO, "Necessary memory  = "<<needMemAll<<" bytes");
            UBLOG(logINFO, "Necessary memory per process = "<<needMem<<" bytes");
            UBLOG(logINFO, "Available memory per process = "<<availMem<<" bytes");
         }

         SPtr<LBMKernel> kernel(new IncompressibleCumulantLBMKernel(blocknx1, blocknx2, blocknx3, IncompressibleCumulantLBMKernel::NORMAL));

         SPtr<BCProcessor> bcProcessor(new BCProcessor());


         kernel->setBCProcessor(bcProcessor);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel>0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         intHelper.setBC();

         grid->accept(bcVisitor);

         mu::Parser fctRoh;
         fctRoh.SetExpr("(x1max-x1)/l*dp*3.0");
         fctRoh.DefineConst("dp", dp_LB);
         fctRoh.DefineConst("x1max", d_maxX1);
         fctRoh.DefineConst("l", d_maxX1-d_minX1);

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         initVisitor.setVx1(fct);
         //initVisitor.setRho(fctRoh);
         grid->accept(initVisitor);

         //Postrozess
         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         WriteBoundaryConditionsSPtr<CoProcessor> ppgeo(
            new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         ppgeo->process(0);
         ppgeo.reset();;

         if (myid==0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         UBLOG(logINFO, "SetConnectors - start, id="<<myid);

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         //D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         SPtr<ConnectorFactory> cFactory(new Block3DConnectorFactory());
         ConnectorBlockVisitor setConnsVisitor(comm, nuLB, iProcessor, cFactory);
         grid->accept(setConnsVisitor);

         UBLOG(logINFO, "SetConnectors - end, id="<<myid);
      }

      SPtr<UbScheduler> stepSch(new UbScheduler(outstep));
      //stepSch->addSchedule(10000, 0, 1000000);
      WriteMacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      const SPtr<ConcreteCalculatorFactory> calculatorFactory = std::make_shared<ConcreteCalculatorFactory>(stepSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endstep, calculatorFactory, CalculatorType::HYBRID));

      if (myid==0)
         UBLOG(logINFO, "Simulation-start");

      calculation->calculate();

      if (myid==0)
         UBLOG(logINFO, "Simulation-end");

   }
   catch (std::exception& e)
   {
      cerr<<e.what()<<endl<<flush;
   }
   catch (std::string& s)
   {
      cerr<<s<<endl;
   }
   catch (...)
   {
      cerr<<"unknown exception"<<endl;
   }

}

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   if (argv!=NULL)
   {
      if (argv[1]!=NULL)
      {
         run(string(argv[1]));
      }
      else
      {
         cout<<"Configuration file is missing!"<<endl;
      }
   }
}


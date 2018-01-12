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
      double          uLB = config.getDouble("uLB");
      double          endTime = config.getDouble("endTime");
      double          outTime = config.getDouble("outTime");
      double          availMem = config.getDouble("availMem");
      int             refineLevel = config.getInt("refineLevel");
      double          Re = config.getDouble("Re");
      double          dx = config.getDouble("dx");
      vector<double>  length = config.getVector<double>("length");
      bool            logToFile = config.getBool("logToFile");
      double          restartStep = config.getDouble("restartStep");
      double          cpStart = config.getValue<double>("cpStart");
      double          cpStep = config.getValue<double>("cpStep");
      bool            newStart = config.getValue<bool>("newStart");

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

      //Sleep(30000);

      LBMReal dLB = length[1] / dx;
      LBMReal rhoLB = 0.0;
      LBMReal nuLB = (uLB*dLB) / Re;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      //BC Adapter
      //////////////////////////////////////////////////////////////////////////////
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      //noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new ThinWallNoSlipBCAlgorithm()));
      //SPtr<BCAdapter> denBCAdapter(new DensityBCAdapter(rhoLB));
      //denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new EqDensityBCAlgorithm()));

      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapter(new DensityBCAdapter(rhoLB));
      denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

      //mu::Parser fct;
      //fct.SetExpr("U");
      //fct.DefineConst("U", uLB);
      //SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      //velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingVelocityBCAlgorithm()));


      //////////////////////////////////////////////////////////////////////////////////
      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(denBCAdapter);
      //bcVisitor.addBC(velBCAdapter);
      
      SPtr<Grid3D> grid(new Grid3D(comm));
      
      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
      //RestartCoProcessor rp(grid, rSch, comm, pathname, RestartCoProcessor::TXT);
      MPIIORestart1CoProcessor rcp(grid, rSch, pathname, comm);
      //////////////////////////////////////////////////////////////////////////

      if (newStart)
      {

         //bounding box
         double g_minX1 = 0.0;
         double g_minX2 = -length[1] / 2.0;
         double g_minX3 = -length[2] / 2.0;

         double g_maxX1 = length[0];
         double g_maxX2 = length[1] / 2.0;
         double g_maxX3 = length[2] / 2.0;

         //geometry
         //x
         //SPtr<GbObject3D> cylinder(new GbCylinder3D(g_minX1 - 2.0*dx, 0.0, 0.0, g_maxX1 + 2.0*dx, 0.0, 0.0, dLB / 2.0));
         //y
         SPtr<GbObject3D> cylinder(new GbCylinder3D(g_minX1 - 2.0*dx, 0.0, 0.0, g_maxX1 + 2.0*dx, 0.0, 0.0, dLB / 2.0));
         GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


         double blockLength = blocknx[0] * dx;



         if (myid == 0)
         {
            UBLOG(logINFO, "uLb = " << uLB);
            UBLOG(logINFO, "rho = " << rhoLB);
            UBLOG(logINFO, "nuLb = " << nuLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << dx);
            UBLOG(logINFO, "Preprocess - start");
         }

         grid->setDeltaX(dx);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

         grid->setPeriodicX1(false);
         grid->setPeriodicX2(true);
         grid->setPeriodicX3(false);

         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         ////inflow
         //GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1 - 2.0*dx, g_minX2 - dx, g_minX3 - dx, g_minX1, g_maxX2, g_maxX3 + dx));
         //if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         ////outflow
         //GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2, g_minX3 - dx, g_maxX1 + 2.0*dx, g_maxX2 + dx, g_maxX3 + dx));
         //if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         WriteBlocksSPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         ppblocks->process(0);

         SPtr<D3Q27Interactor> cylinderInt(new D3Q27Interactor(cylinder, grid, noSlipBCAdapter, Interactor3D::INVERSESOLID));

         double r = dynamicPointerCast<GbCylinder3D>(cylinder)->getRadius();
         double cx1 = g_minX1;
         double cx2 = cylinder->getX2Centroid();
         double cx3 = cylinder->getX3Centroid();

         mu::Parser fct;
         fct.SetExpr("vx1*(1-((x2-y0)^2+(x3-z0)^2)/(R^2))");
         fct.DefineConst("x2Vmax", 0.0); //x2-Pos fuer vmax
         fct.DefineConst("x3Vmax", 0.0); //x3-Pos fuer vmax
         fct.DefineConst("R", r);
         fct.DefineConst("vx1", uLB);
         fct.DefineConst("x0", cx1);
         fct.DefineConst("y0", cx2);
         fct.DefineConst("z0", cx3);
         fct.DefineConst("nue", nuLB);
         SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
         //velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));
         velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));


         SPtr<D3Q27Interactor> inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //outflow
         SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));


         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(cylinderInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();


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

         SPtr<LBMKernel> kernel;

         kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], IncompressibleCumulantLBMKernel::NORMAL));
         //kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], CompressibleCumulantLBMKernel::NORMAL));

         //
         SPtr<BCProcessor> bcProc(new BCProcessor());
         //SPtr<BCProcessor> bcProc(new ThinWallBCProcessor());

         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         intHelper.setBC();

         bcVisitor.addBC(velBCAdapter);
         grid->accept(bcVisitor);

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         //initVisitor.setVx1(fct);
         //initVisitor.setVx1(uLB);
         grid->accept(initVisitor);

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         //InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         //SPtr<ConnectorFactory> factory(new Block3DConnectorFactory());
         //ConnectorBlockVisitor setConnsVisitor(comm, nuLB, iProcessor, factory);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //boundary conditions grid
         {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            WriteBoundaryConditionsSPtr<CoProcessor> ppgeo(
               new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
            ppgeo->process(0);
            ppgeo.reset();
         }

         if (myid == 0) UBLOG(logINFO, "Preprocess - end");
      }
      else
      {
         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "uLb = " << uLB);
            UBLOG(logINFO, "rho = " << rhoLB);
            UBLOG(logINFO, "nuLb = " << nuLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << dx);
            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
         }

         rcp.restart((int)restartStep);
         grid->setTimeStep(restartStep);

         SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter());
         //velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));
         velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));
         bcVisitor.addBC(velBCAdapter);
         grid->accept(bcVisitor);

         //set connectors
         //InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      WriteMacroscopicQuantitiesCoProcessor pp(grid, visSch, pathname, WbWriterVtkXmlASCII::getInstance(), conv, comm);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch,CalculationManager::MPI));


      const SPtr<ConcreteCalculatorFactory> calculatorFactory = std::make_shared<ConcreteCalculatorFactory>(visSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, calculatorFactory, CalculatorType::PREPOSTBC));
      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
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


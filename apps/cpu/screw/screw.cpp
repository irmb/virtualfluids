#include <VirtualFluids.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   try
   {
      // Verify input arguments
      if (argc != 2)
      {
         std::cout << "Usage: " << argv[0] << " <config file> " << std::endl;
         return EXIT_FAILURE;
      }

      vf::basics::ConfigurationFile   config;
      config.load(argv[1]);

      string          pathname     = config.getString("pathname");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx      = config.getVector<int>("blocknx");
      double          uLB          = config.getValue<double>("uLB");
      double          Re           = config.getValue<double>("Re");
      double          dx           = config.getValue<double>("dx");
      double          endTime      = config.getValue<double>("endTime");
      double          outTime      = config.getValue<double>("outTime");
      double          availMem     = config.getValue<double>("availMem");
      int             restartStep  = config.getValue<int>("restartStep");


      SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      double r = 150e-3;
      double h = 45e-3;

      //bounding box
      double g_minX1 = 0;
      double g_minX2 = 0;
      double g_minX3 = 0;

      double g_maxX1 = 2.0*sqrt(2.0*r*h-h*h);
      double g_maxX2 = 45e-3;
      double g_maxX3 = 10e-3;

      

      double blockLength = (double)blocknx[0]*dx;

      double nuLB = (uLB*(h/dx))/Re;
      double rhoLB = 0.0;

      //bc
      mu::Parser fctVx;
      fctVx.SetExpr("omega*(r-x2)");
      fctVx.DefineConst("omega", uLB);
      fctVx.DefineConst("r", r);

      mu::Parser fctVy;
      fctVy.SetExpr("omega*(x1-k)");
      fctVy.DefineConst("omega", uLB);
      fctVy.DefineConst("k", g_maxX1*0.5);

      mu::Parser fctVz;
      fctVz.SetExpr("0.0");

      SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, true, true, fctVx,fctVy,fctVz, 0, BCFunction::INFCONST));
      velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));

      SPtr<BCAdapter> slipBCAdapter(new SlipBCAdapter());
      slipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SlipBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(slipBCAdapter);
      bcVisitor.addBC(velBCAdapter);


      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setPeriodicX1(false);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(false);
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> rSch(new UbScheduler(restartStep));
      RestartCoProcessor rp(grid, rSch, comm, pathname, RestartCoProcessor::TXT);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep() == 0)
      {
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "h = " << h);
            UBLOG(logINFO, "rho = " << rhoLB);
            UBLOG(logINFO, "nue = " << nuLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << dx);

            //UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
            UBLOG(logINFO, "Preprozess - start");
         }

         //BC
         SPtr<GbObject3D> cylinder(new GbCylinder3D(g_maxX1*0.5, r, g_minX3, g_maxX1*0.5, r, g_maxX3, r));
         GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_maxX2, g_minX3-blockLength, g_maxX1+blockLength,  g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         WriteBlocksSPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         //interactors
         SPtr<D3Q27Interactor> cylinderInt(new D3Q27Interactor(cylinder, grid, velBCAdapter, Interactor3D::INVERSESOLID));
         SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallZmax, grid, slipBCAdapter, Interactor3D::SOLID));


         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(cylinderInt);

         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolator());
         //InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolator());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         ppblocks->process(0);
         ppblocks.reset();

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nod = nob * (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nodg = nob * (blocknx[0] + gl) * (blocknx[1] + gl) * (blocknx[1] + gl);
         double needMemAll = double(nodg*(27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
         double needMem = needMemAll / double(comm->getNumberOfProcesses());

         if (myid == 0)
         {
            UBLOG(logINFO, "Number of blocks = " << nob);
            UBLOG(logINFO, "Number of nodes  = " << nod);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for (int level = minInitLevel; level <= maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO, "Number of blocks for level " << level << " = " << nobl);
               UBLOG(logINFO, "Number of nodes for level " << level << " = " << nobl*nodb);
            }
            UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
            UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
            UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
         }

         int kernelType = 2;
         SPtr<LBMKernel> kernel;
         kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], IncompressibleCumulantLBMKernel::NORMAL));
         //kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], CompressibleCumulantLBMKernel::NORMAL));
         
         SPtr<BCProcessor> bcProc(new BCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         //BC
         intHelper.setBC();

         grid->accept(bcVisitor);

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         initVisitor.setVx1(fctVx);
         initVisitor.setVx2(fctVy);
         initVisitor.setVx3(fctVz);
         grid->accept(initVisitor);

         //Postrozess
         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         WriteBoundaryConditionsSPtr<CoProcessor> ppgeo(
            new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         ppgeo->process(0);
         ppgeo.reset();

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         grid->accept(bcVisitor);

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolator());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      SPtr<UbScheduler> stepSch(new UbScheduler(outTime));

      WriteMacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm);

      const SPtr<ConcreteCalculatorFactory> calculatorFactory = std::make_shared<ConcreteCalculatorFactory>(stepSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, calculatorFactory, CalculatorType::HYBRID));
      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
      if (myid == 0) UBLOG(logINFO, "Simulation-end");

      return EXIT_SUCCESS;
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

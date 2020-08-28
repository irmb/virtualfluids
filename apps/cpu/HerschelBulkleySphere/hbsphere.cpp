#include <iostream>
#include <string>

#include <VirtualFluids.h>

using namespace std;


void bflow(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          outputPath = config.getValue<string>("outputPath");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
      vector<double>  boundingBox = config.getVector<double>("boundingBox");
      //double          nuLB = config.getValue<double>("nuLB");
      double          endTime = config.getValue<double>("endTime");
      double          outTime = config.getValue<double>("outTime");
      double          availMem = config.getValue<double>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      bool            logToFile = config.getValue<bool>("logToFile");
      double          restartStep = config.getValue<double>("restartStep");
      double          deltax = config.getValue<double>("deltax");
      double          radius = config.getValue<double>("radius");
      double          cpStep = config.getValue<double>("cpStep");
      double          cpStart = config.getValue<double>("cpStart");
      bool            newStart = config.getValue<bool>("newStart");
      double          velocity = config.getValue<double>("velocity");
      double          n = config.getValue<double>("n");
      double          Re = config.getValue<double>("Re");
      double          Bn = config.getValue<double>("Bn");

      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (logToFile)
      {
#if defined(__unix__)
         if (myid == 0)
         {
            const char* str = outputPath.c_str();
            mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         }
#endif 

         if (myid == 0)
         {
            stringstream logFilename;
            logFilename << outputPath + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }

      LBMReal rhoLB = 0.0;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      //bounding box

      //double g_minX1 = 0;
      //double g_minX2 = 0;
      //double g_minX3 = 0;

      //double g_maxX1 = boundingBox[0];
      //double g_maxX2 = boundingBox[1];
      //double g_maxX3 = boundingBox[2];

      double g_minX1 = -boundingBox[0]/2.0;
      double g_minX2 = -boundingBox[1] / 2.0;
      double g_minX3 = -boundingBox[2]/2.0;

      double g_maxX1 = boundingBox[0]/2.0;
      double g_maxX2 = boundingBox[1]/2.0;
      double g_maxX3 = boundingBox[2]/2.0;

      double blockLength = 3.0 * deltax;

      double d = 2.0 * radius;
      double U = velocity;
      double Gamma = U / d;

      double k = (U * d) / (Re * std::pow(Gamma, n - 1));

      double tau0 = Bn * k * std::pow(Gamma, n);

      SPtr<Thixotropy> thix = Thixotropy::getInstance();
      thix->setPowerIndex(n);
      thix->setViscosityParameter(k);
      thix->setYieldStress(tau0);

      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new HerschelBulkleyModelNoSlipBCAlgorithm()));
      //noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new BinghamModelNoSlipBCAlgorithm()));

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", velocity);
      SPtr<BCAdapter> velocityBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      velocityBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SimpleVelocityBCAlgorithm()));

      SPtr<BCAdapter> densityBCAdapter(new DensityBCAdapter());
      densityBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(velocityBCAdapter);
      bcVisitor.addBC(densityBCAdapter);

      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new BCProcessor());

      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new HerschelBulkleyModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new BinghamModelLBMKernel());
      kernel->setBCProcessor(bcProc);

      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setPeriodicX1(false);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(false);
      grid->setDeltaX(deltax);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), outputPath + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

      //sphere
      SPtr<GbObject3D> sphere(new GbSphere3D(0, 0, 0, radius));
      GbSystem3D::writeGeoObject(sphere.get(), outputPath + "/geo/sphere", WbWriterVtkXmlBinary::getInstance());
      SPtr<D3Q27Interactor> sphereInt(new D3Q27Interactor(sphere, grid, noSlipBCAdapter, Interactor3D::SOLID));

      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIOMigrationBECoProcessor> restartCoProcessor(new MPIIOMigrationBECoProcessor(grid, mSch, outputPath, comm));
      restartCoProcessor->setLBMKernel(kernel);
      restartCoProcessor->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////

      if (myid == 0)
      {
         UBLOG(logINFO, "Parameters:");
         //UBLOG(logINFO, "forcing = " << forcing);
         UBLOG(logINFO, "rho = " << rhoLB);
         UBLOG(logINFO, "U = " << U);
         UBLOG(logINFO, "Re = " << Re);
         UBLOG(logINFO, "Bn = " << Bn);
         UBLOG(logINFO, "k = " << k);
         UBLOG(logINFO, "n = " << n);
         UBLOG(logINFO, "tau0 = " << tau0);
         UBLOG(logINFO, "deltax = " << deltax);
         //UBLOG(logINFO, "number of levels = " << refineLevel + 1);
         UBLOG(logINFO, "number of threads = " << numOfThreads);
         UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
         UBLOG(logINFO, "Preprozess - start");
      }

      if (newStart)
      {
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         if (refineLevel > 0)
         {
            GbCuboid3DPtr refCube(new GbCuboid3D(-10, -10, -10, 0, 0, 0));
            if (myid == 0) GbSystem3D::writeGeoObject(refCube.get(), outputPath + "/geo/refCube", WbWriterVtkXmlASCII::getInstance());
            
            if (myid == 0) UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
            //refineHelper.addGbObject(sphere, refineLevel);
            refineHelper.addGbObject(refCube, refineLevel);
            refineHelper.refine();
            if (myid == 0) UBLOG(logINFO, "Refinement - end");
         }

         //walls
         GbCuboid3DPtr wallZmin(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(wallZmin.get(), outputPath + "/geo/wallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr wallZmax(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_maxX3, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(wallZmax.get(), outputPath + "/geo/wallZmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr wallYmin(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_minX2, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(wallYmin.get(), outputPath + "/geo/wallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr wallYmax(new GbCuboid3D(g_minX1 - blockLength, g_maxX2, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(wallYmax.get(), outputPath + "/geo/wallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr wallXmin(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_minX1, g_maxX2 + blockLength, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(wallXmin.get(), outputPath + "/geo/wallXmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr wallXmax(new GbCuboid3D(g_maxX1, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(wallXmax.get(), outputPath + "/geo/wallXmax", WbWriterVtkXmlASCII::getInstance());

         //wall interactors
         SPtr<D3Q27Interactor> wallZminInt(new D3Q27Interactor(wallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> wallZmaxInt(new D3Q27Interactor(wallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));

         SPtr<D3Q27Interactor> wallYminInt(new D3Q27Interactor(wallYmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> wallYmaxInt(new D3Q27Interactor(wallYmax, grid, noSlipBCAdapter, Interactor3D::SOLID));

         SPtr<D3Q27Interactor> wallXminInt(new D3Q27Interactor(wallXmin, grid, velocityBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> wallXmaxInt(new D3Q27Interactor(wallXmax, grid, densityBCAdapter, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(wallZminInt);
         intHelper.addInteractor(wallZmaxInt);
         intHelper.addInteractor(wallYminInt);
         intHelper.addInteractor(wallYmaxInt);
         intHelper.addInteractor(wallXminInt);
         intHelper.addInteractor(wallXmaxInt);
         //intHelper.addInteractor(sphereInt);
         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         SPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), outputPath, WbWriterVtkXmlBinary::getInstance(), comm));
         ppblocks->process(0);

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nod = nob * (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nodg = nob * (blocknx[0] + gl) * (blocknx[1] + gl) * (blocknx[1] + gl);
         double needMemAll = double(nodg * (27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
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
               UBLOG(logINFO, "Number of nodes for level " << level << " = " << nobl * nodb);
            }
            UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
            UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
            UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
         }

         SetKernelBlockVisitor kernelVisitor(kernel, k, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //BC
         intHelper.setBC();

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor;
         grid->accept(initVisitor);


         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         restartCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);
      }
      
      omp_set_num_threads(numOfThreads);

      //set connectors
      InterpolationProcessorPtr iProcessor(new ThixotropyInterpolationProcessor());
      SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, k, iProcessor);
      grid->accept(setConnsVisitor);

      grid->accept(bcVisitor);

      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      WriteBoundaryConditionsCoProcessor ppgeo = WriteBoundaryConditionsCoProcessor(grid, geoSch, outputPath, WbWriterVtkXmlASCII::getInstance(), comm);
      ppgeo.process(0);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      SPtr<CoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      //write data for visualization of macroscopic quantities
      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      //SPtr<UbScheduler> visSch(new UbScheduler(10,1));
      SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, outputPath, WbWriterVtkXmlASCII::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
      writeMQCoProcessor->process(0);

      double area = 4*UbMath::PI*radius*radius;
      SPtr<UbScheduler> forceSch(new UbScheduler(1000));
      SPtr<CalculateForcesCoProcessor> fp = make_shared<CalculateForcesCoProcessor>(grid, forceSch, outputPath + "/forces/forces.txt", comm, velocity, area);
      fp->addInteractor(sphereInt);

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
      calculator->addCoProcessor(npr);
      calculator->addCoProcessor(fp);
      calculator->addCoProcessor(writeMQCoProcessor);
      calculator->addCoProcessor(restartCoProcessor);

      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
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

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   if (argv != NULL)
   {
      if (argv[1] != NULL)
      {
         //pflowForcing(string(argv[1]));
         bflow(string(argv[1]));
      }
      else
      {
         cout << "Configuration file is missing!" << endl;
      }
   }

   return 0;
}

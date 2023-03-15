#include <iostream>
#include <string>

#include <VirtualFluids.h>

using namespace std;


void bflow(string configname)
{
    using namespace vf::lbm::dir;

   try
   {
      vf::basics::ConfigurationFile config;
      config.load(configname);

      string          outputPath = config.getValue<string>("outputPath");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
      vector<real>  boundingBox = config.getVector<real>("boundingBox");
      //double          nuLB = config.getValue<double>("nuLB");
      real          endTime = config.getValue<real>("endTime");
      real          outTime = config.getValue<real>("outTime");
      real          availMem = config.getValue<real>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      bool            logToFile = config.getValue<bool>("logToFile");
      real          restartStep = config.getValue<real>("restartStep");
      real          deltax = config.getValue<real>("deltax");
      real          radius = config.getValue<real>("radius");
      real          cpStep = config.getValue<real>("cpStep");
      real          cpStart = config.getValue<real>("cpStart");
      bool            newStart = config.getValue<bool>("newStart");
      real          velocity = config.getValue<real>("velocity");
      real          n = config.getValue<real>("n");
      real          Re = config.getValue<real>("Re");
      real          Bn = config.getValue<real>("Bn");
      vector<real>  sphereCenter = config.getVector<real>("sphereCenter");

      SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
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

      real rhoLB = 0.0;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      //bounding box

      real g_minX1 = 0;
      real g_minX2 = 0;
      real g_minX3 = 0;

      real g_maxX1 = boundingBox[0];
      real g_maxX2 = boundingBox[1];
      real g_maxX3 = boundingBox[2];

      //double g_minX1 = -boundingBox[0]/2.0;
      //double g_minX2 = -boundingBox[1] / 2.0;
      //double g_minX3 = -boundingBox[2]/2.0;

      //double g_maxX1 = boundingBox[0]/2.0;
      //double g_maxX2 = boundingBox[1]/2.0;
      //double g_maxX3 = boundingBox[2]/2.0;

      real blockLength = 3.0 * deltax;

      real d = 2.0 * radius;
      real U = velocity;
      real Gamma = U / d;

      real k = (U * d) / (Re * std::pow(Gamma, n - 1));
      real tau0 = Bn * k * std::pow(Gamma, n);

      //double k = 0.05; // (U * d) / (Re * std::pow(Gamma, n - 1));
      //double tau0 = 3e-6; //Bn * k * std::pow(Gamma, n);

      //double forcing = 8e-7;

      real omegaMin = 1.0e-8;

      SPtr<Rheology> thix = Rheology::getInstance();
      thix->setPowerIndex(n);
      thix->setViscosityParameter(k);
      thix->setYieldStress(tau0);
      thix->setOmegaMin(omegaMin);

      SPtr<BC> noSlipBC(new NoSlipBC());
      //noSlipBC->setBCStrategy(SPtr<BCStrategy>(new NoSlipBCStrategy()));
      noSlipBC->setBCStrategy(SPtr<BCStrategy>(new RheologyHerschelBulkleyModelNoSlipBCStrategy()));
      //noSlipBC->setBCStrategy(SPtr<BCStrategy>(new RheologyBinghamModelNoSlipBCStrategy()));

      SPtr<BC> slipBC(new SlipBC());
      slipBC->setBCStrategy(SPtr<BCStrategy>(new SimpleSlipBCStrategy()));

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", velocity);
      //double H = boundingBox[1];
      //mu::Parser fct;
      //fct.SetExpr("16*U*x2*x3*(H-x2)*(H-x3)/H^4");
      //fct.DefineConst("U", U);
      //fct.DefineConst("H", H);
      SPtr<BC> velocityBC(new VelocityBC(true, false, false, fct, 0, BCFunction::INFCONST));
      velocityBC->setBCStrategy(SPtr<BCStrategy>(new SimpleVelocityBCStrategy()));
      //velocityBC->setBCStrategy(SPtr<BCStrategy>(new VelocityWithDensityBCStrategy()));

      SPtr<BC> densityBC(new DensityBC());
      densityBC->setBCStrategy(SPtr<BCStrategy>(new NonEqDensityBCStrategy()));
      //densityBC->setBCStrategy(SPtr<BCStrategy>(new NonReflectingOutflowBCStrategy()));


      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBC);
      bcVisitor.addBC(slipBC);
      bcVisitor.addBC(velocityBC);
      bcVisitor.addBC(densityBC);
      
      SPtr<BCSet> bcProc;
      bcProc = SPtr<BCSet>(new BCSet());

      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CumulantLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new RheologyK17LBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new HerschelBulkleyModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new BinghamModelLBMKernel());
      kernel->setBCSet(bcProc);
      //kernel->setForcingX1(forcing);
      //kernel->setWithForcing(true);

      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setPeriodicX1(false);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(true);
      grid->setDeltaX(deltax);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), outputPath + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

      //sphere
      SPtr<GbObject3D> sphere(new GbSphere3D(sphereCenter[0], sphereCenter[1], sphereCenter[2], radius));
      GbSystem3D::writeGeoObject(sphere.get(), outputPath + "/geo/sphere", WbWriterVtkXmlBinary::getInstance());
      SPtr<D3Q27Interactor> sphereInt(new D3Q27Interactor(sphere, grid, noSlipBC, Interactor3D::SOLID));

      ////////////////////////////////////////////
      //METIS
      SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, DIR_MMM, MetisPartitioner::KWAY));
      ////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIOMigrationSimulationObserver> restartSimulationObserver(new MPIIOMigrationSimulationObserver(grid, mSch, metisVisitor, outputPath, comm));
      restartSimulationObserver->setLBMKernel(kernel);
      restartSimulationObserver->setBCSet(bcProc);
      //restartSimulationObserver->setNu(k);
      //////////////////////////////////////////////////////////////////////////

      if (myid == 0)
      {
         UBLOG(logINFO, "Parameters:");
         //UBLOG(logINFO, "forcing = " << forcing);
         UBLOG(logINFO, "rho = " << rhoLB);
         UBLOG(logINFO, "U = " << U);
         UBLOG(logINFO, "Re = " << (U * d) / (k * std::pow(Gamma, n - 1)));
         UBLOG(logINFO, "Bn = " << tau0 /(k * std::pow(Gamma, n)));
         UBLOG(logINFO, "k = " << k);
         UBLOG(logINFO, "n = " << n);
         UBLOG(logINFO, "tau0 = " << tau0);
         UBLOG(logINFO, "deltax = " << deltax);
         UBLOG(logINFO, "number of levels = " << refineLevel + 1);
         UBLOG(logINFO, "number of threads = " << numOfThreads);
         UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
         UBLOG(logINFO, "blocknx = " << blocknx[0] << " " << blocknx[1] << " " << blocknx[2]);
         UBLOG(logINFO, "boundingBox = " << boundingBox[0] << " " << boundingBox[1] << " " << boundingBox[2]);
         UBLOG(logINFO, "sphereCenter = " << sphereCenter[0] << " " << sphereCenter[1] << " " << sphereCenter[2]);
         UBLOG(logINFO, "output path = " << outputPath);
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
         SPtr<D3Q27Interactor> wallZminInt(new D3Q27Interactor(wallZmin, grid, slipBC, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> wallZmaxInt(new D3Q27Interactor(wallZmax, grid, slipBC, Interactor3D::SOLID));
                                                                               
         SPtr<D3Q27Interactor> wallYminInt(new D3Q27Interactor(wallYmin, grid, slipBC, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> wallYmaxInt(new D3Q27Interactor(wallYmax, grid, slipBC, Interactor3D::SOLID));

         SPtr<D3Q27Interactor> wallXminInt(new D3Q27Interactor(wallXmin, grid, velocityBC, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> wallXmaxInt(new D3Q27Interactor(wallXmax, grid, densityBC, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, DIR_MMM, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(wallXminInt);
         intHelper.addInteractor(wallXmaxInt);
         intHelper.addInteractor(wallYminInt);
         intHelper.addInteractor(wallYmaxInt);
         intHelper.addInteractor(wallZminInt);
         intHelper.addInteractor(wallZmaxInt);
         intHelper.addInteractor(sphereInt);
         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), outputPath, WbWriterVtkXmlBinary::getInstance(), comm));
         ppblocks->process(0);

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nod = nob * (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nodg = nob * (blocknx[0] + gl) * (blocknx[1] + gl) * (blocknx[1] + gl);
         real needMemAll = real(nodg * (27 * sizeof(real) + sizeof(int) + sizeof(float) * 4));
         real needMem = needMemAll / real(comm->getNumberOfProcesses());

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
         restartSimulationObserver->restart((int)restartStep);
         grid->setTimeStep(restartStep);
         SetBcBlocksBlockVisitor v(sphereInt);
         grid->accept(v);
         sphereInt->initInteractor();
      }
      
      omp_set_num_threads(numOfThreads);

      //set connectors
      //InterpolationProcessorPtr iProcessor(new ThixotropyInterpolationProcessor());
      //static_pointer_cast<ThixotropyInterpolationProcessor>(iProcessor)->setOmegaMin(thix->getOmegaMin());
      //SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, k, iProcessor);
      //grid->accept(setConnsVisitor);

      OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
      grid->accept(setConnsVisitor);

      grid->accept(bcVisitor);

      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      WriteBoundaryConditionsSimulationObserver ppgeo = WriteBoundaryConditionsSimulationObserver(grid, geoSch, outputPath, WbWriterVtkXmlASCII::getInstance(), comm);
      ppgeo.process(0);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      SPtr<SimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

      //write data for visualization of macroscopic quantities
      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      //SPtr<UbScheduler> visSch(new UbScheduler(10,1));
      SPtr<WriteMacroscopicQuantitiesSimulationObserver> writeMQSimulationObserver(new WriteMacroscopicQuantitiesSimulationObserver(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
      //writeMQSimulationObserver->process(0);

      real area = UbMath::PI*radius*radius;
      SPtr<UbScheduler> forceSch(new UbScheduler(100));
      SPtr<CalculateForcesSimulationObserver> fp = make_shared<CalculateForcesSimulationObserver>(grid, forceSch, outputPath + "/forces/forces.txt", comm, velocity, area);
      fp->addInteractor(sphereInt);

      SPtr<WriteThixotropyQuantitiesSimulationObserver> writeThixotropicMQSimulationObserver(new WriteThixotropyQuantitiesSimulationObserver(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
      calculator->addSimulationObserver(npr);
      calculator->addSimulationObserver(fp);
      calculator->addSimulationObserver(writeMQSimulationObserver);
      calculator->addSimulationObserver(writeThixotropicMQSimulationObserver);
      calculator->addSimulationObserver(restartSimulationObserver);

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

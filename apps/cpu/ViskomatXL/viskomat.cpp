#include <iostream>
#include <string>

#include <VirtualFluids.h>

using namespace std;


void bflow(string configname)
{
   try
   {
      vf::basics::ConfigurationFile   config;
      config.load(configname);

      string          outputPath = config.getValue<string>("outputPath");
      string          geoPath = config.getValue<string>("geoPath");
      string          geoFile = config.getValue<string>("geoFile");
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
      double          cpStep = config.getValue<double>("cpStep");
      double          cpStart = config.getValue<double>("cpStart");
      bool            newStart = config.getValue<bool>("newStart");
      double          OmegaLB = config.getValue<double>("OmegaLB");
      double          tau0 = config.getValue<double>("tau0");
      double          scaleFactor = config.getValue<double>("scaleFactor");
      double          resolution = config.getValue<double>("resolution");

      vf::basics::ConfigurationFile   viscosity;
      //viscosity.load(viscosityPath + "/viscosity.cfg");
      //double nuLB = viscosity.getValue<double>("nuLB");

      //outputPath = outputPath + "/rheometerBingham_" + config.getValue<string>("resolution") + "_" + config.getValue<string>("OmegaLB");

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

      LBMReal rhoLB = 0.0;

      //akoustic
       //OmegaLB /= scaleFactor;
       //nuLB *=scaleFactor;
       //endTime *= scaleFactor;
       ////outTime = endTime;
       //cpStart = endTime;
       //cpStep  = endTime;

//diffusive
      //OmegaLB /= scaleFactor * scaleFactor;
      //tau0 /= scaleFactor * scaleFactor;
      //endTime *= scaleFactor * scaleFactor;
      //outTime = endTime;
      //cpStart = endTime;
      //cpStep = endTime;

      //double Re = 1.38230076758;
      double N  = 80; //rpm
      double Omega = 2 * UbMath::PI / 60.0 * N; //rad/s
      double mu    = 1; //Pa s
      double R     = 0.165 / 2.0; //m
      double rho   = 970; //kg/m^3
      double Re    = Omega * R * R * rho / mu;

      double nuLB = OmegaLB * R * 1e3 * R * 1e3 / Re;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());
      // double uWorld = (N * PI) / 30.0; //0.0037699111843
      // double rhoWorld = 2350.0; //kg/m^3
      //double R0 = boundingBox[0] * 0.5;

      //SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter(deltax, uWorld*R0, rhoWorld, 1.0, uLB));
      //if (myid == 0) UBLOG(logINFO, conv->toString());

      //bounding box

      double g_minX1 = boundingBox[0];
      double g_maxX1 = boundingBox[1];

      double g_minX2 = boundingBox[2];
      double g_maxX2 = boundingBox[3];
      
      double g_minX3 = boundingBox[4];
      double g_maxX3 = boundingBox[5];

      //double g_minX1 = -boundingBox[0]/2.0;
      //double g_minX2 = -boundingBox[1] / 2.0;
      //double g_minX3 = -boundingBox[2]/2.0;

      //double g_maxX1 = boundingBox[0]/2.0;
      //double g_maxX2 = boundingBox[1]/2.0;
      //double g_maxX3 = boundingBox[2]/2.0;

//      double blockLength = 3.0 * deltax;

      // double d = 2.0 * radius;
      // double U = uLB;
      // double Gamma = U / d;

      // double muWorld = 20; //Pa*s
      // double k = 0.0015; // muWorld / rhoWorld * conv->getFactorViscosityWToLb(); //(U * d) / (Re);

      // //double k = (U * d) / (Re * std::pow(Gamma, n - 1));
      // double yielStressWorld = 20; //Pa
      // double tau0 = 1e-6;// 3e-6;//yielStressWorld * conv->getFactorPressureWToLb(); //Bn * k * std::pow(Gamma, n);

      //double k = 0.05; // (U * d) / (Re * std::pow(Gamma, n - 1));
      //double tau0 = 3e-6; //Bn * k * std::pow(Gamma, n);

      //double forcing = 8e-7;

      //double omegaMin = 1.0e-8;

      SPtr<Rheology> thix = Rheology::getInstance();
      //thix->setPowerIndex(n);
      //thix->setViscosityParameter(k);
      thix->setYieldStress(tau0);
      //thix->setOmegaMin(omegaMin);

      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
      //noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new RheologyHerschelBulkleyModelNoSlipBCAlgorithm()));
      //noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new RheologyBinghamModelNoSlipBCAlgorithm()));

      SPtr<BCAdapter> slipBCAdapter(new SlipBCAdapter());
      slipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SimpleSlipBCAlgorithm()));

      //// rotation around X-axis
      mu::Parser fctVy;
      fctVy.SetExpr("-Omega*(x3-z0-r)");
      fctVy.DefineConst("Omega", OmegaLB);
      fctVy.DefineConst("r", 0.5 * (g_maxX3 - g_minX3));
      fctVy.DefineConst("z0", g_minX3);

      mu::Parser fctVz;
      fctVz.SetExpr("Omega*(x2-y0-r)");
      fctVz.DefineConst("Omega", OmegaLB);
      fctVz.DefineConst("r", 0.5 * (g_maxX2 - g_minX2));
      fctVz.DefineConst("y0", g_minX2);

      mu::Parser fctVx;
      fctVx.SetExpr("0.0");

      // rotation around Y-axis
      //mu::Parser fctVz;
      //// fctVx.SetExpr("omega*(r-x2)");
      //fctVz.SetExpr("Omega*(x1-r)");
      //fctVz.DefineConst("Omega", OmegaLB);
      //fctVz.DefineConst("r", 0.5 * (g_maxX1 - g_minX1));

      //mu::Parser fctVx;
      //fctVx.SetExpr("-Omega*(x3-r)");
      //fctVx.DefineConst("Omega", OmegaLB);
      //fctVx.DefineConst("r", 0.5 * (g_maxX1 - g_minX1));

      //mu::Parser fctVy;
      //fctVy.SetExpr("0.0");

      SPtr<BCAdapter> velocityBCAdapter(new VelocityBCAdapter(true, true, true, fctVx, fctVy, fctVz, 0, BCFunction::INFCONST));
      velocityBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));
      //velocityBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SimpleVelocityBCAlgorithm()));
      //velocityBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));
      //velocityBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new RheologyBinghamModelVelocityBCAlgorithm()));

      //SPtr<BCAdapter> densityBCAdapter(new DensityBCAdapter());
      //densityBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));
      ////densityBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));


      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(slipBCAdapter);
      bcVisitor.addBC(velocityBCAdapter);
      //bcVisitor.addBC(densityBCAdapter);
      
      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new BCProcessor());

      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new BGKLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CumulantLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel()); 
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CumulantK17LBMKernel()); 
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new RheologyBinghamModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new HerschelBulkleyModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new BinghamModelLBMKernel());
      kernel->setBCProcessor(bcProc);
      //kernel->setForcingX1(forcing);
      //kernel->setWithForcing(true);

      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setPeriodicX1(false);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(false);
      grid->setDeltaX(deltax);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), outputPath + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIOMigrationCoProcessor> restartCoProcessor(new MPIIOMigrationCoProcessor(grid, mSch, outputPath, comm));
      //SPtr<MPIIORestartCoProcessor> restartCoProcessor(new MPIIORestartCoProcessor(grid, mSch, outputPath, comm));
      restartCoProcessor->setLBMKernel(kernel);
      restartCoProcessor->setBCProcessor(bcProc);
      //restartCoProcessor->setNu(k);
      //////////////////////////////////////////////////////////////////////////

      ////stator
      // rotation around X-axis 
      SPtr<GbObject3D> stator(new GbCylinder3D(g_minX1 - 3.0 * deltax, g_minX2 + 0.5 * (g_maxX2 - g_minX2),
                                               g_minX3 + 0.5 * (g_maxX3 - g_minX3), g_maxX1 + 3.0 * deltax,
          g_minX2 + 0.5 * (g_maxX2 - g_minX2), g_minX3 + 0.5 * (g_maxX3 - g_minX3), 0.5 * (g_maxX3 - g_minX3) * 0.5));

       // rotation around Y-axis 
      //SPtr<GbObject3D> stator(new GbCylinder3D(g_minX1 + 0.5 * (g_maxX1 - g_minX1), g_minX2 - 3.0 * deltax, 
      //                                         g_minX3 + 0.5 * (g_maxX3 - g_minX3), g_minX1 + 0.5 * (g_maxX1 - g_minX1),
      //                                         g_maxX2 + 3.0 * deltax, g_minX3 + 0.5 * (g_maxX3 - g_minX3),
      //                                         0.5 * (g_maxX3 - g_minX3) * 0.5));

      SPtr<D3Q27Interactor> statorInt =
          SPtr<D3Q27Interactor>(new D3Q27Interactor(stator, grid, noSlipBCAdapter, Interactor3D::SOLID));
      
      //SPtr<GbTriFaceMesh3D> stator = make_shared<GbTriFaceMesh3D>();
      //stator->readMeshFromSTLFileBinary(geoPath + "/" + geoFile, false);
      //stator->translate(4.0, -73.0, -6.0);
      GbSystem3D::writeGeoObject(stator.get(), outputPath + "/geo/stator", WbWriterVtkXmlBinary::getInstance());
      
      //SPtr<D3Q27Interactor> statorInt = SPtr<D3Q27TriFaceMeshInteractor>(
      //    new D3Q27TriFaceMeshInteractor(stator, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));

      ////rotor (cylinder)
      // rotation around X-axis 
      SPtr<GbObject3D> rotor(new GbCylinder3D(
          g_minX1 - 3.0 * deltax, g_minX2 + 0.5 * (g_maxX2 - g_minX2),
                                              g_minX3 + 0.5 * (g_maxX3 - g_minX3), g_maxX1 + 3.0 * deltax,
          g_minX2 + 0.5 * (g_maxX2 - g_minX2), g_minX3 + 0.5 * (g_maxX3 - g_minX3), 0.5 * (g_maxX3 - g_minX3)));
      // rotation around Y-axis
      //SPtr<GbObject3D> rotor(new GbCylinder3D(g_minX1 + 0.5 * (g_maxX1 - g_minX1), g_minX2 - 3.0 * deltax,
      //                                        g_minX3 + 0.5 * (g_maxX3 - g_minX3), g_minX1 + 0.5 * (g_maxX1 - g_minX1),
      //                                        g_maxX2 + 3.0 * deltax, g_minX3 + 0.5 * (g_maxX3 - g_minX3),
      //                                        0.5 * (g_maxX3 - g_minX3)));

      GbSystem3D::writeGeoObject(rotor.get(), outputPath + "/geo/rotor", WbWriterVtkXmlBinary::getInstance());

      SPtr<D3Q27Interactor> rotorInt =
          SPtr<D3Q27Interactor>(new D3Q27Interactor(rotor, grid, velocityBCAdapter, Interactor3D::INVERSESOLID));

      //walls
      GbCuboid3DPtr wallXmin(new GbCuboid3D(g_minX1 - deltax, g_minX2 - deltax, g_minX3 - deltax, g_minX1,
          g_maxX2 + deltax, g_maxX3 + deltax));
      if (myid == 0) GbSystem3D::writeGeoObject(wallXmin.get(), outputPath + "/geo/wallXmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr wallXmax(new GbCuboid3D(g_maxX1, g_minX2 - deltax, g_minX3 - deltax, g_maxX1 +  (double)blocknx[0]*deltax,
          g_maxX2 + deltax, g_maxX3 + deltax));
      if (myid == 0) GbSystem3D::writeGeoObject(wallXmax.get(), outputPath + "/geo/wallXmax", WbWriterVtkXmlASCII::getInstance());

      //wall interactors
      SPtr<D3Q27Interactor> wallXminInt(new D3Q27Interactor(wallXmin, grid, slipBCAdapter, Interactor3D::SOLID));
      SPtr<D3Q27Interactor> wallXmaxInt(new D3Q27Interactor(wallXmax, grid, slipBCAdapter, Interactor3D::SOLID));

      if (myid == 0)
      {
         UBLOG(logINFO, "Parameters:");
         //UBLOG(logINFO, "forcing = " << forcing);
         UBLOG(logINFO, "N = " << N << " rpm");
         UBLOG(logINFO, "Omega = " << Omega << " rad/s");
         UBLOG(logINFO, "Re = " << Re);
         UBLOG(logINFO, "rho = " << rhoLB);
         UBLOG(logINFO, "uLB = " << OmegaLB);
         UBLOG(logINFO, "nuLB = " << nuLB);
         // UBLOG(logINFO, "Re = " << (U * d) / (k * std::pow(Gamma, n - 1)));
         // UBLOG(logINFO, "Bn = " << tau0 /(k * std::pow(Gamma, n)));
         // UBLOG(logINFO, "k = " << k);
         // UBLOG(logINFO, "n = " << n);
         UBLOG(logINFO, "tau0 = " << tau0);
         UBLOG(logINFO, "scaleFactor = " << scaleFactor);
         UBLOG(logINFO, "deltax = " << deltax);
         UBLOG(logINFO, "number of levels = " << refineLevel + 1);
         UBLOG(logINFO, "number of threads = " << numOfThreads);
         UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
         UBLOG(logINFO, "blocknx = " << blocknx[0] << " " << blocknx[1] << " " << blocknx[2]);
         UBLOG(logINFO, "resolution = " << resolution);
         // UBLOG(logINFO, "boundingBox = " << boundingBox[0] << " " << boundingBox[1] << " " << boundingBox[2]);
         // UBLOG(logINFO, "sphereCenter = " << sphereCenter[0] << " " << sphereCenter[1] << " " << sphereCenter[2]);
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


         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(wallXminInt);
         intHelper.addInteractor(wallXmaxInt);
         intHelper.addInteractor(statorInt);
         intHelper.addInteractor(rotorInt);
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

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
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

         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         WriteBoundaryConditionsCoProcessor ppgeo = WriteBoundaryConditionsCoProcessor(grid, geoSch, outputPath, WbWriterVtkXmlBinary::getInstance(), comm);
         ppgeo.process(0);

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         restartCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);
         
         //SetBcBlocksBlockVisitor v1(wallXminInt);
         //grid->accept(v1);
         //wallXminInt->initInteractor();
         //
         //SetBcBlocksBlockVisitor v2(wallXmaxInt);
         //grid->accept(v2);
         //wallXmaxInt->initInteractor();
         
         SetBcBlocksBlockVisitor v3(statorInt);
         grid->accept(v3);
         statorInt->initInteractor();

         SetBcBlocksBlockVisitor v4(rotorInt);
         grid->accept(v4);
         rotorInt->initInteractor();


      }
      
      omp_set_num_threads(numOfThreads);

      //set connectors
      //InterpolationProcessorPtr iProcessor(new ThixotropyInterpolationProcessor());
      //static_pointer_cast<ThixotropyInterpolationProcessor>(iProcessor)->setOmegaMin(thix->getOmegaMin());
      //SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
      //grid->accept(setConnsVisitor);

      OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
      grid->accept(setConnsVisitor);

      grid->accept(bcVisitor);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      SPtr<CoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      //write data for visualization of macroscopic quantities
      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      //SPtr<UbScheduler> visSch(new UbScheduler(10,1));
      SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
      //writeMQCoProcessor->process(100);

      SPtr<UbScheduler> forceSch(new UbScheduler(100));
      SPtr<CalculateTorqueCoProcessor> fp = make_shared<CalculateTorqueCoProcessor>(grid, forceSch, outputPath + "/torque/TorqueRotor.csv", comm);
      fp->addInteractor(rotorInt);
      SPtr<CalculateTorqueCoProcessor> fp2 = make_shared<CalculateTorqueCoProcessor>(grid, forceSch, outputPath + "/torque/TorqueStator.csv", comm);
      fp2->addInteractor(statorInt);

      //SPtr<WriteThixotropyQuantitiesCoProcessor> writeThixotropicMQCoProcessor(new WriteThixotropyQuantitiesCoProcessor(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
      calculator->addCoProcessor(npr);
      calculator->addCoProcessor(fp);
      calculator->addCoProcessor(fp2);
      calculator->addCoProcessor(writeMQCoProcessor);
      //calculator->addCoProcessor(writeThixotropicMQCoProcessor);
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

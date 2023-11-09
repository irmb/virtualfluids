#include <iostream>
#include <string>

#include "VirtualFluids.h"
#include "NonNewtonianFluids/NonNewtonianFluids.h"

using namespace std;


void bflow(string configname)
{
    using namespace vf::lbm::dir;

   try
   {
      vf::basics::ConfigurationFile   config;
      config.load(configname);

      string          outputPath = config.getValue<string>("outputPath");
      string          geoPath = config.getValue<string>("geoPath");
      string          geoFile = config.getValue<string>("geoFile");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
      vector<real>  boundingBox = config.getVector<real>("boundingBox");
      real          endTime = config.getValue<real>("endTime");
      real          outTime = config.getValue<real>("outTime");
      real          availMem = config.getValue<real>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      bool            logToFile = config.getValue<bool>("logToFile");
      real          restartStep = config.getValue<real>("restartStep");
      real          deltax = config.getValue<real>("deltax");
      real          cpStep = config.getValue<real>("cpStep");
      real          cpStart = config.getValue<real>("cpStart");
      bool            newStart = config.getValue<bool>("newStart");
      real          OmegaLB = config.getValue<real>("OmegaLB");
      real          tau0 = config.getValue<real>("tau0");
      real          N = config.getValue<real>("N");
      real          mu = config.getValue<real>("mu");


      vf::basics::ConfigurationFile   viscosity;

      std::shared_ptr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
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

      //double N  = 70; //rpm
      real Omega = 2 * UbMath::PI / 60.0 * N; //rad/s
      //double mu    = 5; //Pa s
      real R     = 0.165 / 2.0; //m
      real rho = 2150;// 970; //kg/m^3
      real Re    = Omega * R * R * rho / mu;

      //double nuLB = OmegaLB * R * 1e3 * R * 1e3 / Re;

      real dx = deltax * 1e-3;
      real nuLB = OmegaLB * (R / dx)*(R / dx) / Re;

      real Bm = tau0/(mu*Omega);
      real tau0LB = Bm*nuLB*OmegaLB;


      //double dx = 1.0 * 1e-3;
      //double nuLB = OmegaLB * (R / dx)*(R / dx) / Re;

      //acustic scaling
      // OmegaLB /= 2.0;
      // nuLB    *= 2.0;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());
      //SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter(1, 1461, 970, 1e3));
      //UBLOG(logINFO, conv->toString());

      //bounding box

      real g_minX1 = boundingBox[0];
      real g_maxX1 = boundingBox[1];

      real g_minX2 = boundingBox[2];
      real g_maxX2 = boundingBox[3];
      
      real g_minX3 = boundingBox[4];
      real g_maxX3 = boundingBox[5];

      SPtr<Rheology> thix = Rheology::getInstance();
      //thix->setPowerIndex(n);
      //thix->setViscosityParameter(k);
      thix->setYieldStress(tau0LB);
      //thix->setOmegaMin(omegaMin);

      SPtr<BC> noSlipBC(new NoSlipBC());
      noSlipBC->setBCStrategy(SPtr<BCStrategy>(new NoSlipBCStrategy()));
      //noSlipBC->setBCStrategy(SPtr<BCStrategy>(new RheologyHerschelBulkleyModelNoSlipBCStrategy()));
      noSlipBC->setBCStrategy(SPtr<BCStrategy>(new RheologyBinghamModelNoSlipBCStrategy()));

      SPtr<BC> slipBC(new SlipBC());
      slipBC->setBCStrategy(SPtr<BCStrategy>(new SimpleSlipBCStrategy()));
      //slipBC->setBCStrategy(SPtr<BCStrategy>(new SlipBCStrategy()));

      //// rotation around X-axis
      mu::Parser fctVy;
      fctVy.SetExpr("-Omega*(x3-z0-r)/deltax");
      fctVy.DefineConst("Omega", OmegaLB);
      fctVy.DefineConst("r", 0.5 * (g_maxX3 - g_minX3));
      fctVy.DefineConst("z0", g_minX3);
      fctVy.DefineConst("deltax", deltax);

      mu::Parser fctVz;
      fctVz.SetExpr("Omega*(x2-y0-r)/deltax");
      fctVz.DefineConst("Omega", OmegaLB);
      fctVz.DefineConst("r", 0.5 * (g_maxX2 - g_minX2));
      fctVz.DefineConst("y0", g_minX2);
      fctVz.DefineConst("deltax", deltax);

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

      SPtr<BC> velocityBC(new VelocityBC(true, true, true, fctVx, fctVy, fctVz, 0, BCFunction::INFCONST));
      //velocityBC->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));
      velocityBC->setBCStrategy(SPtr<BCStrategy>(new SimpleVelocityBCStrategy()));
      //velocityBC->setBCStrategy(SPtr<BCStrategy>(new VelocityWithDensityBCStrategy()));
      //velocityBC->setBCStrategy(SPtr<BCStrategy>(new RheologyBinghamModelVelocityBCStrategy()));

      //SPtr<BC> densityBC(new DensityBC());
      //densityBC->setBCStrategy(SPtr<BCStrategy>(new NonEqDensityBCStrategy()));
      ////densityBC->setBCStrategy(SPtr<BCStrategy>(new NonReflectingOutflowBCStrategy()));


      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      //bcVisitor.addBC(noSlipBC);
      bcVisitor.addBC(slipBC);
      bcVisitor.addBC(velocityBC);
      //bcVisitor.addBC(densityBC);
      
      SPtr<BCSet> bcProc;
      bcProc = SPtr<BCSet>(new BCSet());

      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new BGKLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CumulantLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel()); 
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CumulantK17LBMKernel()); 
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new RheologyBinghamModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new HerschelBulkleyModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new BinghamModelLBMKernel());
      kernel->setBCSet(bcProc);
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

      ////////////////////////////////////////////
      //METIS
      SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, dMMM, MetisPartitioner::RECURSIVE));
      ////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIOMigrationSimulationObserver> restartSimulationObserver(new MPIIOMigrationSimulationObserver(grid, mSch, metisVisitor, outputPath, comm));
      //SPtr<MPIIORestartSimulationObserver> restartSimulationObserver(new MPIIORestartSimulationObserver(grid, mSch, outputPath, comm));
      restartSimulationObserver->setLBMKernel(kernel);
      restartSimulationObserver->setBCSet(bcProc);
      //restartSimulationObserver->setNu(k);
      //////////////////////////////////////////////////////////////////////////

      ////stator
      // rotation around X-axis 
       //SPtr<GbObject3D> stator(new GbCylinder3D(g_minX1 - 3.0 * deltax, g_minX2 + 0.5 * (g_maxX2 - g_minX2),
       //                                         g_minX3 + 0.5 * (g_maxX3 - g_minX3), g_maxX1 + 3.0 * deltax,
       //    g_minX2 + 0.5 * (g_maxX2 - g_minX2), g_minX3 + 0.5 * (g_maxX3 - g_minX3), 0.5 * (g_maxX3 - g_minX3) * 0.5));

      // SPtr<GbObject3D> stator(new GbCylinder3D(g_minX1 - 4.0 * deltax, g_minX2 + 0.5 * (g_maxX2 - g_minX2),
      //                                          g_minX3 + 0.5 * (g_maxX3 - g_minX3), g_maxX1 + 3.0 * deltax,
      //     g_minX2 + 0.5 * (g_maxX2 - g_minX2), g_minX3 + 0.5 * (g_maxX3 - g_minX3), 12.0*0.5));

      ////  // rotation around Y-axis 
      //// //SPtr<GbObject3D> stator(new GbCylinder3D(g_minX1 + 0.5 * (g_maxX1 - g_minX1), g_minX2 - 3.0 * deltax, 
      //// //                                         g_minX3 + 0.5 * (g_maxX3 - g_minX3), g_minX1 + 0.5 * (g_maxX1 - g_minX1),
      //// //                                         g_maxX2 + 3.0 * deltax, g_minX3 + 0.5 * (g_maxX3 - g_minX3),
      //// //                                         0.5 * (g_maxX3 - g_minX3) * 0.5));

      // SPtr<D3Q27Interactor> statorInt =
      //    SPtr<D3Q27Interactor>(new D3Q27Interactor(stator, grid, noSlipBC, Interactor3D::SOLID));
      
      SPtr<GbTriFaceMesh3D> stator = make_shared<GbTriFaceMesh3D>();
      stator->readMeshFromSTLFileBinary(geoPath + "/" + geoFile, false);
      //stator->scale(2.0, 2.0, 2.0);
      //stator->translate(16.0, 0.0, 0.0);
      //stator->translate(4.0, -73.0, -6.0);

      SPtr<D3Q27Interactor> statorInt = SPtr<D3Q27TriFaceMeshInteractor>(
         new D3Q27TriFaceMeshInteractor(stator, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::EDGES));

      GbSystem3D::writeGeoObject(stator.get(), outputPath + "/geo/stator", WbWriterVtkXmlBinary::getInstance());

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
          SPtr<D3Q27Interactor>(new D3Q27Interactor(rotor, grid, velocityBC, Interactor3D::INVERSESOLID));

      //walls
      GbCuboid3DPtr wallXmin(new GbCuboid3D(g_minX1 - deltax, g_minX2 - deltax, g_minX3 - deltax, g_minX1,
          g_maxX2 + deltax, g_maxX3 + deltax));
      if (myid == 0) GbSystem3D::writeGeoObject(wallXmin.get(), outputPath + "/geo/wallXmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr wallXmax(new GbCuboid3D(g_maxX1, g_minX2 - deltax, g_minX3 - deltax, g_maxX1 +  (real)blocknx[0]*deltax,
          g_maxX2 + deltax, g_maxX3 + deltax));
      if (myid == 0) GbSystem3D::writeGeoObject(wallXmax.get(), outputPath + "/geo/wallXmax", WbWriterVtkXmlASCII::getInstance());

      //wall interactors
      SPtr<D3Q27Interactor> wallXminInt(new D3Q27Interactor(wallXmin, grid, slipBC, Interactor3D::SOLID));
      SPtr<D3Q27Interactor> wallXmaxInt(new D3Q27Interactor(wallXmax, grid, slipBC, Interactor3D::SOLID));

      if (myid == 0)
      {
         UBLOG(logINFO, "Parameters:");
         UBLOG(logINFO, "N = " << N << " rpm");
         UBLOG(logINFO, "Omega = " << Omega << " rad/s");
         UBLOG(logINFO, "mu = " << mu << " Pa s");
         UBLOG(logINFO, "tau0 = " << tau0<< " Pa");
         UBLOG(logINFO, "rho = " << rho<< " kg/m^3");
         UBLOG(logINFO, "Re = " << Re);
         UBLOG(logINFO, "Bm = " << Bm);
         UBLOG(logINFO, "rhoLB = " << rhoLB);
         UBLOG(logINFO, "uLB = " << OmegaLB);
         UBLOG(logINFO, "nuLB = " << nuLB);
         UBLOG(logINFO, "tau0LB = " << tau0LB);
         UBLOG(logINFO, "deltax = " << deltax << " mm");
         UBLOG(logINFO, "number of levels = " << refineLevel + 1);
         UBLOG(logINFO, "number of threads = " << numOfThreads);
         UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
         UBLOG(logINFO, "blocknx = " << blocknx[0] << " " << blocknx[1] << " " << blocknx[2]);
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



         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(wallXmaxInt);
         intHelper.addInteractor(statorInt);
         intHelper.addInteractor(rotorInt);
         intHelper.addInteractor(wallXminInt);
         
         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), outputPath, WbWriterVtkXmlBinary::getInstance(), comm));
         ppblocks->update(0);

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
         WriteBoundaryConditionsSimulationObserver ppgeo = WriteBoundaryConditionsSimulationObserver(grid, geoSch, outputPath, WbWriterVtkXmlBinary::getInstance(), comm);
         ppgeo.update(0);

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         restartSimulationObserver->restart((int)restartStep);
         
         //restartSimulationObserver->readBlocks((int)restartStep);
         //restartSimulationObserver->readDataSet((int)restartStep);
         ////restartSimulationObserver->readBoundaryConds((int)restartStep);
         //grid->setTimeStep((int)restartStep);
         
         SetBcBlocksBlockVisitor v2(wallXmaxInt);
         grid->accept(v2);
         wallXmaxInt->initInteractor();

         SetBcBlocksBlockVisitor v3(statorInt);
         grid->accept(v3);
         statorInt->initInteractor();

         SetBcBlocksBlockVisitor v4(rotorInt);
         grid->accept(v4);
         rotorInt->initInteractor();

         SetBcBlocksBlockVisitor v1(wallXminInt);
         grid->accept(v1);
         wallXminInt->initInteractor();

         SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), outputPath,
                                                               WbWriterVtkXmlBinary::getInstance(), comm));
         ppblocks->update(1);
      }
      
      //omp_set_num_threads(numOfThreads);

      //set connectors
      //InterpolationProcessorPtr iProcessor(new ThixotropyInterpolationProcessor());
      //static_pointer_cast<ThixotropyInterpolationProcessor>(iProcessor)->setOmegaMin(thix->getOmegaMin());
      //SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
      //grid->accept(setConnsVisitor);

      OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
      grid->accept(setConnsVisitor);

      grid->accept(bcVisitor);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      SPtr<SimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

      //write data for visualization of macroscopic quantities
      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      //SPtr<UbScheduler> visSch(new UbScheduler(10,1));
      SPtr<WriteMacroscopicQuantitiesSimulationObserver> writeMQSimulationObserver(new WriteMacroscopicQuantitiesSimulationObserver(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
      //writeMQSimulationObserver->update(100);

      SPtr<UbScheduler> forceSch(new UbScheduler(100));
      SPtr<CalculateTorqueSimulationObserver> fp = make_shared<CalculateTorqueSimulationObserver>(grid, forceSch, outputPath + "/torque/TorqueRotor.csv", comm);
      fp->addInteractor(rotorInt);
      SPtr<CalculateTorqueSimulationObserver> fp2 = make_shared<CalculateTorqueSimulationObserver>(grid, forceSch, outputPath + "/torque/TorqueStator.csv", comm);
      fp2->addInteractor(statorInt);

      //SPtr<WriteThixotropyQuantitiesSimulationObserver> writeThixotropicMQSimulationObserver(new WriteThixotropyQuantitiesSimulationObserver(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Simulation> simulation(new Simulation(grid, stepGhostLayer, endTime));
      simulation->addSimulationObserver(npr);
      //simulation->addSimulationObserver(fp);
      simulation->addSimulationObserver(fp2);
      //simulation->addSimulationObserver(writeMQSimulationObserver);
      //simulation->addSimulationObserver(writeThixotropicMQSimulationObserver);
      simulation->addSimulationObserver(restartSimulationObserver);

      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      simulation->run();
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

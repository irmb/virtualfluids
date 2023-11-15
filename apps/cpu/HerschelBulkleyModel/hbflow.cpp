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

      string          pathname = config.getValue<string>("pathname");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
      vector<real>  boundingBox = config.getVector<real>("boundingBox");
      real          nuLB = config.getValue<real>("nuLB");
      real          endTime = config.getValue<real>("endTime");
      real          outTime = config.getValue<real>("outTime");
      real          availMem = config.getValue<real>("availMem");
      //int             refineLevel = config.getValue<int>("refineLevel");
      bool            logToFile = config.getValue<bool>("logToFile");
      //double          restartStep = config.getValue<double>("restartStep");
      real          deltax = config.getValue<real>("deltax");
      //double          cpStep = config.getValue<double>("cpStep");
      //double          cpStepStart = config.getValue<double>("cpStepStart");
      //bool            newStart = config.getValue<bool>("newStart");
      real          forcing = config.getValue<real>("forcing");
      //double          n = config.getValue<double>("n");
      //double          k = config.getValue<double>("k");
      real          tau0 = config.getValue<real>("tau0");
      real          velocity = config.getValue<real>("velocity");
      real          n = config.getValue<real>("n");
//      double          Re = config.getValue<double>("Re");
//      double          Bn = config.getValue<double>("Bn");
      real          scaleFactor = config.getValue<real>("scaleFactor");

      SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
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

      real rhoLB = 0.0;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      //bounding box
      //double g_minX1 = -0.707107+1.0;
      //double g_minX2 = -1.41421+1.0;
      //double g_minX3 = 1.0;

      //double g_maxX1 = boundingBox[0];
      //double g_maxX2 = boundingBox[1];
      //double g_maxX3 = boundingBox[2]+1.0;

      real g_minX1 = 0.0;
      real g_minX2 = -boundingBox[1]/2.0;
      real g_minX3 = -boundingBox[2]/2.0;

      real g_maxX1 = boundingBox[0];
      real g_maxX2 = boundingBox[1]/2.0;
      real g_maxX3 = boundingBox[2]/2.0;

      

      real blockLength = 3.0 * deltax;

//      double h = (g_maxX2) / 2.0;
//      double dex = g_maxX1;

      //LBMReal tau0 = 1.2e-7;
      //LBMReal k = nuLB;
      //LBMReal n = 0.4;


      real d = boundingBox[1];
      real U = velocity;
      real Gamma = U / d;

      //double scaleFactor = 2.0;

      // Diffusive Scaling

      //double k = 0.005; // (U * d) / (Re * std::pow(Gamma, n - 1));
      //double tau0 = 3e-5 / (scaleFactor * scaleFactor);// *4.0);//*4.0);// Bn* k* std::pow(Gamma, n);
      //forcing /= scaleFactor * scaleFactor * scaleFactor; 
      //endTime *= scaleFactor * scaleFactor;
      //deltax /= scaleFactor;

      // Acoustic Scaling

      real k = nuLB * scaleFactor;
      //double tau0 = 3e-5; 
      forcing /= scaleFactor;
      endTime *= scaleFactor;
      deltax /= scaleFactor;

      //outTime = endTime;

      real beta = 14;
      real c = 10; // 1.0 / 6.0;
      real mu0 = 1e-4;

      SPtr<Rheology> thix = Rheology::getInstance();
      //Herschel-Bulkley
      thix->setPowerIndex(n);
      thix->setViscosityParameter(k);
      thix->setYieldStress(tau0);
      //Powell-Eyring
      thix->setBeta(beta);
      thix->setC(c);
      thix->setMu0(mu0);

      SPtr<BC> noSlipBC(new NoSlipBC());
      //noSlipBC->setBCStrategy(SPtr<BCStrategy>(new RheologyHerschelBulkleyModelNoSlipBCStrategy()));
      //noSlipBC->setBCStrategy(SPtr<BCStrategy>(new RheologyPowellEyringModelNoSlipBCStrategy()));
      //noSlipBC->setBCStrategy(SPtr<BCStrategy>(new RheologyBinghamModelNoSlipBCStrategy()));

      mu::Parser fctVx;
      fctVx.SetExpr("u");
      fctVx.DefineConst("u", 0.001);
 

      SPtr<BC> velocityBC(new VelocityBC(true, false, false, fctVx, 0, BCFunction::INFCONST));
      //velocityBC->setBCStrategy(SPtr<BCStrategy>(new VelocityInterpolated()));
      velocityBC->setBCStrategy(SPtr<BCStrategy>(new RheologyBinghamModelVelocityBCStrategy()));

      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      //bcVisitor.addBC(noSlipBC);
      //bcVisitor.addBC(velocityBC);

      SPtr<BCSet> bcProc;
      bcProc = SPtr<BCSet>(new BCSet());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new PowellEyringModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new HerschelBulkleyModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new RheologyK17LBMKernel());
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new RheologyBinghamModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new BinghamModelLBMKernel());
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      
      //double forcingXY = forcing / sqrt(2.0);
      //kernel->setForcingX1(forcingXY);
      //kernel->setForcingX2(forcingXY);
      
      kernel->setForcingX1(forcing);
      kernel->setWithForcing(true);
      kernel->setBCSet(bcProc);


      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(true);
      grid->setDeltaX(deltax);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      if (myid == 0)
      {
         UBLOG(logINFO, "Parameters:");
         UBLOG(logINFO, "forcing = " << forcing);
         UBLOG(logINFO, "rho = " << rhoLB);
         //UBLOG(logINFO, "nu = " << nuLB);
         UBLOG(logINFO, "U = " << U);
         UBLOG(logINFO, "Re = " << (U * d) / (k * std::pow(Gamma, n - 1)));
         UBLOG(logINFO, "Bn = " << tau0 / (k * std::pow(Gamma, n)));
         UBLOG(logINFO, "k = " << k);
         UBLOG(logINFO, "n = " << n);
         UBLOG(logINFO, "tau0 = " << tau0);
         UBLOG(logINFO, "deltax = " << deltax);
         //UBLOG(logINFO, "number of levels = " << refineLevel + 1);
         UBLOG(logINFO, "numOfThreads = " << numOfThreads);
         UBLOG(logINFO, "Preprozess - start");
      }

      //walls
      //GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_minX3));
      //if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname + "/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

      //GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_maxX3, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
      //if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname + "/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_minX2, g_maxX3 + blockLength));
      if (myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname + "/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1 - blockLength, g_maxX2, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
      if (myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname + "/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

      //GbCuboid3DPtr addWallXY(new GbCuboid3D(g_minX1 - 0.5*deltax, g_minX2 - 0.5 * deltax, g_minX3 - 0.5 * deltax, g_maxX1 + 0.5 * deltax, g_maxX2 + 0.5 * deltax, g_maxX3 + 0.5 * deltax));
      //if (myid == 0) GbSystem3D::writeGeoObject(addWallXY.get(), pathname + "/geo/addWallXY", WbWriterVtkXmlASCII::getInstance());

      //SPtr<GbTriFaceMesh3D> wall45(new GbTriFaceMesh3D());
      //wall45->readMeshFromSTLFile("d:/Projects/TRR277/Project/WP1/PF45/Geo/wall45_ASCII.stl", true);


      //wall interactors
      //SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, noSlipBC, Interactor3D::SOLID));
      //SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, noSlipBC, Interactor3D::SOLID));

      SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, noSlipBC, Interactor3D::SOLID));
      SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, noSlipBC, Interactor3D::SOLID));

      //SPtr<D3Q27TriFaceMeshInteractor> wall45Int(new D3Q27TriFaceMeshInteractor(wall45, grid, noSlipBC, Interactor3D::SOLID));

      ////////////////////////////////////////////
      //METIS
      SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, dMMM, MetisPartitioner::RECURSIVE));
      ////////////////////////////////////////////
      /////delete solid blocks
      if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
      InteractorsHelper intHelper(grid, metisVisitor);
      //intHelper.addInteractor(addWallZminInt);
      //intHelper.addInteractor(addWallZmaxInt);
      intHelper.addInteractor(addWallYminInt);
      intHelper.addInteractor(addWallYmaxInt);
      //intHelper.addInteractor(wall45Int);
      intHelper.selectBlocks();
      if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
      //////////////////////////////////////

      SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
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

      SetKernelBlockVisitor kernelVisitor(kernel, k, availMem, needMem);
      grid->accept(kernelVisitor);

      //BC
      intHelper.setBC();

      //initialization of distributions
      InitDistributionsBlockVisitor initVisitor;
      grid->accept(initVisitor);


      if (myid == 0) UBLOG(logINFO, "Preprozess - end");


      //set connectors
      //InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolator());
      //SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, k, iProcessor);
      //grid->accept(setConnsVisitor);

      OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
      grid->accept(setConnsVisitor);

      grid->accept(bcVisitor);

      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      WriteBoundaryConditionsSimulationObserver ppgeo = WriteBoundaryConditionsSimulationObserver(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm);
      ppgeo.update(1);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      SPtr<SimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> forceSch(new UbScheduler(1000));
      //real dummy = 1;
      SPtr<CalculateTorqueSimulationObserver> fp = std::make_shared<CalculateTorqueSimulationObserver>(grid, forceSch, pathname + "/forces/forces.csv", comm);
      fp->addInteractor(addWallYminInt);

      SPtr<CalculateTorqueSimulationObserver> fp2 = std::make_shared<CalculateTorqueSimulationObserver>(grid, forceSch, pathname + "/forces/forces2.csv", comm);
      fp2->addInteractor(addWallYmaxInt);

      //write data for visualization of macroscopic quantities
      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      SPtr<WriteMacroscopicQuantitiesSimulationObserver> writeMQSimulationObserver(new WriteMacroscopicQuantitiesSimulationObserver(grid, visSch, pathname,
         WbWriterVtkXmlASCII::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

      SPtr<WriteThixotropyQuantitiesSimulationObserver> writeThixotropicMQSimulationObserver(new WriteThixotropyQuantitiesSimulationObserver(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(outTime));
      SPtr<Simulation> simulation(new Simulation(grid, stepGhostLayer, endTime));
      simulation->addSimulationObserver(npr);
      simulation->addSimulationObserver(writeMQSimulationObserver);
      simulation->addSimulationObserver(writeThixotropicMQSimulationObserver);
      simulation->addSimulationObserver(fp);
      simulation->addSimulationObserver(fp2);
      //simulation->addSimulationObserver(migSimulationObserver);
      //simulation->addSimulationObserver(restartSimulationObserver);

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

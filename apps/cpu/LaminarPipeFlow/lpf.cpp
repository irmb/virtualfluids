#include <iostream>
#include <string>
//#include <omp.h>

#include "VirtualFluids.h"

using namespace std;


void run(string configname)
{
    using namespace vf::lbm::dir;

   try
   {
      vf::basics::ConfigurationFile   config;
      config.load(configname);

      string          pathname = config.getValue<string>("pathname");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
      real          uLB = config.getValue<real>("uLB");
      real          endTime = config.getValue<real>("endTime");
      real          outTime = config.getValue<real>("outTime");
      real          availMem = config.getValue<real>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      real          Re = config.getValue<real>("Re");
      real          dx = config.getValue<real>("dx");
      vector<real>  length = config.getVector<real>("length");
      bool            logToFile = config.getValue<bool>("logToFile");
      real          restartStep = config.getValue<real>("restartStep");
      real          cpStart = config.getValue<real>("cpStart");
      real          cpStep = config.getValue<real>("cpStep");
      bool            newStart = config.getValue<bool>("newStart");

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

            vf::logging::Logger::changeLogPath(pathname);
         }
      }

      //Sleep(30000);

      real dLB = length[1] / dx;
      real rhoLB = 0.0;
      real nuLB = (uLB*dLB) / Re;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      //BC Adapter
      //////////////////////////////////////////////////////////////////////////////
      SPtr<BC> noSlipBC(new NoSlipBC());
      //noSlipBC->setBCStrategy(SPtr<BCStrategy>(new ThinWallNoSlip()));
      //SPtr<BC> denBC(new PressureBC(rhoLB));
      //denBC->setBCStrategy(SPtr<BCStrategy>(new EqDensityBCStrategy()));

      noSlipBC->setBCStrategy(SPtr<BCStrategy>(new NoSlipInterpolated()));

      SPtr<BC> denBC(new PressureBC(rhoLB));
      denBC->setBCStrategy(SPtr<BCStrategy>(new OutflowNonReflecting()));
      //denBC->setBCStrategy(SPtr<BCStrategy>(new PressureNonEquilibrium()));

      //double startTime = 5;
      mu::Parser fct1;
      fct1.SetExpr("U");
      fct1.DefineConst("U", uLB);
      SPtr<BC> velBC1(new VelocityBC(true, false, false, fct1, 0, BCFunction::INFCONST));
      //velBC1->setBCStrategy(SPtr<BCStrategy>(new VelocityInterpolated()));
      velBC1->setBCStrategy(SPtr<BCStrategy>(new VelocityWithPressureInterpolated()));

      //mu::Parser fct2;
      //fct2.SetExpr("U");
      //fct2.DefineConst("U", uLB);
      //SPtr<BC> velBC2(new VelocityBC(true, false, false, fct2, startTime, BCFunction::INFCONST));

      //////////////////////////////////////////////////////////////////////////////////
      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;

      SPtr<Grid3D> grid(new Grid3D(comm));

      SPtr<BCSet> bcProc;
      bcProc = SPtr<BCSet>(new BCSet());

      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new K17CompressibleNavierStokes());
      //double bulckViscosity = 3700*nuLB;
      //dynamicPointerCast<K17CompressibleNavierStokes>(kernel)->setBulkViscosity(bulckViscosity);
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new K17CompressibleNavierStokes());
      kernel->setBCSet(bcProc);
      kernel->setBCSet(bcProc);

      //////////////////////////////////////////////////////////////////////////
      SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, d00M));
      //restart
      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      //SPtr<MPIIOMigrationSimulationObserver> migSimulationObserver(new MPIIOMigrationSimulationObserver(grid, mSch, metisVisitor, pathname + "/mig", comm));
      SPtr<MPIIOMigrationBESimulationObserver> migSimulationObserver(new MPIIOMigrationBESimulationObserver(grid, mSch, metisVisitor, pathname + "/mig", comm));
      migSimulationObserver->setLBMKernel(kernel);
      migSimulationObserver->setBCSet(bcProc);
      migSimulationObserver->setNu(nuLB);
      migSimulationObserver->setNuLG(0.01, 0.01);
      migSimulationObserver->setDensityRatio(1);
      //////////////////////////////////////////////////////////////////////////

      SPtr<D3Q27Interactor> inflowInt;

      if (newStart)
      {

         //bounding box
         real g_minX1 = 0.0;
         real g_minX2 = -length[1] / 2.0;
         real g_minX3 = -length[2] / 2.0;

         real g_maxX1 = length[0];
         real g_maxX2 = length[1] / 2.0;
         real g_maxX3 = length[2] / 2.0;

         //geometry
         //x
         //SPtr<GbObject3D> cylinder(new GbCylinder3D(g_minX1 - 2.0*dx, 0.0, 0.0, g_maxX1 + 2.0*dx, 0.0, 0.0, dLB / 2.0));
         //y
         SPtr<GbObject3D> cylinder(new GbCylinder3D(g_minX1 - 2.0*dx, 0.0, 0.0, g_maxX1 + 2.0*dx, 0.0, 0.0, dLB / 2.0));
         GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


         real blockLength = blocknx[0] * dx;



         if (myid == 0)
         {
            VF_LOG_INFO("uLb = {}", uLB);
            VF_LOG_INFO("rho = {}", rhoLB);
            VF_LOG_INFO("nuLb = {}", nuLB);
            VF_LOG_INFO("Re = {}", Re);
            VF_LOG_INFO("dx = {}", dx);
            VF_LOG_INFO("Preprocess - start");
         }

         grid->setDeltaX(dx);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

         grid->setPeriodicX1(false);
         grid->setPeriodicX2(false);
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

         SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         ppblocks->update(0);

         SPtr<D3Q27Interactor> cylinderInt(new D3Q27Interactor(cylinder, grid, noSlipBC, Interactor3D::INVERSESOLID));

         //double r = dynamicPointerCast<GbCylinder3D>(cylinder)->getRadius();
         //double cx1 = g_minX1;
         //double cx2 = cylinder->getX2Centroid();
         //double cx3 = cylinder->getX3Centroid();

         //mu::Parser fct;
         //fct.SetExpr("vx1");
         ////fct.SetExpr("vx1*(1-((x2-y0)^2+(x3-z0)^2)/(R^2))");
         ////fct.DefineConst("x2Vmax", 0.0); //x2-Pos fuer vmax
         ////fct.DefineConst("x3Vmax", 0.0); //x3-Pos fuer vmax
         ////fct.DefineConst("R", r);
         //fct.DefineConst("vx1", uLB);
         ////fct.DefineConst("x0", cx1);
         ////fct.DefineConst("y0", cx2);
         ////fct.DefineConst("z0", cx3);
         ////fct.DefineConst("nue", nuLB);
         //SPtr<BC> velBC(new VelocityBC(true, false, false, fct, 0, BCFunction::INFCONST));
         //velBC->setBCStrategy(SPtr<BCStrategy>(new VelocityInterpolated()));
         //velBC->setBCStrategy(SPtr<BCStrategy>(new VelocityWithPressureInterpolated()));
         
         inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBC1, Interactor3D::SOLID));
         //inflowInt->addBC(velBC2);


         //outflow
         SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, denBC, Interactor3D::SOLID));

         //SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::DIR_00M));
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(cylinderInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();

         ppblocks->update(0);
         ppblocks.reset();

         unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
         int ghostLayer = 3;
         unsigned long long numberOfNodesPerBlock = (unsigned long long)(blocknx[0])* (unsigned long long)(blocknx[1])* (unsigned long long)(blocknx[2]);
         unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
         unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
         real needMemAll = real(numberOfNodesPerBlockWithGhostLayer*(27 * sizeof(real) + sizeof(int) + sizeof(float) * 4));
         real needMem = needMemAll / real(comm->getNumberOfProcesses());

         if (myid == 0)
         {
            VF_LOG_INFO("Number of blocks = {}", numberOfBlocks);
            VF_LOG_INFO("Number of nodes  = {}", numberOfNodes);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for (int level = minInitLevel; level <= maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               VF_LOG_INFO("Number of blocks for level {} = {}", level, nobl);
               VF_LOG_INFO("Number of nodes for level {} = {}", level, nobl*numberOfNodesPerBlock);
            }
            VF_LOG_INFO("Necessary memory  = {} bytes", needMemAll);
            VF_LOG_INFO("Necessary memory per process = {} bytes", needMem);
            VF_LOG_INFO("Available memory per process = {} bytes", availMem);
         }

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         intHelper.setBC();

         //bcVisitor.addBC(velBC1);
         grid->accept(bcVisitor);

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor;
         //initVisitor.setVx1(fct);
         //initVisitor.setVx1(uLB);
         grid->accept(initVisitor);

         //boundary conditions grid
         {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            SPtr<SimulationObserver> ppgeo(new WriteBoundaryConditionsSimulationObserver(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
            ppgeo->update(0);
            ppgeo.reset();
         }

         if (myid == 0) VF_LOG_INFO("Preprocess - end");
      }
      else
      {
         if (myid == 0)
         {
            VF_LOG_INFO("Parameters:");
            VF_LOG_INFO("uLb = {}", uLB);
            VF_LOG_INFO("rho = {}", rhoLB);
            VF_LOG_INFO("nuLb = {}", nuLB);
            VF_LOG_INFO("Re = {}", Re);
            VF_LOG_INFO("dx = {}", dx);
            VF_LOG_INFO("number of levels = {}", refineLevel + 1);
            VF_LOG_INFO("numOfThreads = {}", numOfThreads);
            VF_LOG_INFO("path = {}", pathname);
         }

         migSimulationObserver->restart((int)restartStep);
         grid->setTimeStep(restartStep);

         if (myid == 0) VF_LOG_INFO("Restart - end");
      }

      OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
      grid->accept(setConnsVisitor);

      SPtr<Interpolator> iProcessor(new CompressibleOffsetMomentsInterpolator());
      SetInterpolationConnectorsBlockVisitor setInterConnsVisitor(comm, nuLB, iProcessor);
      grid->accept(setInterConnsVisitor);

      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      SPtr<SimulationObserver> pp(new WriteMacroscopicQuantitiesSimulationObserver(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      SPtr<UbScheduler> nupsSch(new UbScheduler(100, 100, 100000000));
      SPtr<SimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

      //SPtr<UbScheduler> timeBCSch(new UbScheduler(1, startTime, startTime));
      //auto timeDepBC = make_shared<TimeDependentBCSimulationObserver>(TimeDependentBCSimulationObserver(grid, timeBCSch));
      //timeDepBC->addInteractor(inflowInt);
#ifdef _OPENMP
      omp_set_num_threads(numOfThreads);
#endif
      numOfThreads = 1;
      SPtr<UbScheduler> stepGhostLayer(visSch);
      SPtr<Simulation> simulation(new Simulation(grid, stepGhostLayer, int(endTime)));
      simulation->addSimulationObserver(npr);
      simulation->addSimulationObserver(pp);
      simulation->addSimulationObserver(migSimulationObserver);
      //simulation->addSimulationObserver(timeDepBC);

      if (myid == 0) VF_LOG_INFO("Simulation-start");
      simulation->run();
      if (myid == 0) VF_LOG_INFO("Simulation-end");
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


int main(int argc, char *argv[])
{
    try {
        vf::logging::Logger::initializeLogger();

        VF_LOG_INFO("Starting VirtualFluids...");

        if (argc > 1)
            run(std::string(argv[1]));
        else
            VF_LOG_CRITICAL("Configuration file is missing!");

        VF_LOG_INFO("VirtualFluids is finished.");

    } catch (const spdlog::spdlog_ex &ex) {
        std::cout << "Log initialization failed: " << ex.what() << std::endl;
    }
}

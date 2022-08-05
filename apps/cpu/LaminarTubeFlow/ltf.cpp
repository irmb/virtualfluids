#include <iostream>
#include <string>
//#include <omp.h>

#include "VirtualFluids.h"

using namespace std;


void run(string configname)
{
   try
   {
      vf::basics::ConfigurationFile   config;
      config.load(configname);

      string          pathname = config.getValue<string>("pathname");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
      double          uLB = config.getValue<double>("uLB");
      double          endTime = config.getValue<double>("endTime");
      double          outTime = config.getValue<double>("outTime");
      double          availMem = config.getValue<double>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      double          Re = config.getValue<double>("Re");
      double          dx = config.getValue<double>("dx");
      vector<double>  length = config.getVector<double>("length");
      bool            logToFile = config.getValue<bool>("logToFile");
      double          restartStep = config.getValue<double>("restartStep");
      double          cpStart = config.getValue<double>("cpStart");
      double          cpStep = config.getValue<double>("cpStep");
      bool            newStart = config.getValue<bool>("newStart");

      SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
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

      LBMReal dLB = length[1] / dx;
      LBMReal rhoLB = 0.0;
      LBMReal nuLB = (uLB*dLB) / Re;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      //BC Adapter
      //////////////////////////////////////////////////////////////////////////////
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      //noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new ThinWallNoSlipBCAlgorithm()));
      //SPtr<BCAdapter> denBCAdapter(new DensityBCAdapter(rhoLB));
      //denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new EqDensityBCAlgorithm()));

      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapter(new DensityBCAdapter(rhoLB));
      denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));
      //denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

      //double startTime = 5;
      mu::Parser fct1;
      fct1.SetExpr("U");
      fct1.DefineConst("U", uLB);
      SPtr<BCAdapter> velBCAdapter1(new VelocityBCAdapter(true, false, false, fct1, 0, BCFunction::INFCONST));
      //velBCAdapter1->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));
      velBCAdapter1->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

      //mu::Parser fct2;
      //fct2.SetExpr("U");
      //fct2.DefineConst("U", uLB);
      //SPtr<BCAdapter> velBCAdapter2(new VelocityBCAdapter(true, false, false, fct2, startTime, BCFunction::INFCONST));

      //////////////////////////////////////////////////////////////////////////////////
      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(denBCAdapter);
      //bcVisitor.addBC(velBCAdapter1);

      SPtr<Grid3D> grid(new Grid3D(comm));

      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new BCProcessor());

      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      //double bulckViscosity = 3700*nuLB;
      //dynamicPointerCast<CompressibleCumulant4thOrderViscosityLBMKernel>(kernel)->setBulkViscosity(bulckViscosity);
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CumulantK17LBMKernel());
      kernel->setBCProcessor(bcProc);
      kernel->setBCProcessor(bcProc);

      //////////////////////////////////////////////////////////////////////////
      SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::DIR_00M));
      //restart
      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      //SPtr<MPIIOMigrationCoProcessor> migCoProcessor(new MPIIOMigrationCoProcessor(grid, mSch, metisVisitor, pathname + "/mig", comm));
      SPtr<MPIIOMigrationBECoProcessor> migCoProcessor(new MPIIOMigrationBECoProcessor(grid, mSch, metisVisitor, pathname + "/mig", comm));
      migCoProcessor->setLBMKernel(kernel);
      migCoProcessor->setBCProcessor(bcProc);
      migCoProcessor->setNu(nuLB);
      migCoProcessor->setNuLG(0.01, 0.01);
      migCoProcessor->setDensityRatio(1);
      //////////////////////////////////////////////////////////////////////////

      SPtr<D3Q27Interactor> inflowInt;

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

         SPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         ppblocks->process(0);

         SPtr<D3Q27Interactor> cylinderInt(new D3Q27Interactor(cylinder, grid, noSlipBCAdapter, Interactor3D::INVERSESOLID));

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
         //SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
         //velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));
         //velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));
         
         inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter1, Interactor3D::SOLID));
         //inflowInt->addBCAdapter(velBCAdapter2);


         //outflow
         SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));

         //SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::DIR_00M));
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

         bcVisitor.addBC(velBCAdapter1);
         grid->accept(bcVisitor);

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor;
         //initVisitor.setVx1(fct);
         //initVisitor.setVx1(uLB);
         grid->accept(initVisitor);

         //boundary conditions grid
         {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            SPtr<CoProcessor> ppgeo(new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
            ppgeo->process(0);
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

         migCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);

         if (myid == 0) VF_LOG_INFO("Restart - end");
      }

      OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
      grid->accept(setConnsVisitor);

      SPtr<InterpolationProcessor> iProcessor(new CompressibleOffsetMomentsInterpolationProcessor());
      SetInterpolationConnectorsBlockVisitor setInterConnsVisitor(comm, nuLB, iProcessor);
      grid->accept(setInterConnsVisitor);

      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      SPtr<CoProcessor> pp(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      SPtr<UbScheduler> nupsSch(new UbScheduler(100, 100, 100000000));
      SPtr<CoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      //SPtr<UbScheduler> timeBCSch(new UbScheduler(1, startTime, startTime));
      //auto timeDepBC = make_shared<TimeDependentBCCoProcessor>(TimeDependentBCCoProcessor(grid, timeBCSch));
      //timeDepBC->addInteractor(inflowInt);

      omp_set_num_threads(numOfThreads);
      numOfThreads = 1;
      SPtr<UbScheduler> stepGhostLayer(visSch);
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, int(endTime)));
      calculator->addCoProcessor(npr);
      calculator->addCoProcessor(pp);
      calculator->addCoProcessor(migCoProcessor);
      //calculator->addCoProcessor(timeDepBC);

      if (myid == 0) VF_LOG_INFO("Simulation-start");
      calculator->calculate();
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
        vf::logging::Logger::initalizeLogger();

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

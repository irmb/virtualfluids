#include <VirtualFluids.h>
#include <set>
#include <map>
using namespace std;


////////////////////////////////////////////////////////////////////////
void run(string configname)
{
    using namespace vf::lbm::dir;

   try
   {
      SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();

      int myid = comm->getProcessID();

      //vf::basics::ConfigurationFile config;
      //config.load(configname);

      //string pathname = config.getValue<string>("path");
      //double availMem = config.getValue<double>("memory");
      //string metafile = config.getValue<string>("metafile");
      //double outstep  = config.getValue<double>("outstep");
      //double endstep        = config.getValue<double>("endstep");
      //int numOfThreads      = config.getValue<int>("threads");
      //const int refineLevel = config.getValue<int>("level");

      string outputPath = "d:/temp/sphereBlock_5_SBB";
      real availMem = 8e9;
      real outstep = 10000;
      real endstep = 1e6;
      int numOfThreads = 4;
      omp_set_num_threads(numOfThreads);
      int refineLevel = 0;

      real radius = 5;
      real uLB = 1e-3;
      real Re = 1;
      real rhoLB = 0.0;
      real nuLB = (uLB*2.0*radius)/Re;

      real dp_LB = 1e-6;
//      double rhoLBinflow = dp_LB*3.0;

      SPtr<BC> noSlipBC(new NoSlipBC());
      noSlipBC->setBCStrategy(SPtr<BCStrategy>(new NoSlipBCStrategy()));
      SPtr<BC> slipBC(new SlipBC());
      slipBC->setBCStrategy(SPtr<BCStrategy>(new SimpleSlipBCStrategy()));
      
      real H = 50;
      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);
      //mu::Parser fct;
      //fct.SetExpr("16*U*x2*x3*(H-x2)*(H-x3)/H^4");
      //fct.DefineConst("U", uLB);
      //fct.DefineConst("H", H);
      SPtr<BC> velBC(new VelocityBC(true, false, false, fct, 0, BCFunction::INFCONST));
      //velBC->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));
      velBC->setBCStrategy(SPtr<BCStrategy>(new SimpleVelocityBCStrategy()));

      SPtr<BC> denBC(new DensityBC(rhoLB));
      denBC->setBCStrategy(SPtr<BCStrategy>(new NonEqDensityBCStrategy()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBC);
      bcVisitor.addBC(slipBC);
      bcVisitor.addBC(velBC);
      bcVisitor.addBC(denBC);

      real dx = 1;

      const int blocknx1 = 50;
      const int blocknx2 = 50;
      const int blocknx3 = 50;

      const int gridNx1 = 150;
      const int gridNx2 = H;
      const int gridNx3 = H;

      real L1, L2, L3;
      L1 = gridNx1;
      L2 = gridNx2;
      L3 = gridNx3;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);

      //sphere
      //SPtr<GbObject3D> sphere(new GbSphere3D(L1 * 0.5, L2 * 0.5, L3 * 0.5, radius));
      SPtr<GbObject3D> sphere(new GbSphere3D(75, 25, 25, radius));
      GbSystem3D::writeGeoObject(sphere.get(), outputPath + "/geo/sphere", WbWriterVtkXmlBinary::getInstance());
      SPtr<D3Q27Interactor> sphereInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(sphere, grid, noSlipBC, Interactor3D::SOLID));

      if (true)
      {
         //bounding box
         real d_minX1 = 0.0;
         real d_minX2 = 0.0;
         real d_minX3 = 0.0;

         real d_maxX1 = L1;
         real d_maxX2 = L2;
         real d_maxX3 = L3;

         real blockLength = blocknx1*dx;

         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "uLB = " << uLB);
            UBLOG(logINFO, "rhoLB = " << rhoLB);
            UBLOG(logINFO, "nueLB = " << nuLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << dx);
            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "Preprozess - start");
         }

         SPtr<GbObject3D> gridCube(new GbCuboid3D(d_minX1, d_minX2, d_minX3, d_maxX1, d_maxX2, d_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), outputPath + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         real off = 0.0;
         SPtr<GbObject3D> refCube(new GbCuboid3D(sphere->getX1Minimum() - off, sphere->getX2Minimum() - off, sphere->getX3Minimum(),
            sphere->getX1Maximum() + off, sphere->getX2Maximum() + off, sphere->getX3Maximum()));
         if (myid == 0) GbSystem3D::writeGeoObject(refCube.get(), outputPath + "/geo/refCube", WbWriterVtkXmlBinary::getInstance());

         if (refineLevel > 0)
         {
            if (myid == 0) UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
            //refineHelper.addGbObject(sphere, refineLevel);
            refineHelper.addGbObject(refCube, refineLevel);
            refineHelper.refine();
            if (myid == 0) UBLOG(logINFO, "Refinement - end");
         }

         //walls
         GbCuboid3DPtr addWallYmin(new GbCuboid3D(d_minX1 - 4.0*blockLength, d_minX2 - 4.0*blockLength, d_minX3 - 4.0*blockLength, d_maxX1 + 4.0*blockLength, d_minX2, d_maxX3 + 4.0*blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), outputPath + "/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin(new GbCuboid3D(d_minX1 - 4.0*blockLength, d_minX2 - 4.0*blockLength, d_minX3 - 4.0*blockLength, d_maxX1 + 4.0*blockLength, d_maxX2 + 4.0*blockLength, d_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), outputPath + "/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax(new GbCuboid3D(d_minX1 - 4.0*blockLength, d_maxX2, d_minX3 - 4.0*blockLength, d_maxX1 + 4.0*blockLength, d_maxX2 + 4.0*blockLength, d_maxX3 + 4.0*blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), outputPath + "/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(d_minX1 - 4.0*blockLength, d_minX2 - 4.0*blockLength, d_maxX3, d_maxX1 + 4.0*blockLength, d_maxX2 + 4.0*blockLength, d_maxX3 + 4.0*blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), outputPath + "/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(d_minX1 - 4.0*blockLength, d_minX2 - 4.0*blockLength, d_minX3 - 4.0*blockLength, d_minX1, d_maxX2 + 4.0*blockLength, d_maxX3 + 4.0*blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), outputPath + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(d_maxX1, d_minX2 - 4.0*blockLength, d_minX3 - 4.0*blockLength, d_maxX1 + 4.0*blockLength, d_maxX2 + 4.0*blockLength, d_maxX3 + 4.0*blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), outputPath + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), outputPath, WbWriterVtkXmlBinary::getInstance(), comm));

         //walls
         SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, slipBC, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, slipBC, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, slipBC, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, slipBC, Interactor3D::SOLID));

         //inflow
         SPtr<D3Q27Interactor> inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBC, Interactor3D::SOLID));

         //D3Q27BoundaryConditionAdapterPtr denBCInflow(new D3Q27DensityBC(rhoLBinflow));
         //denBCInflow->setSecondaryBcOption(0);
         //SPtr<D3Q27Interactor> inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, denBCInflow, Interactor3D::SOLID));

         //outflow
         SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, denBC, Interactor3D::SOLID));

         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, DIR_00M));
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(sphereInt);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();

         ////domain decomposition for threads
         //PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         //grid->accept(pqPartVisitor);

         ppblocks->update(0);
         ppblocks.reset();

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nod = nob * (blocknx1 + gl) * (blocknx2 + gl) * (blocknx3 + gl);

         real needMemAll = real(nod*(27 * sizeof(real) + sizeof(int) + sizeof(float) * 4));
         real needMem = needMemAll / real(comm->getNumberOfProcesses());

         if (myid == 0)
         {
            UBLOG(logINFO, "Number of blocks = " << nob);
            UBLOG(logINFO, "Number of nodes  = " << nod);
            UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
            UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
            UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
         }

         SPtr<LBMKernel> kernel(new IncompressibleCumulantLBMKernel());
         //SPtr<LBMKernel> kernel(new CompressibleCumulantLBMKernel());

         SPtr<BCSet> bcSet(new BCSet());


         kernel->setBCSet(bcSet);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
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
         fctRoh.DefineConst("l", d_maxX1 - d_minX1);

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor;
         initVisitor.setVx1(fct);
         //initVisitor.setRho(fctRoh);
         grid->accept(initVisitor);

         //Postrozess
         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         SPtr<SimulationObserver> ppgeo(
            new WriteBoundaryConditionsSimulationObserver(grid, geoSch, outputPath, WbWriterVtkXmlBinary::getInstance(), comm));
         ppgeo->update(0);
         ppgeo.reset();

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         UBLOG(logINFO, "SetConnectors - start, id=" << myid);
      }

      

      UBLOG(logINFO, "SetConnectors - start, id=" << myid);
      //set connectors
      //SPtr<InterpolationProcessor> iProcessor(new  IncompressibleOffsetInterpolationProcessor());
      //SPtr<CompressibleOffsetMomentsInterpolationProcessor> iProcessor(new  CompressibleOffsetMomentsInterpolationProcessor());
      //SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);

      OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
      grid->accept(setConnsVisitor);

      SPtr<InterpolationProcessor> iProcessor(new CompressibleOffsetMomentsInterpolationProcessor());
      SetInterpolationConnectorsBlockVisitor setInterConnsVisitor(comm, nuLB, iProcessor);
      grid->accept(setInterConnsVisitor);

      UBLOG(logINFO, "SetConnectors - end, id=" << myid);

      SPtr<UbScheduler> stepSch(new UbScheduler(outstep));
      //stepSch->addSchedule(10000, 0, 1000000);
      SPtr<WriteMacroscopicQuantitiesSimulationObserver> writeMQSimulationObserver(new WriteMacroscopicQuantitiesSimulationObserver(grid, stepSch, outputPath, WbWriterVtkXmlBinary::getInstance(), conv,comm));

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      SPtr<SimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

      real area = UbMath::PI * radius * radius;
      SPtr<UbScheduler> forceSch(new UbScheduler(100));
      SPtr<CalculateForcesSimulationObserver> fp = make_shared<CalculateForcesSimulationObserver>(grid, forceSch, outputPath + "/forces/forces.txt", comm, uLB, area);
      fp->addInteractor(sphereInt);

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Simulation> simulation(new Simulation(grid, stepGhostLayer, endstep));
      simulation->addSimulationObserver(npr);
      simulation->addSimulationObserver(fp);
      simulation->addSimulationObserver(writeMQSimulationObserver);


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
      //if (argv[1] != NULL)
      //{
      //   run(string(argv[1]));
      //}
      //else
      //{
      //   cout << "Configuration file is missing!" << endl;
      //}
      run("sphere");
   }
}


#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;


void run()
{
    using namespace vf::lbm::dir;

   try
   {
      SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      int    numOfThreads = 4;
      real availMem = 5e9;

      

      //////////////////////////////////////////////////////////////////////////
      //DLR-F16 test
      //dx_coarse = 0.003 mm

      string  pathname = "d:/temp/ConvectionOfVortex_0.003_4th";
      int     endTime = 10000;
      real  outTime = 10;
      real dx =  0.003;
      real rhoLB = 0.0;
      real nuLB = 8.66025e-6;
      real yFactor = 1.0;

      //string  pathname = "d:/temp/ConvectionOfVortex_0.003_square";
      //int     endTime = 20;
      //double  outTime = 10;
      //LBMReal dx =  0.003;
      //LBMReal rhoLB = 0.0;
      //LBMReal nuLB = 8.66025e-6;

      //////////////////////////////////////////////////////////////////////////
      ////dx_coarse = 0.0015 mm
      //string  pathname = "d:/temp/ConvectionOfVortex_0.0015";
      //double  endTime = 40;
      //double  outTime = 40;
      //LBMReal dx =  0.0015;
      //LBMReal rhoLB = 0.0;
      //LBMReal nuLB = 8.66025e-6*2.0;
      ////////////////////////////////////////////////////////////////////////////
      //dx_coarse = 0.00075 mm

      //string  pathname = "d:/temp/ConvectionOfVortex_0.00075_4th_moments";
      //double  endTime = 2000;
      //double  outTime = 10;
      //LBMReal dx =  0.00075;
      //LBMReal rhoLB = 0.0;
      //LBMReal nuLB = 8.66025e-6*4.0;
      //double yFactor = 4.0;

      //string  pathname = "d:/temp/ConvectionOfVortex_0.00075_moments";
      //double  endTime = 160;
      //double  outTime = 160;
      //LBMReal dx =  0.00075;
      //LBMReal rhoLB = 0.0;
      //LBMReal nuLB = 8.66025e-6*4.0;

      //////////////////////////////////////////////////////////////////////////
      ////dx_coarse = 0.000375 mm
      //string  pathname = "d:/temp/ConvectionOfVortex_0.000375";
      //double  endTime = 80;
      //double  outTime = 80;
      //LBMReal dx =  0.00075;
      //LBMReal rhoLB = 0.0;
      //LBMReal nuLB = 8.66025e-6*8.0;
      //////////////////////////////////////////////////////////////////////////

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      //int baseLevel = 0;
      int refineLevel = 1;

      //bounding box
      real g_minX1 = -0.045;
      real g_minX2 = -0.015/yFactor;
      real g_minX3 = -0.06;

      real g_maxX1 = 0.045;
      real g_maxX2 = 0.015/yFactor;
      real g_maxX3 = 0.06;

      vector<int>  blocknx(3);
      blocknx[0] = 10;
      blocknx[1] = 10;
      blocknx[2] = 10;

      //geometry
      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


      real blockLength = blocknx[0] * dx;

      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(true);


      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      SPtr<BC> outflowBC(new DensityBC(rhoLB));
      outflowBC->setBCStrategy(SPtr<BCStrategy>(new NonReflectingOutflowBCStrategy()));

      BoundaryConditionsBlockVisitor bcVisitor;

      SPtr<BCSet> bcProc;
      bcProc = SPtr<BCSet>(new BCSet());

      SPtr<GbObject3D> refCube(new GbCuboid3D(g_minX1-blockLength,-0.02,-0.02,g_maxX1+blockLength,0.02,0.02));
      if (myid==0) GbSystem3D::writeGeoObject(refCube.get(), pathname+"/geo/refCube", WbWriterVtkXmlBinary::getInstance());

      if (refineLevel>0)
      {
         if (myid==0) UBLOG(logINFO, "Refinement - start");
         RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
         refineHelper.addGbObject(refCube, refineLevel);
         refineHelper.refine();
         if (myid==0) UBLOG(logINFO, "Refinement - end");
      }

      SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

      //outflow
      GbCuboid3DPtr geoOutflow1(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
      if (myid==0) GbSystem3D::writeGeoObject(geoOutflow1.get(), pathname+"/geo/geoOutflow1", WbWriterVtkXmlASCII::getInstance());
      SPtr<D3Q27Interactor> outflowIntr1 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow1, grid, outflowBC, Interactor3D::SOLID));

      GbCuboid3DPtr geoOutflow2(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      if (myid==0) GbSystem3D::writeGeoObject(geoOutflow2.get(), pathname+"/geo/geoOutflow2", WbWriterVtkXmlASCII::getInstance());
      SPtr<D3Q27Interactor> outflowIntr2 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow2, grid, outflowBC, Interactor3D::SOLID));
      
      GbCuboid3DPtr geoOutflow3(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
      if (myid==0) GbSystem3D::writeGeoObject(geoOutflow3.get(), pathname+"/geo/geoOutflow3", WbWriterVtkXmlASCII::getInstance());
      SPtr<D3Q27Interactor> outflowIntr3 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow3, grid, outflowBC, Interactor3D::SOLID));

      GbCuboid3DPtr geoOutflow4(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      if (myid==0) GbSystem3D::writeGeoObject(geoOutflow4.get(), pathname+"/geo/geoOutflow4", WbWriterVtkXmlASCII::getInstance());
      SPtr<D3Q27Interactor> outflowIntr4 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow4, grid, outflowBC, Interactor3D::SOLID));

      SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, d00M));
      InteractorsHelper intHelper(grid, metisVisitor);
      //intHelper.addInteractor(outflowIntr1);
      //intHelper.addInteractor(outflowIntr2);
      //intHelper.addInteractor(outflowIntr3);
      //intHelper.addInteractor(outflowIntr4);
      intHelper.selectBlocks();

      ppblocks->update(0);
      ppblocks.reset();

      //set connectors  
      //SPtr<Interpolator> iProcessor(new CompressibleOffsetInterpolator());
      //SPtr<Interpolator> iProcessor(new CompressibleOffsetMomentsInterpolator());
      //dynamicPointerCast<CompressibleOffsetMomentsInterpolator>(iProcessor)->setBulkOmegaToOmega(true);
      //SPtr<Interpolator> iProcessor(new CompressibleOffsetSquarePressureInterpolator());

      OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
      grid->accept(setConnsVisitor);

      SPtr<Interpolator> iProcessor(new CompressibleOffsetMomentsInterpolator());
      SetInterpolationConnectorsBlockVisitor setInterConnsVisitor(comm, nuLB, iProcessor);
      grid->accept(setInterConnsVisitor);

      UBLOG(logINFO, "SetConnectorsBlockVisitor:start");
      grid->accept(setConnsVisitor);
      UBLOG(logINFO, "SetConnectorsBlockVisitor:end");

      unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
      int ghostLayer = 3;
      unsigned long long numberOfNodesPerBlock = (unsigned long long)(blocknx[0])* (unsigned long long)(blocknx[1])* (unsigned long long)(blocknx[2]);
      unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
      unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
      real needMemAll = real(numberOfNodesPerBlockWithGhostLayer*(27 * sizeof(real) + sizeof(int) + sizeof(float) * 4));
      real needMem = needMemAll / real(comm->getNumberOfProcesses());

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


      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      //dynamicPointerCast<CompressibleCumulant4thOrderViscosityLBMKernel>(kernel)->setBulkViscosity(10.0*nuLB);
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
      //dynamicPointerCast<CompressibleCumulantLBMKernel>(kernel)->setBulkOmegaToOmega(true);
      //
      SPtr<BCSet> bcSet(new BCSet());

      kernel->setBCSet(bcSet);

      SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
      grid->accept(kernelVisitor);

      if (refineLevel>0)
      {
         SetUndefinedNodesBlockVisitor undefNodesVisitor;
         grid->accept(undefNodesVisitor);
      }

      intHelper.setBC();

      real Ma = 0.005;

      mu::Parser initRho, initVx1, initVx2; 
      initRho.SetExpr("rhoLB + (-(rho0*epsilon^2)/2) * exp(1-(scaleFactor*(x1^2+x3^2))/R^2) + (1/(2*gamma*rho0)) * ((-(rho0*epsilon^2)/2) * exp(1-(scaleFactor*(x1^2+x3^2))/R^2))^2");
      initRho.DefineConst("rhoLB", rhoLB);
      initRho.DefineConst("rho0", rhoLB+1.0);
      initRho.DefineConst("R", 0.1);
      initRho.DefineConst("gamma", 0.1);
      initRho.DefineConst("epsilon", 0.14);
      initRho.DefineConst("scaleFactor", 277.777777779);

      initVx1.SetExpr("-epsilon*c0*((x3*scaleFactor1)/R)*exp(0.5*(1-scaleFactor*(x1^2+x3^2)/R^2))");
      initVx1.DefineConst("c0", 1.0/std::sqrt(3.0));
      initVx1.DefineConst("scaleFactor", 277.777777779);
      initVx1.DefineConst("scaleFactor1", 16.6666666667);
      initVx1.DefineConst("epsilon", 0.14);
      initVx1.DefineConst("R", 0.1);

      initVx2.SetExpr("V0 + epsilon*c0*((x1*scaleFactor1)/R)*exp(0.5*(1-scaleFactor*(x1^2+x3^2)/R^2))");
      initVx2.DefineConst("c0", 1.0/std::sqrt(3.0));
      initVx2.DefineConst("scaleFactor", 277.777777779);
      initVx2.DefineConst("scaleFactor1", 16.6666666667);
      initVx2.DefineConst("epsilon", 0.14);
      initVx2.DefineConst("R", 0.1);
      initVx2.DefineConst("V0", -Ma/(1.0/std::sqrt(3.0)));


      //initialization of distributions
      InitDistributionsBlockVisitor initVisitor;
      initVisitor.setRho(initRho);
      initVisitor.setVx1(initVx1);
      initVisitor.setVx3(initVx2);
      grid->accept(initVisitor);

      grid->accept(bcVisitor);

      //Postrozess
      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      SPtr<SimulationObserver> ppgeo(new WriteBoundaryConditionsSimulationObserver(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
      ppgeo->update(0);
      ppgeo.reset();

      if (myid==0) UBLOG(logINFO, "Preprozess - end");

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      SPtr<WriteMacroscopicQuantitiesSimulationObserver> writeMQSimulationObserver(new WriteMacroscopicQuantitiesSimulationObserver(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
      writeMQSimulationObserver->update(0);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      std::shared_ptr<NUPSCounterSimulationObserver> nupsSimulationObserver(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

      //SPtr<UbScheduler> tavSch(new UbScheduler(1, 0, endTime));
      //SPtr<TimeAveragedValuesSimulationObserver> tav(new TimeAveragedValuesSimulationObserver(grid, pathname, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,
      //   TimeAveragedValuesSimulationObserver::Density | TimeAveragedValuesSimulationObserver::Velocity | TimeAveragedValuesSimulationObserver::Fluctuations));
      //tav->setWithGhostLayer(true);

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Simulation> simulation(new Simulation(grid, stepGhostLayer, endTime));
      simulation->addSimulationObserver(nupsSimulationObserver);
      simulation->addSimulationObserver(writeMQSimulationObserver);
      //simulation->addSimulationObserver(tav);

      //omp_set_num_threads(1);

      if (myid==0) UBLOG(logINFO, "Simulation-start");
      simulation->run();
      if (myid==0) UBLOG(logINFO, "Simulation-end");
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
   run();
}


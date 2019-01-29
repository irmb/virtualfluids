#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;


void run()
{
   try
   {
      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      int    numOfThreads = 4;
      double availMem = 5e9;

      

      //////////////////////////////////////////////////////////////////////////
      //DLR-F16 test
      //dx_coarse = 0.003 mm
      //string  pathname = "d:/temp/ConvectionOfVortex_0.003";
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
      string  pathname = "d:/temp/ConvectionOfVortex_0.00075_moments";
      double  endTime = 80;
      double  outTime = 80;
      LBMReal dx =  0.00075;
      LBMReal rhoLB = 0.0;
      LBMReal nuLB = 8.66025e-6*4.0;
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

      int baseLevel = 0;
      int refineLevel = 1;

      //bounding box
      double g_minX1 = -0.045;
      double g_minX2 = -0.015;
      double g_minX3 = -0.06;

      double g_maxX1 = 0.045;
      double g_maxX2 = 0.015;
      double g_maxX3 = 0.06;

      vector<int>  blocknx(3);
      blocknx[0] = 10;
      blocknx[1] = 10;
      blocknx[2] = 10;

      //geometry
      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


      double blockLength = blocknx[0] * dx;

      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(true);


      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

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

      SPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

      SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
      InteractorsHelper intHelper(grid, metisVisitor);
      intHelper.selectBlocks();

      ppblocks->process(0);
      ppblocks.reset();

      //set connectors  
      //SPtr<InterpolationProcessor> iProcessor(new CompressibleOffsetInterpolationProcessor());
      SPtr<InterpolationProcessor> iProcessor(new CompressibleOffsetMomentsInterpolationProcessor());
      //dynamicPointerCast<CompressibleOffsetMomentsInterpolationProcessor>(iProcessor)->setBulkOmegaToOmega(true);
      SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);

      UBLOG(logINFO, "SetConnectorsBlockVisitor:start");
      grid->accept(setConnsVisitor);
      UBLOG(logINFO, "SetConnectorsBlockVisitor:end");

      //domain decomposition for threads
      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      grid->accept(pqPartVisitor);


      unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
      int ghostLayer = 3;
      unsigned long long numberOfNodesPerBlock = (unsigned long long)(blocknx[0])* (unsigned long long)(blocknx[1])* (unsigned long long)(blocknx[2]);
      unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
      unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
      double needMemAll = double(numberOfNodesPerBlockWithGhostLayer*(27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
      double needMem = needMemAll / double(comm->getNumberOfProcesses());

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
      SPtr<BCProcessor> bcProcessor(new BCProcessor());

      kernel->setBCProcessor(bcProcessor);

      SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
      grid->accept(kernelVisitor);

      if (refineLevel>0)
      {
         SetUndefinedNodesBlockVisitor undefNodesVisitor;
         grid->accept(undefNodesVisitor);
      }

      double Ma = 0.15;

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

      //Postrozess
      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      SPtr<CoProcessor> ppgeo(new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
      ppgeo->process(0);
      ppgeo.reset();

      if (myid==0) UBLOG(logINFO, "Preprozess - end");

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      SPtr<UbScheduler> visSch(new UbScheduler(outTime));
      SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
      writeMQCoProcessor->process(0);

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      std::shared_ptr<NUPSCounterCoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> tavSch(new UbScheduler(1, 0, endTime));
      SPtr<TimeAveragedValuesCoProcessor> tav(new TimeAveragedValuesCoProcessor(grid, pathname, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,
         TimeAveragedValuesCoProcessor::Density | TimeAveragedValuesCoProcessor::Velocity | TimeAveragedValuesCoProcessor::Fluctuations));
      tav->setWithGhostLayer(true);

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(writeMQCoProcessor);
      calculator->addCoProcessor(tav);

      //omp_set_num_threads(1);

      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
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


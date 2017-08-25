#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;


void run(string configname)
{
   try
   {
      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      int    numOfThreads = 4;
      double availMem = 5e9;

      //40
      //string  pathname = "d:/temp/AcousticPulse40Cube2y";
      //double  endTime = 20;
      //double  outTime = 20;
      //LBMReal dx =  0.05;

      //80
      //string  pathname = "d:/temp/AcousticPulse80Cube2y";
      //double  endTime = 40;
      //double  outTime = 40;
      //LBMReal dx = 0.025;

      //160
      string  pathname = "d:/temp/AcousticPulse160Cube2y";
      double  endTime = 80;
      double  outTime = 80;
      LBMReal dx = 0.0125;

      //LBMReal dx = 0.1; 
      //LBMReal dx = 1.66666666667e-2; //120
      
      LBMReal rhoLB = 0.0;
      LBMReal nuLB = 3.97e-7;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      int baseLevel = 0;
      int refineLevel = 1;

      //bounding box
      double g_minX1 = -1.0;
      double g_minX2 = -1.0;
      double g_minX3 = -1.0;

      double g_maxX1 = 1.0;
      double g_maxX2 = 1.0;
      double g_maxX3 = 1.0;

      //double g_minX1 = 0.0;
      //double g_minX2 = 0.0;
      //double g_minX3 = 0.0;

      //double g_maxX1 = 5.0;
      //double g_maxX2 = 5.0;
      //double g_maxX3 = dx;

      vector<int>  blocknx(3);
      blocknx[0] = 10;
      blocknx[1] = 10;
      blocknx[2] = 10;

      //geometry
      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


      double blockLength = blocknx[0] * dx;

      Grid3DPtr grid(new Grid3D(comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      Grid3DPtr oldGrid(new Grid3D(comm));
      UbSchedulerPtr rSch(new UbScheduler(10));
      //MPIIORestartCoProcessor rcp(oldGrid, rSch, pathname, comm);
      //rcp.restart(0);
      //////////////////////////////////////////////////////////////////////////

      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(true);


      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      GbObject3DPtr refCube(new GbCuboid3D(-0.4,-0.4,-0.4,0.4,0.4,0.4));
      if (myid==0) GbSystem3D::writeGeoObject(refCube.get(), pathname+"/geo/refCube", WbWriterVtkXmlBinary::getInstance());

      if (refineLevel>0)
      {
         if (myid==0) UBLOG(logINFO, "Refinement - start");
         RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
         refineHelper.addGbObject(refCube, refineLevel);
         refineHelper.refine();
         if (myid==0) UBLOG(logINFO, "Refinement - end");
      }

      WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

      Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
      InteractorsHelper intHelper(grid, metisVisitor);
      intHelper.selectBlocks();

      ppblocks->process(0);
      ppblocks.reset();

      //set connectors
      //InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
      InterpolationProcessorPtr iProcessor(new CompressibleOffsetMomentsInterpolationProcessor());
      boost::dynamic_pointer_cast<CompressibleOffsetMomentsInterpolationProcessor>(iProcessor)->setBulkOmegaToOmega(true);
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


      LBMKernelPtr kernel = LBMKernelPtr(new CompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], CompressibleCumulantLBMKernel::NORMAL));
      boost::dynamic_pointer_cast<CompressibleCumulantLBMKernel>(kernel)->setBulkOmegaToOmega(true);
      //
      BCProcessorPtr bcProcessor(new BCProcessor());

      kernel->setBCProcessor(bcProcessor);

      SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
      grid->accept(kernelVisitor);

      if (refineLevel>0)
      {
         SetUndefinedNodesBlockVisitor undefNodesVisitor;
         grid->accept(undefNodesVisitor);
      }

      mu::Parser fctRoh;
      //z
      //fctRoh.SetExpr("epsilon*exp(-alpha*(x1*x1+x2*x2))");
      //x
      //fctRoh.SetExpr("epsilon*exp(-alpha*(x3*x3+x2*x2))");
      //y
      fctRoh.SetExpr("epsilon*exp(-alpha*(x3*x3+x1*x1))");

      fctRoh.DefineConst("epsilon", 1e-3);
      fctRoh.DefineConst("alpha", log(2.0)/(0.01));
      //fctRoh.SetExpr("x1*0.001");

      //initialization of distributions
      InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
      initVisitor.setRho(fctRoh);
      grid->accept(initVisitor);

      //InitDistributionsWithCoarseGridBlockVisitor initVisitor(oldGrid, grid, iProcessor, nuLB);
      //grid->accept(initVisitor);


      //Postrozess
      UbSchedulerPtr geoSch(new UbScheduler(1));
      WriteBoundaryConditionsCoProcessorPtr ppgeo(
         new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
      ppgeo->process(0);
      ppgeo.reset();

      if (myid==0) UBLOG(logINFO, "Preprozess - end");

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      UbSchedulerPtr visSch(new UbScheduler(outTime));
      WriteMacroscopicQuantitiesCoProcessor pp(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
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
int main(int argc, char* argv[])
{
   //if (argv != NULL)
   //{
      //if (argv[1] != NULL)
      //{
      //   run(string(argv[1]));
      //}
      //else
      //{
      //   cout << "Configuration file is missing!" << endl;
      //}
   //}

   run("hallo");
}


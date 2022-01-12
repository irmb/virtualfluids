#include <iostream>
#include <string>
#include <omp.h>

#include "VirtualFluids.h"

using namespace std;


void run()
{
   try
   {
      SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      int    numOfThreads = 4;
      double availMem = 5e9;

      //40
      string  pathname = "d:/temp/InterfaceTest";
      int  endTime = 2000;
      double  outTime = 100;
      LBMReal dx =  0.05;
      
      LBMReal rhoLB = 0.0;
      LBMReal nuLB = 3.97e-7;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      int baseLevel = 0;
      int refineLevel = 1;

      //bounding box
      double g_minX1 = -1.0;
      double g_minX2 = -1.0;
      double g_minX3 = -1.0;

      double g_maxX1 = 1.0;
      double g_maxX2 = 1.0;
      double g_maxX3 = 1.0;

      vector<int>  blocknx(3);
      blocknx[0] = 10;
      blocknx[1] = 10;
      blocknx[2] = 10;

      //geometry
      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


      double blockLength = blocknx[0] * dx;

      SPtr<Grid3D> grid(new Grid3D(comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<Grid3D> oldGrid(new Grid3D(comm));
      SPtr<UbScheduler> rSch(new UbScheduler(10));
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

      SPtr<GbObject3D> refCube(new GbCuboid3D(-0.4,-0.4,-0.4,0.4,0.4,0.4));
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
      //InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
      InterpolationProcessorPtr iProcessor(new CompressibleOffsetMomentsInterpolationProcessor());
      dynamicPointerCast<CompressibleOffsetMomentsInterpolationProcessor>(iProcessor)->setBulkOmegaToOmega(true);
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


      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
      dynamicPointerCast<CompressibleCumulantLBMKernel>(kernel)->setRelaxationParameter(CompressibleCumulantLBMKernel::NORMAL);
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

      //initialization of distributions
      InitDistributionsBlockVisitor initVisitor;
      initVisitor.setVx1(0.001);
      initVisitor.setVx2(0.001);
      initVisitor.setVx3(0.001);
      grid->accept(initVisitor);

      //InitDistributionsWithCoarseGridBlockVisitor initVisitor(oldGrid, grid, iProcessor, nuLB);
      //grid->accept(initVisitor);


      //Postrozess
      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      SPtr<CoProcessor> ppgeo(new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
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

      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      SPtr<NUPSCounterCoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      omp_set_num_threads(numOfThreads);
      SPtr<Calculator> calculator(new BasicCalculator(grid, visSch, endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(writeMQCoProcessor);

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

   run();
}


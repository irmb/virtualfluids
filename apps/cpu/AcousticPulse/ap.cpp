#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;


void run()
{
   try
   {
      SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      int    numOfThreads = 4;
      double availMem = 5e9;

      //40
      //string  pathname = "d:/temp/AcousticPulse40Cube2y_test";
      //double  endTime = 20;
      //double  outTime = 20;
      //LBMReal dx =  0.05;

      //80
      //string  pathname = "d:/temp/AcousticPulse80Cube2y";
      //double  endTime = 40;
      //double  outTime = 40;
      //LBMReal dx = 0.025;

      //160
      //string  pathname = "d:/temp/AcousticPulse160Cube2y";
      //double  endTime = 80;
      //double  outTime = 80;
      //LBMReal dx = 0.0125;

      //LBMReal dx = 0.1; 
      //LBMReal dx = 1.66666666667e-2; //120
      
      //LBMReal rhoLB = 0.0;
      //LBMReal nuLB = 3.97e-7;

      //////////////////////////////////////////////////////////////////////////
      //DLR-F16 test
      ////dx_coarse = 0.003 mm
      string  pathname = "d:/temp/AcousticPulseXZ-4th-0.003";
      int     endTime = 20;
      double  outTime = 20;
      LBMReal dx =  0.003;
      LBMReal rhoLB = 0.0;
      LBMReal nuLB = 8.66025e-6;
      //////////////////////////////////////////////////////////////////////////
      ////dx_coarse = 0.0015 mm
      //string  pathname = "d:/temp/AcousticPulseXZ-4th-0.0015";
      //double  endTime = 40;
      //double  outTime = 40;
      //LBMReal dx =  0.0015;
      //LBMReal rhoLB = 0.0;
      //LBMReal nuLB = 8.66025e-6*2.0;
      ////////////////////////////////////////////////////////////////////////////
      ////dx_coarse = 0.00075 mm
      //string  pathname = "d:/temp/AcousticPulseXZ-4th-0.00075";
      //double  endTime = 80;
      //double  outTime = 80;
      //LBMReal dx =  0.00075;
      //LBMReal rhoLB = 0.0;
      //LBMReal nuLB = 8.66025e-6*4.0;
      //////////////////////////////////////////////////////////////////////////

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      int baseLevel = 0;
      int refineLevel = 1;

      //bounding box
      double g_minX1 = -0.06;
      double g_minX2 = -0.06;
      double g_minX3 = -0.06;

      double g_maxX1 = 0.06;
      double g_maxX2 = 0.06;
      double g_maxX3 = 0.06;

      //double g_minX1 = -1;
      //double g_minX2 = -1;
      //double g_minX3 = -1;

      //double g_maxX1 = 1;
      //double g_maxX2 = 1;
      //double g_maxX3 = 1;

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

      SPtr<GbObject3D> refCube(new GbCuboid3D(-0.02,-0.02,-0.02,0.02,0.02,0.02));
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

      double bulckViscosity = 10.0*nuLB;
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      //dynamicPointerCast<CompressibleCumulant4thOrderViscosityLBMKernel>(kernel)->setBulkViscosity(bulckViscosity);
      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new K15CompressibleNavierStokes());
      //dynamicPointerCast<K15CompressibleNavierStokes>(kernel)->setBulkOmegaToOmega(true);
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

      //set connectors  
     //SPtr<Interpolator> iProcessor(new CompressibleOffsetInterpolator());
      SPtr<Interpolator> iProcessor(new CompressibleOffsetMomentsInterpolator());
      //dynamicPointerCast<CompressibleOffsetMomentsInterpolator>(iProcessor)->setBulkViscosity(nuLB, bulckViscosity);
      SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);

      UBLOG(logINFO, "SetConnectorsBlockVisitor:start");
      grid->accept(setConnsVisitor);
      UBLOG(logINFO, "SetConnectorsBlockVisitor:end");

      mu::Parser fctRoh;
      //z
      //fctRoh.SetExpr("epsilon*exp(-alpha*(x1*x1+x2*x2))");
      //x
      //fctRoh.SetExpr("epsilon*exp(-alpha*(x3*x3+x2*x2))");
      //y
      fctRoh.SetExpr("epsilon*exp(-alpha*scaleFactor*(x3*x3+x1*x1))");
      //fctRoh.SetExpr("epsilon*exp(-alpha*(x3*x3+x1*x1))");

      fctRoh.DefineConst("epsilon", 1e-3);
      fctRoh.DefineConst("alpha", log(2.0)/(0.01));
      fctRoh.DefineConst("scaleFactor", 277.777777779);
      //fctRoh.SetExpr("x1*0.001");

      //initialization of distributions
      InitDistributionsBlockVisitor initVisitor;
      initVisitor.setRho(fctRoh);
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

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(writeMQCoProcessor);

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


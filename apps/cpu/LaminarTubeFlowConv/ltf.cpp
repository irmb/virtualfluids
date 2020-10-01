#include <iostream>
#include <string>

#include "vfluids.h"

using namespace std;


int x[3] = { 120, 240, 480 };
int y[3] = { 20, 40, 80 };
int z[3] = { 20, 40, 80 };

//int x[3] = { 120, 120, 120 };
//int y[3] = { 20, 20, 20 };
//int z[3] = { 20, 20, 20 };
double nuLB = 0.001;
double dp[3] = { 25000.0, 100000.0, 400000.0 };
double tout[3] = { 4000.0, 16000.0, 64000.0 };
double tend[3] = { 100001.0, 400001.0, 1600001.0 };
//double deltax[3] = { 1.0, 0.5, 0.25 };
double deltax[3] = { 1.0, 1.0, 1.0 };


void run(int tn)
{
   try
   {
      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
      int numOfThreads = 1;
      double availMem = 0;

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if(machine == "BOMBADIL") 
      {
         pathname = "d:/temp/ltfc" + UbSystem::toString(tn);
         numOfThreads = 1;
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         pathname = "/work/koskuche/scratch/ltfc"+UbSystem::toString(tn);
         numOfThreads = 8;
         availMem = 12.0e9;

#if defined(__unix__)
         if (myid == 0)
         {
            const char* str = pathname.c_str();
            int status = mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         }
#endif 

         if(myid ==0)
         {
            stringstream logFilename;
            logFilename << pathname + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + "_" + UbSystem::toString(myid) + ".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      double dx = deltax[tn];

      double L1 = x[tn];
      double L2 = y[tn];
      double L3 = z[tn];

      LBMReal dLB = L2;
      LBMReal rhoLB = 0.0;
      LBMReal l = L2 / dx;


      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = 0;

      //bounding box
      double g_minX1 = 0.0;
      double g_minX2 = -L2 / 2.0;
      double g_minX3 = -L3 / 2.0;

      double g_maxX1 = L1;
      double g_maxX2 = L2 / 2.0;
      double g_maxX3 = L3 / 2.0;

      //obstacle
      GbObject3DPtr cylinder(new GbCylinder3D(g_minX1-2.0*dx, 0.0, 0.0, g_maxX1+2.0*dx, 0.0, 0.0, dLB/2.0));
      GbSystem3D::writeGeoObject(cylinder.get(),pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

      double offs = dx;
      //GbObject3DPtr gridCube(new GbCuboid3D(g_minX1-offs, g_minX2-offs, g_minX3-offs, g_maxX1+offs, g_maxX2+offs, g_maxX3+offs));
      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

      const int blocknx1 = 10;
      const int blocknx2 = 10;
      const int blocknx3 = 10;
      
      double blockLength = blocknx1*dx;

      Grid3DPtr grid(new Grid3D(comm));
      UbSchedulerPtr rSch(new UbScheduler(100000, 100000));
      RestartPostprocessor rp(grid, rSch, comm, pathname, RestartPostprocessor::BINARY);

         if(myid ==0)
         {
            UBLOG(logINFO,"L = " << l );
            UBLOG(logINFO,"lLB = " << L1 );
            UBLOG(logINFO,"rho = " << rhoLB );
            UBLOG(logINFO,"nue = " << nuLB );
            UBLOG(logINFO,"dx = " << dx );
            UBLOG(logINFO,"Preprozess - start");
         }

         grid->setDeltaX(dx);
         grid->setBlockNX(blocknx1, blocknx2, blocknx3);
         //grid->setPeriodicX3(true);

         if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
      
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1 - 2.0*dx, g_minX2 - dx, g_minX3 - dx, g_minX1, g_maxX2, g_maxX3 + dx));
          if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow (new GbCuboid3D(g_maxX1, g_minX2, g_minX3-dx, g_maxX1+2.0*dx, g_maxX2+dx, g_maxX3+dx));
          if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         ppblocks->update(0);
      
         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
         D3Q27InteractorPtr cylinderInt( new D3Q27Interactor(cylinder, grid, bcObst,Interactor3D::INVERSESOLID));

         D3Q27BoundaryConditionAdapterPtr denBCAdapter1(new D3Q27DensityBCAdapter(1.0/dp[tn]));
         D3Q27InteractorPtr inflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, denBCAdapter1,Interactor3D::SOLID));
         grid->addAndInitInteractor(inflowInt);

         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapter2(new D3Q27DensityBCAdapter(-1.0/dp[tn]));
         denBCAdapter2->setSecondaryBcOption(0);
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, denBCAdapter2, Interactor3D::SOLID));
         grid->addAndInitInteractor(outflowInt);

         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(cylinderInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         ppblocks->update(0);
         ppblocks.reset();

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx1)* (blocknx2)* (blocknx3);
         unsigned long nod = nob * (blocknx1)* (blocknx2)* (blocknx3);
         unsigned long nodg = nob * (blocknx1 + gl) * (blocknx2 + gl) * (blocknx3 + gl);
         double needMemAll = double(nodg*(27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
         double needMem = needMemAll / double(comm->getNumberOfProcesses());

         if (myid == 0)
         {
            UBLOG(logINFO, "Number of blocks = " << nob);
            UBLOG(logINFO, "Number of nodes  = " << nod);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for (int level = minInitLevel; level <= maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO, "Number of blocks for level " << level << " = " << nob);
               UBLOG(logINFO, "Number of nodes for level " << level << " = " << nob*nodb);
            }
            UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
            UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
            UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
         }

         LBMKernel3DPtr kernel;
         rhoLB = 0.0;
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));

         //
         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);

         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         intHelper.setBC();

         //initialization of distributions
         D3Q27ETInitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);

         mu::Parser vx1;
         vx1.DefineConst("dp", dp[tn]);
         vx1.DefineConst("L", x[tn]);
         vx1.DefineConst("nu", nuLB);
         vx1.DefineConst("r", z[tn] * 0.5); 
         vx1.DefineConst("c", 0.0);
         vx1.SetExpr("(2/dp)/(3.0*4.0*L*nu)*(r^2-((c-x2)^2+(c-x3)^2))");

         initVisitor.setVx1(vx1);
         grid->accept(initVisitor);


         //Postrozess
         {
            UbSchedulerPtr geoSch(new UbScheduler(1));
            D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
               new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
            ppgeo->update(0);
            ppgeo.reset();
         }

      if (myid == 0) UBLOG(logINFO, "Preprozess - end");

      double outTime = tout[tn];
      UbSchedulerPtr visSch(new UbScheduler(outTime));

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname, WbWriterVtkXmlASCII::getInstance(), conv);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 10, 30));
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      double endTime = tend[tn];
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      if(myid == 0) UBLOG(logINFO,"Simulation-start");
      calculation->calculate();
      if(myid == 0) UBLOG(logINFO,"Simulation-end");
   }
   catch(std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch(std::string& s)
   {
      cerr << s << endl;
   }
   catch(...)
   {
      cerr << "unknown exception" << endl;
   }

}
int main(int argc, char* argv[])
{

   run(UbSystem::stringTo<int>(argv[1]));

   return 0;
}


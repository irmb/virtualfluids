#include <iostream>
#include <string>

#include "vfluids.h"

using namespace std;


void run(const char *cstr)
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
         pathname = "d:/temp/ltf_ref_bc_long2";
         numOfThreads = 4;
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         pathname = "/work/koskuche/scratch/ltf2";
         numOfThreads = 8;
         availMem = 12.0e9;

         //if(myid ==0)
         //{
         //   stringstream logFilename;
         //   logFilename <<  pathname + "/logfile.txt";
         //   UbLog::output_policy::setStream(logFilename.str());
         //}
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      const double dx = 1.0;

      double L1 = 32*2;
      double L2 = 32;
      double L3 = L2;

      LBMReal dLB = L2;
      LBMReal uLB = 0.05;
      LBMReal Re = 10.0;
      LBMReal rhoLB = 0.0;
      LBMReal l = L2 / dx;
      LBMReal nuLB = (uLB*dLB)/Re;

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

      const int blocknx1 = 8;
      const int blocknx2 = 8;
      const int blocknx3 = 8;
      
      double blockLength = blocknx1*dx;

      Grid3DPtr grid(new Grid3D(comm));
      UbSchedulerPtr rSch(new UbScheduler(100000, 100000));
      RestartPostprocessor rp(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY);

         if(myid ==0)
         {
            UBLOG(logINFO,"L = " << l );
            UBLOG(logINFO,"v = " << uLB );
            UBLOG(logINFO,"rho = " << rhoLB );
            UBLOG(logINFO,"nue = " << nuLB );
            UBLOG(logINFO,"Re = " << Re );
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

         double r = boost::dynamic_pointer_cast<GbCylinder3D>(cylinder)->getRadius();
         double cx1 = g_minX1;
         double cx2 = cylinder->getX2Centroid();
         double cx3 = cylinder->getX3Centroid();
         mu::Parser fct;
         fct.SetExpr("vx1*(1-((x2-y0)^2+(x3-z0)^2)/(R^2))");
         fct.DefineConst("x2Vmax", 0.0); //x2-Pos fuer vmax
         fct.DefineConst("x3Vmax", 0.0); //x3-Pos fuer vmax
         fct.DefineConst("R", r);
         fct.DefineConst("vx1", uLB);
         fct.DefineConst("x0", cx1);
         fct.DefineConst("y0", cx2);
         fct.DefineConst("z0", cx3);
         fct.DefineConst("nue", nuLB);
         //mu::Parser fct;
         //fct.SetExpr("-vx1*((x2+h)*(x2-h))/(h^2)");
         //fct.DefineConst("vx1"  , uLB           );
         //fct.DefineConst("h"  ,    (g_maxX2 - g_minX2)/2.0        );

         D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter(true, false, false, fct, 0, D3Q27BCFunction::INFCONST));
         //velBCAdapter->setSecondaryBcOption(2);

         D3Q27InteractorPtr inflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));
         //D3Q27InteractorPtr inflowInt( new D3Q27Interactor(geoInflow, grid, bcObst,Interactor3D::SOLID));

         //D3Q27BoundaryConditionAdapterPtr denBCAdapter1(new D3Q27DensityBCAdapter(rhoLB*1.001));
         //D3Q27InteractorPtr inflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, denBCAdapter1,Interactor3D::SOLID));
         grid->addAndInitInteractor(inflowInt);

         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         denBCAdapter->setSecondaryBcOption(0);
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));
         //D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, bcObst,Interactor3D::SOLID));
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
         //initVisitor.setVx1(fct);
         //initVisitor.setNu(nuLB);
         //initVisitor.setVx1(0.01);
         //initVisitor.setVx2(0.02);
         //initVisitor.setVx3(0.03);
         grid->accept(initVisitor);

         //Postrozess
         //UbSchedulerPtr geoSch(new UbScheduler(1));
         //D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
         //   new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, comm, true));
         //ppgeo->update(0);
         //ppgeo.reset();

         {
            UbSchedulerPtr geoSch(new UbScheduler(1));
            //D3Q27MacroscopicQuantitiesPostprocessor ppgeo(grid,geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv,  comm, true);
            D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
               new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
            //grid->addObserver(ppgeo);
            grid->doPostProcess(0);
            //grid->notifyObservers(0);
            //grid->removeObserver(ppgeo);
         }

      if (myid == 0) UBLOG(logINFO, "Preprozess - end");

      double outTime = 1000;
      UbSchedulerPtr visSch(new UbScheduler(outTime));
      //visSch->addSchedule(1, 1, 10000);
      //visSch->addSchedule(10000, 10000, 50000);
      //visSch->addSchedule(1000, 1000, 100000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname, WbWriterVtkXmlASCII::getInstance(), conv);

      double endTime = 1000001.0;
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

   run(argv[1]);

   return 0;
}


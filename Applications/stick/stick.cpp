#include <iostream>
#include <string>

#include <boost/pointer_cast.hpp>

#include "vfluids.h"

using namespace std;



void main()
{

   try
   {
      string machine = QUOTEME(CAB_MACHINE);
      string pathname = "d:/temp/stick"; 
      int numOfThreads = 4;
      double availMem = 10e9;

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      double dx = 1;

      const int blocknx1 = 10;
      const int blocknx2 = 10;
      const int blocknx3 = 10;

      const int gridNx1 = 60;
      const int gridNx2 = 1;
      const int gridNx3 = 8;

      double L1 = gridNx1*blocknx1;
      double L2, L3;
      L2 = gridNx2*blocknx1;
      L3 = gridNx3*blocknx1;

      LBMReal radius = 1.0*dx;
      LBMReal uLB = 0.07;
      LBMReal Re = 1000.0;
      LBMReal rhoLB = 0.0;
      LBMReal nueLB = (uLB*1.0*radius)/Re;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      Grid3DPtr grid(new Grid3D(comm));
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);
      grid->setPeriodicX1(false);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(false);

      const int baseLevel = 0;
      const int refineLevel = 0;

      //bounding box
      double d_minX1 = 0.0;
      double d_minX2 = 0.0;
      double d_minX3 = 0.0;

      double d_maxX1 = L1;
      double d_maxX2 = L2;
      double d_maxX3 = L3;

      double blockLength = blocknx1*dx;

      if(myid ==0)
      {
         UBLOG(logINFO,"Parameters:");
         UBLOG(logINFO,"uLB = " << uLB );
         UBLOG(logINFO,"rhoLB = " << rhoLB );
         UBLOG(logINFO,"nueLB = " << nueLB );
         UBLOG(logINFO,"Re = " << Re );
         UBLOG(logINFO,"dx = " << dx );
         UBLOG(logINFO,"number of levels = " << refineLevel+1 );
         UBLOG(logINFO,"numOfThreads = " << numOfThreads );
         UBLOG(logINFO,"Preprozess - start");
      }

      GbObject3DPtr gridCube(new GbCuboid3D(d_minX1, d_minX2, d_minX3, d_maxX1, d_maxX2, d_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance()); 

      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      //cylinder
      //GbObject3DPtr cylinder(new GbCylinder3D(L1/4.0, -2.0, radius, L1/4.0, L2+2.0, radius, radius));
      //GbSystem3D::writeGeoObject(cylinder.get(),pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

      GbCuboid3DPtr stick (new GbCuboid3D(L1/4.0, -2.0, 0.0, L1/4.0+150.0, L2+2.0, radius*3.0));
      if(myid == 0) GbSystem3D::writeGeoObject(stick.get(), pathname+"/geo/stick", WbWriterVtkXmlASCII::getInstance());

      //walls
      GbCuboid3DPtr addWallZmin (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_minX3));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallZmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_maxX3, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

      //inflow
      GbCuboid3DPtr geoInflow (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_minX1, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

      //outflow
      GbCuboid3DPtr geoOutflow (new GbCuboid3D(d_maxX1, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));


      //cylinder
      int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
      D3Q27BoundaryConditionAdapterPtr noSlip(new D3Q27NoSlipBCAdapter(bbOption));
      D3Q27InteractorPtr cylinderInt = D3Q27InteractorPtr ( new D3Q27Interactor(stick, grid, noSlip,Interactor3D::SOLID));

      //walls
      D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, noSlip,Interactor3D::SOLID));
      D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, noSlip,Interactor3D::SOLID));

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);

      //inflow
      D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
      velBCAdapter->setSecondaryBcOption(2);
      D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

      //outflow
      D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
      D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));

      Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
      InteractorsHelper intHelper(grid, metisVisitor);
      intHelper.addInteractor(cylinderInt);
      intHelper.addInteractor(addWallZminInt);
      intHelper.addInteractor(addWallZmaxInt);
      intHelper.addInteractor(inflowInt);
      intHelper.addInteractor(outflowInt);
      intHelper.selectBlocks();

      //domain decomposition for threads
      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      grid->accept(pqPartVisitor);

      ppblocks->update(0);
      ppblocks.reset();

      //set connectors
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
      grid->accept( setConnsVisitor );

      unsigned long nob = grid->getNumberOfBlocks();
      int gl = 3;
      unsigned long nod = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);

      double needMemAll  = double(nod*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
      double needMem  = needMemAll / double(comm->getNumberOfProcesses());

      if(myid == 0)
      {
         UBLOG(logINFO,"Number of blocks = " << nob);
         UBLOG(logINFO,"Number of nodes  = " << nod);
         UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
         UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
         UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
      }            

      LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));

      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
      grid->accept(kernelVisitor);

      intHelper.setBC();

      //initialization of distributions
      D3Q27ETInitDistributionsBlockVisitor initVisitor(nueLB, rhoLB);
      initVisitor.setVx1(fct);
      grid->accept(initVisitor);

      //Postrozess
      UbSchedulerPtr geoSch(new UbScheduler(1));
      D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
         new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
      ppgeo->update(0);
      ppgeo.reset();

      if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 

      UbSchedulerPtr stepSch(new UbScheduler(10000));
      //stepSch->addSchedule(1000, 0, 1000000);
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv);

      //InSituVTKPostprocessor isp(grid, stepSch, "d:/Data/insituDemo/metafile.csv", conv);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);


      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, 1000000, stepSch));

      if(myid == 0) 
         UBLOG(logINFO,"Simulation-start");

      calculation->calculate();

      if(myid == 0) 
         UBLOG(logINFO,"Simulation-end");

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



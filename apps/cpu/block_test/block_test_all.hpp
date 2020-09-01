#include <iostream>
#include <string>

#include "geometry3d/CoordinateTransformation3D.h"
#include "Grid3D.h"
#include "GenBlocksGridVisitor.h"
#include "geometry3d/GbSystem3D.h"
#include "geometry3d/GbCuboid3D.h"
#include "geometry3d/GbCylinder3D.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "OverlapBlockVisitor.h"
#include "RefineInterGbObjectsVisitor.h"
#include "SetKernelBlockVisitor.h"
#include "LBMKernelETD3Q27Cascaded.h"
#include "D3Q27MacroscopicQuantitiesPostprocessor.h"
#include "MPICommunicator.h"
#include "D3Q27ETBCProcessor.h"
#include "SimulationParameters.h"
#include "D3Q27SetUndefinedNodesBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include "D3Q27SetConnectorsBlockVisitor.h"
#include "NullCommunicator.h"
#include "D3Q27ETInitDistributionsBlockVisitor.h"
#include "CalculationManager.h"
#include "PQueuePartitioningGridVisitor.h"
#include "MetisPartitioningGridVisitor.h"
#include "D3Q27Interactor.h"
#include "D3Q27NoSlipBCAdapter.h"
#include "D3Q27VelocityBCAdapter.h"
#include "D3Q27DensityBCAdapter.h"
#include "D3Q27BoundaryConditionAdapter.h"
#include "StringUtil.hpp"
#include "D3Q27OffsetInterpolationProcessor.h"
#include "D3Q27CompactInterpolationProcessor.h"
#include "D3Q27PressureDifferencePostprocessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "RestartPostprocessor.h"
#include "SolidBlocksHelper.h"
#include "NUPSCounterPostprocessor.h"
#include "BlocksPostprocessor.h"
#include "LBMKernelETD3Q27BGK.h"
#include "EmergencyExitPostprocessor.h"
#include "D3Q27ForcesPostprocessor.h"
#include "ConfigFileReader.h"
#include "MemoryUtil.h"
#include "LBMKernelETD3Q27CCLB.h"


using namespace std;


void block_test_all(const char *cstr)
{
   try
   {
      ConfigFileReader cf(cstr);
      if ( !cf.read() )
      {
         std::string exceptionText = "Unable to read configuration file\n";
         throw exceptionText;
      }

      string machine = QUOTEME(CAB_MACHINE);
      string pathname = cf.getValue("path"); 
      int numOfThreads = UbSystem::stringTo<int>(cf.getValue("numOfThreads"));
      double availMem = 0;

      CommunicatorPtr comm(new MPICommunicator());
      int myid = comm->getProcessID();

      if(machine == "BOMBADIL") 
      {
         //pathname = "c:/temp/block_test";
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         //pathname = "/work/koskuche/scratch/block_test";
         availMem = 12.0e9;

         if(myid ==0)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      double dx = 1;

      const int blocknx1 = UbSystem::stringTo<int>(cf.getValue("blocknx1")); //16;
      const int blocknx2 = UbSystem::stringTo<int>(cf.getValue("blocknx2"));//16;
      const int blocknx3 = UbSystem::stringTo<int>(cf.getValue("blocknx3"));//16;

      const int gridNx1 = UbSystem::stringTo<int>(cf.getValue("gridNx1"));//3;
      const int gridNx2 = UbSystem::stringTo<int>(cf.getValue("gridNx2"));//3;
      const int gridNx3 = UbSystem::stringTo<int>(cf.getValue("gridNx3"));//3;


      double L1 = gridNx1*blocknx1;
      double L2, L3, H;
      L2 = L3 = H = gridNx2*blocknx1;

      LBMReal radius = 3;
      LBMReal uLB = 0.05;
      LBMReal Re = 300.0;
      LBMReal rhoLB = 1.0;
      LBMReal l = L2 / dx;
      LBMReal nueLB = (((4.0/9.0)*uLB)*2.0*(radius/dx))/Re;


      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = UbSystem::stringTo<int>(cf.getValue("refineLevel"));

      //obstacle
      GbObject3DPtr cylinder(new GbCylinder3D(L1*0.5, L2*0.5, 0, L1*0.5, L2*0.5, L3, radius));
      GbSystem3D::writeGeoObject(cylinder.get(),pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

      D3Q27InteractorPtr cylinderInt;

      //bounding box
      double d_minX1 = 0.0;
      double d_minX2 = 0.0;
      double d_minX3 = 0.0;

      double d_maxX1 = L1;
      double d_maxX2 = L2;
      double d_maxX3 = L3;

      double offs = dx;

      //double g_minX1 = d_minX1-offs-0.499999*dx;
      double g_minX1 = d_minX1-offs;
      double g_minX2 = d_minX2-offs;
      double g_minX3 = d_minX3-offs;

      double g_maxX1 = d_maxX1+offs;
      double g_maxX2 = d_maxX2+offs;
      double g_maxX3 = d_maxX3+offs;

      double blockLength = blocknx1*dx;

      //refinement area
      double off = 1;
      GbObject3DPtr refineCube(new  GbCuboid3D(cylinder->getX1Minimum()-off, cylinder->getX2Minimum()-off, cylinder->getX3Minimum(), 
         cylinder->getX1Maximum()+off, cylinder->getX2Maximum()+off, cylinder->getX3Maximum()));

      Grid3DPtr grid(new Grid3D(blocknx1, blocknx2, blocknx3, gridNx1, gridNx2, gridNx3));


      if(myid ==0)
      {
         UBLOG(logINFO,"Parameters:");
         UBLOG(logINFO,"L = " << L2/dx );
         UBLOG(logINFO,"v = " << uLB );
         UBLOG(logINFO,"rho = " << rhoLB );
         UBLOG(logINFO,"nue = " << nueLB );
         UBLOG(logINFO,"Re = " << Re );
         UBLOG(logINFO,"dx = " << dx );
         UBLOG(logINFO,"number of levels = " << refineLevel+1 );
         UBLOG(logINFO,"numOfThreads = " << numOfThreads );
         UBLOG(logINFO,"Preprozess - start");
      }

      if(myid ==0) GbSystem3D::writeGeoObject(refineCube.get(),pathname + "/geo/refineCube", WbWriterVtkXmlBinary::getInstance());

      //walls
      GbCuboid3DPtr addWallYmin (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_minX2, d_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallZmin (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_minX3));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallYmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_maxX2, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallZmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_maxX3, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

      //inflow
      GbCuboid3DPtr geoInflow (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_minX1, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

      //outflow
      GbCuboid3DPtr geoOutflow (new GbCuboid3D(d_maxX1, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));

      if (refineLevel > 0)
      {
         if(myid == 0) UBLOG(logINFO,"Refinement - start");	
         RefineCrossAndInsideGbObjectBlockVisitor refVisitor(refineCube, baseLevel, refineLevel-1);
         grid->accept(refVisitor);

         RatioBlockVisitor ratioVisitor(refineLevel);
         grid->accept(ratioVisitor);

         RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
         grid->accept(ratioSmoothVisitor);

         OverlapBlockVisitor overlapVisitor(refineLevel);
         grid->accept(overlapVisitor);

         std::vector<int> dirs;
         D3Q27System::getLBMDirections(dirs);
         SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
         grid->accept(interDirsVisitor);
         if(myid == 0) UBLOG(logINFO,"Refinement - end");	
      }

      MetisPartitioningGridVisitor metisVisitor(numOfThreads, D3Q27System::B, comm, true);
      grid->accept( metisVisitor );

      SolidBlocksHelper sd(grid, comm);

      int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
      D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
      cylinderInt = D3Q27InteractorPtr ( new D3Q27Interactor(cylinder, grid, bcObst,Interactor3D::SOLID));

      //walls
      D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcObst,Interactor3D::SOLID));
      D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcObst,Interactor3D::SOLID));
      D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcObst,Interactor3D::SOLID));
      D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcObst,Interactor3D::SOLID));

      mu::Parser fct;
      fct.SetExpr("16*U*x2*x3*(H-x2)*(H-x3)/H^4");
      fct.DefineConst("U", uLB);
      fct.DefineConst("H", H);

      //inflow
      D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
      velBCAdapter->setSecondaryBcOption(2);
      D3Q27BoundaryConditionAdapterPtr denBCAdapterFront(new D3Q27DensityBCAdapter(rhoLB));
      denBCAdapterFront->setSecondaryBcOption(1);
      D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));
      //D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, bcObst, Interactor3D::SOLID));

      //outflow
      D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
      denBCAdapter->setSecondaryBcOption(1);
      D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));
      //D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, bcObst,Interactor3D::SOLID));

      //sd.addInteractor(cylinderInt);
      sd.addInteractor(addWallYminInt);
      sd.addInteractor(addWallZminInt);
      sd.addInteractor(addWallYmaxInt);
      sd.addInteractor(addWallZmaxInt);
      sd.addInteractor(inflowInt);
      sd.addInteractor(outflowInt);

      sd.deleteSolidBlocks();

      grid->accept( metisVisitor );

      //set connectors
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
      grid->accept( setConnsVisitor );

      //domain decomposition for threads
      //PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      //grid->accept(pqPartVisitor);

      ppblocks->update(0);
      ppblocks.reset();

      unsigned long nob = grid->getNumberOfBlocks();
      int gl = 3;
      unsigned long nodb = (blocknx1) * (blocknx2) * (blocknx3);
      unsigned long nod = nob * (blocknx1) * (blocknx2) * (blocknx3);
      unsigned long nodg = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);
      double needMemAll  = double(nodg*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
      double needMem  = needMemAll / double(comm->getNumberOfProcesses());

      if(myid == 0)
      {
         UBLOG(logINFO,"Number of blocks = " << nob);
         UBLOG(logINFO,"Number of nodes  = " << nod);
         int minInitLevel = grid->getCoarsestInitializedLevel();
         int maxInitLevel = grid->getFinestInitializedLevel();
         for(int level = minInitLevel; level<=maxInitLevel; level++)
         {
            int nobl = grid->getNumberOfBlocks(level);
            UBLOG(logINFO,"Number of blocks for level " << level <<" = " << nob);
            UBLOG(logINFO,"Number of nodes for level " << level <<" = " << nob*nodb);
         }
         UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
         UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
         UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
      }            

      int kernelType = UbSystem::stringTo<int>(cf.getValue("kernel"));
      LBMKernel3DPtr kernel;
      if (kernelType == 0)
      {
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));
      }
      else if (kernelType == 1)
      {
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3));
      }
      else if (kernelType == 2)
      {
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3));
      }

      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
      grid->accept(kernelVisitor);

      if (refineLevel > 0)
      {
         D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
         grid->accept(undefNodesVisitor);
      }

      //walls
      grid->addAndInitInteractor(addWallYminInt);
      grid->addAndInitInteractor(addWallZminInt);
      grid->addAndInitInteractor(addWallYmaxInt);
      grid->addAndInitInteractor(addWallZmaxInt);

      //obstacle
      //grid->addAndInitInteractor(cylinderInt);

      //inflow
      grid->addAndInitInteractor(inflowInt);

      //outflow
      grid->addAndInitInteractor(outflowInt);

      //initialization of distributions
      D3Q27ETInitDistributionsBlockVisitor initVisitor(1.0);
      //initVisitor.setVx1(fct);
      grid->accept(initVisitor);

      //Postrozess
      UbSchedulerPtr geoSch(new UbScheduler(1));
      D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
         new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, comm, true));
      ppgeo->update(0);
      ppgeo.reset();

      if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      double outTime = 1.0;
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      //UbSchedulerPtr stepSch(new UbScheduler());
      //stepSch->addSchedule(10, 100, 1000);
      //nodeSch->addSchedule(1000, 1000, 10000);
      //nodeSch->addSchedule(10000, 10000, 50000);
      //stepSch->addSchedule(100, 100, 1000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv, comm);

      UbSchedulerPtr visSch(new UbScheduler());
      //UbSchedulerPtr visSch(stepSch);
      double endTime = UbSystem::stringTo<int>(cf.getValue("endTime"));//10001.0;

      //cout << "PID = " << myid << " Total Physical Memory (RAM): " << MemoryUtil::getTotalPhysMem()<<endl;
      //cout << "PID = " << myid << " Physical Memory currently used: " << MemoryUtil::getPhysMemUsed()<<endl;
      //cout << "PID = " << myid << " Physical Memory currently used by current process: " << MemoryUtil::getPhysMemUsedByMe()<<endl;

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
#include <iostream>
#include <string>

#include <boost/pointer_cast.hpp>

#include "vfluids.h"

using namespace std;


void setup(const char *cstr1, const char *cstr2)
{
   try
   {
      //Sleep(30000);

      ConfigFileReader cf(cstr1);
      if ( !cf.read() )
      {
         std::string exceptionText = "Unable to read configuration file\n";
         throw exceptionText;
      }

      //parameters from config file
      string machine = cf.getValue("machine");
      string pathname = cf.getValue("path");
      string geoFile = cf.getValue("geoFile");
      int numOfThreads = UbSystem::stringTo<int>(cf.getValue("numOfThreads"));
      double availMem = UbSystem::stringTo<double>(cf.getValue("availMem"));
      int refineLevel = UbSystem::stringTo<int>(cf.getValue("refineLevel"));
      int blocknx = UbSystem::stringTo<int>(cf.getValue("blocknx"));

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if(machine == "Bombadil") int dumy=0; 
      else if(machine == "Ludwig" || machine == "HLRN")      
      {
         if(myid ==0)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown machine");

      GbTriFaceMesh3DPtr geo (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(geoFile,"geo"));
      if(myid == 0) GbSystem3D::writeGeoObject(geo.get(), pathname+"/geo/geo", WbWriterVtkXmlASCII::getInstance());

      double dx = (fabs(geo->getX3Maximum()-geo->getX3Minimum())*10e-3)*(double)(1<<refineLevel);
      dx /= 4.0;

      double blockLength = blocknx*dx;

      double offsetX1 = fabs(geo->getX1Maximum()-geo->getX1Minimum());
      double h = fabs(geo->getX3Maximum()-geo->getX3Minimum());
      double offsetX2 = fabs(geo->getX2Maximum()-geo->getX2Minimum())/3.0;
      double offsetX3 = 3.0*h; //30.0*h;

      double g_minX1 = geo->getX1Minimum()-offsetX1;
      double g_minX2 = geo->getX2Minimum()+offsetX2;
      double g_minX3 = geo->getX3Centroid()-offsetX3;

      double g_maxX1 = geo->getX1Maximum()+5.0*offsetX1;
      double g_maxX2 = g_minX2 + 4.0*blockLength; 
      double g_maxX3 = geo->getX3Centroid()+offsetX3;

      //##########################################################################
      //## physical parameters
      //##########################################################################
      double Re            = 1e6;

      double rhoLB    = 0.0;
      double rhoReal  = 1.0;
      double nueReal  = 0.000015;//0.015;

      double lReal    =  3.0;//<-m     ;//Profile laenge in cm(! cm nicht m !)
      double uReal    = Re*nueReal/lReal;

      //##Machzahl:
      //#Ma     = uReal/csReal
      double Ma      = 0.1;//Ma-Real!
      double csReal  = uReal/Ma;
      double hLB     = lReal/dx;

      LBMUnitConverter unitConverter(lReal, csReal, rhoReal, hLB);

      double uLB     = uReal   * unitConverter.getFactorVelocityWToLb();
      double nueLB   = nueReal * unitConverter.getFactorViscosityWToLb();

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;

      ////////////////////////////////////////////////////////////////////////
      //Grid
      //////////////////////////////////////////////////////////////////////////
      Grid3DPtr grid(new Grid3D(comm));
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx, blocknx, blocknx);
      
      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      //gridCube->setCenterCoordinates(geo->getX1Centroid(), geo->getX2Centroid(), geo->getX3Centroid());
      if(myid == 0) GbSystem3D::writeGeoObject(gridCube.get(),pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());
      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      grid->setPeriodicX2(true);
      grid->setPeriodicX3(true);

      double outTime = 1.0;
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      //PostprocessorPtr pp(new D3Q27MacroscopicQuantitiesPostprocessor(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlASCII::getInstance(), conv, comm));

      UbSchedulerPtr rSch(new UbScheduler());
      rSch->addSchedule(50,50,50);
      RestartPostprocessorPtr rp(new RestartPostprocessor(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::TXT));
      

      std::string opt;

      if(cstr2!= NULL)
         opt = std::string(cstr2);

      if/*(cstr== NULL)*/(cstr2!= NULL)
      {
         if(myid==0) UBLOG(logINFO,"Restart step: " << opt);
         grid = rp->restart(UbSystem::stringTo<int>(opt));
         rp->reconnect(grid);

         SetForcingBlockVisitor forcingVisitor(0.0, 0.0, 0.0);
         grid->accept(forcingVisitor);

         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
      }
      else
{
      //rp->addPostprocessor(pp);
      if(myid ==0)
      {
         UBLOG(logINFO,"Parameters:");
         UBLOG(logINFO, "* Re            ="<<Re);
         UBLOG(logINFO, "* Ma            ="<<Ma);
         UBLOG(logINFO, "* uReal         ="<<uReal);
         UBLOG(logINFO, "* nueReal       ="<<nueReal);
         UBLOG(logINFO, "* nue           ="<<nueLB);
         UBLOG(logINFO, "* velocity      ="<<uLB);
         //UBLOG(logINFO, "* LX1 (world/LB)="<<kanallaengeSI<<"/"<<kanallaengeSI/coarseNodeDx);
         //UBLOG(logINFO, "* LX2 (world/LB)="<<kanalbreiteSI<<"/"<<kanalbreiteSI/coarseNodeDx);
         //UBLOG(logINFO, "* LX3 (world/LB)="<<kanalhoeheSI<<"/"<<kanalhoeheSI/coarseNodeDx);
         UBLOG(logINFO, "* dx_base       ="<<dx);
         UBLOG(logINFO, "* dx_refine     ="<<dx/(double)(1<<refineLevel));
         //UBLOG(logINFO, "* nx1/2/3       ="<<nx[0]<<"/"<<nx[1]<<"/"<<nx[2]);
         UBLOG(logINFO, "* blocknx1/2/3  ="<<blocknx<<"/"<<blocknx<<"/"<<blocknx);
         //UBLOG(logINFO, "* x2Periodic    ="<<periodicx2);
         //UBLOG(logINFO, "* x3Periodic    ="<<periodicx3);
         UBLOG(logINFO, "*****************************************");
         UBLOGML(logINFO, "UnitConverter:"<<unitConverter.toString());
         UBLOG(logINFO, "*****************************************");    
         UBLOG(logINFO,"number of levels = " << refineLevel+1 );
         UBLOG(logINFO,"numOfThreads     = " << numOfThreads );
         UBLOG(logINFO,"Preprozess - start");
      }


      //inflow
      GbCuboid3DPtr geoInflow (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_minX3-4.0*blockLength, g_minX1+2.0*dx, g_maxX2+4.0*blockLength, g_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

      //outflow
      GbCuboid3DPtr geoOutflow (new GbCuboid3D(g_maxX1-2.0*dx, g_minX2-4.0*blockLength, g_minX3-4.0*blockLength, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_maxX3+4.0*blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));

      double scaleFactorX = 1.2;
      double scaleFactorZ = 1.2;
      //geo->scale(scaleFactorX, 1.0, scaleFactorZ);
      if(myid == 0) GbSystem3D::writeGeoObject(geo.get(), pathname+"/geo/geo2", WbWriterVtkXmlASCII::getInstance());

      int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
      D3Q27BoundaryConditionAdapterPtr noSlipBCAdapter(new D3Q27NoSlipBCAdapter(bbOption));

      Interactor3DPtr geoIntr = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(geo, grid, noSlipBCAdapter,Interactor3D::SOLID));

      //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(geoIntr)->refineBlockGridToLevel(refineLevel, 0.0, 5.0);

      if (refineLevel > 0)
      {
         if(myid == 0) UBLOG(logINFO,"Refinement - start");	
         //RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
         //refineHelper.addGbObject(geo, refineLevel);
         //refineHelper.refine();
         RefineAroundGbObjectHelper refineHelper(grid, refineLevel, boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(geoIntr), 0.0, 0.5);
         refineHelper.refine();
         if(myid == 0) UBLOG(logINFO,"Refinement - end");	
      }

      ppblocks->update(0);
      ppblocks.reset();
      return;

      //geo->scale(1.0/scaleFactorX, 1.0, 1.0/scaleFactorX);
      //geo = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(geoFile,"geo"));
      if(myid == 0) GbSystem3D::writeGeoObject(geo.get(), pathname+"/geo/geo3", WbWriterVtkXmlASCII::getInstance());

      MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, true, numOfThreads);
      grid->accept( metisVisitor );

      SolidBlocksHelper sd(grid, comm);

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);

      //inflow
      D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
      velBCAdapter->setSecondaryBcOption(2);
      D3Q27InteractorPtr inflowIntr  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

      //outflow
      D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
      denBCAdapter->setSecondaryBcOption(0);
      D3Q27InteractorPtr outflowIntr = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));


      sd.addInteractor(inflowIntr);
      sd.addInteractor(outflowIntr);
      sd.addInteractor(geoIntr);

      sd.deleteSolidBlocks();

      grid->accept( metisVisitor );

      sd.setTransBlocks();

      ppblocks->update(0);
      ppblocks.reset();

      //set connectors
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
      grid->accept( setConnsVisitor );

      //domain decomposition for threads
      //PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      //grid->accept(pqPartVisitor);


      unsigned long nob = grid->getNumberOfBlocks();
      int gl = 3;
      unsigned long nodb = (blocknx) * (blocknx) * (blocknx);
      unsigned long nod = nob * (blocknx) * (blocknx) * (blocknx);
      unsigned long nodg = nob * (blocknx+gl) * (blocknx+gl) * (blocknx+gl);
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
            UBLOG(logINFO,"Number of blocks for level " << level <<" = " << nobl);
            UBLOG(logINFO,"Number of nodes for level " << level <<" = " << nobl*nodb);
         }
         UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
         UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
         UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
      }            

      LBMKernel3DPtr kernel;
      kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx, blocknx, blocknx, LBMKernelETD3Q27CCLB::NORMAL));

      //mu::Parser fctForcingX1;
      //fctForcingX1.SetExpr("Fx1");
      //fctForcingX1.DefineConst("Fx1", 9.99685e-7);

      //kernel->setForcingX1(fctForcingX1);
      //kernel->setWithForcing(true);
      //
      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
      grid->accept(kernelVisitor);

      if (refineLevel > 0)
      {
         D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
         grid->accept(undefNodesVisitor);
      }

	  //UbSchedulerPtr geoSch(new UbScheduler(1));
	  //D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
		 // new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, comm, true));
	  //ppgeo->update(0);
	  //ppgeo.reset();

	  //return;

      //inflow
      grid->addAndInitInteractor(inflowIntr);

      //outflow
      grid->addAndInitInteractor(outflowIntr);

      //geo
      grid->addAndInitInteractor(geoIntr);

      //initialization of distributions
      D3Q27ETInitDistributionsBlockVisitor initVisitor(nueLB, rhoLB);
      initVisitor.setVx1(fct);
      initVisitor.setNu(nueLB);
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
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, true));
         //grid->addObserver(ppgeo);
         grid->doPostProcess(0);
         //grid->notifyObservers(0);
         //grid->removeObserver(ppgeo);
      }

      //grid->notifyObservers(0);

      //UbSchedulerPtr stepSch(new UbScheduler(outTime));
      D3Q27MacroscopicQuantitiesPostprocessorPtr pp(new D3Q27MacroscopicQuantitiesPostprocessor(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlASCII::getInstance(), conv));
      rp->addPostprocessor(pp);

      if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 
}
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

     // double outTime = 3.0;
     // UbSchedulerPtr stepSch(new UbScheduler(outTime));
      //UbSchedulerPtr stepSch(new UbScheduler());
      //stepSch->addSchedule(10, 100, 1000);
      //nodeSch->addSchedule(1000, 1000, 10000);
      //nodeSch->addSchedule(10000, 10000, 50000);
      //stepSch->addSchedule(100, 100, 1000);

      //UbSchedulerPtr st(new UbScheduler(100,50,1000));
      //UbSchedulerPtr rs(new UbScheduler(3));
      //AverageValuesPostprocessor ap(grid, pathname + "/av/av", WbWriterVtkXmlASCII::getInstance(), stepSch, rs, comm);

      //D3Q27ShearStressPostprocessor shs(grid,pathname + "/shs/shs", WbWriterVtkXmlASCII::getInstance(), stepSch, rs, comm);
      //shs.addInteractor(boost::dynamic_pointer_cast<D3Q27Interactor>(addWallZminInt));

      //D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlASCII::getInstance(), conv, comm);

      UbSchedulerPtr visSch(new UbScheduler(1));
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

int main(int argc, char* argv[])
{

   if ( argv != NULL )
   {
      if (argc > 1)
      {
         setup(argv[1], argv[2]);
      }
      else
      {
         cout << "Configuration file must be set!: " <<  argv[0] << " <config file>" << endl << std::flush;
      }
   }

   return 0;
} 


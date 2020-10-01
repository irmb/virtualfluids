

#include <iostream>
#include <string>
#include <math.h> 

#include <vfluids.h>

using namespace std;


void run(string configname)
{
   try
   {
      Configuration   config;
      config.load(configname);

      string          pathname          = config.getString("pathname");
      string          pathGeo           = config.getString("pathGeo");
      int             numOfThreads      = config.getInt("numOfThreads");
      string          diskFilename      = config.getString("diskFilename");
      string          mastFilename      = config.getString("mastFilename");
      vector<int>     blocknx           = config.getVector<int>("blocknx");
      double          restartStep       = config.getDouble("restartStep");
      double          restartStepStart  = config.getDouble("restartStepStart");
      double          endTime           = config.getDouble("endTime");
      double          outTime           = config.getDouble("outTime");
      double          availMem          = config.getDouble("availMem");
      bool            logToFile         = config.getBool("logToFile");
      vector<double>  geoLength         = config.getVector<double>("geoLength");
      int             refineLevel       = config.getInt("refineLevel");
      double          Re                = config.getDouble("Re");
      double          u_LB              = config.getDouble("u_LB");
      double          rho_LB            = config.getDouble("rho_LB");
      double          fineNodeDx        = config.getDouble("fineNodeDx");
      bool            restart           = config.getBool("restart");
      bool            averaging         = config.getBool("averaging");

      //UbLog::reportingLevel() = logDEBUG5;

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();


#if defined(__unix__)
      if (myid==0) 
      {
         const char* str = pathname.c_str();
         int status=mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      }
#endif 

      if(myid == 0 && logToFile)
      {
         stringstream logFilename;
         logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(myid)+".txt";
         UbLog::output_policy::setStream(logFilename.str());
      }

      if(myid==0) UBLOG(logINFO,"Test case: porous disk");

      int baseLevel = 0;

      double coarseNodeDx = fineNodeDx * (double)(1<<refineLevel);//geowerte

      double blockLengthx1 = blocknx[0]*coarseNodeDx; //geowerte
      double blockLengthx2 = blockLengthx1;
      double blockLengthx3 = blockLengthx1;

      bool periodicx1 = false;
      bool periodicx2 = true;
      bool periodicx3 = false;

      int sizeSP= (int)(500.0/coarseNodeDx); //500 mm sponge layer
      mu::Parser spongeLayer;
      spongeLayer.SetExpr("x1>=(sizeX-sizeSP)/dt ? (sizeX-x1)/sizeSP/2.0 + 0.5 : 1.0");
      spongeLayer.DefineConst("sizeX", 5000.0/coarseNodeDx);
      spongeLayer.DefineConst("sizeSP", sizeSP);

      //##########################################################################
      //## physical parameters
      //##########################################################################
      double nu_LB = (u_LB*(geoLength[2]/coarseNodeDx))/Re;

      LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter());
      
      Grid3DPtr grid(new Grid3D(comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(restartStep, restartStepStart));
      RestartPostprocessor rp(grid, rSch, comm, pathname, RestartPostprocessor::BINARY);
      //////////////////////////////////////////////////////////////////////////

      

      if (grid->getTimeStep() == 0)
      {
         if (myid==0) UBLOG(logINFO, "new start..");

         if (myid==0) UBLOG(logINFO, "load geometry start");
         GbTriFaceMesh3DPtr geoMast(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+mastFilename, "mast"));
         if (myid == 0) GbSystem3D::writeGeoObject(geoMast.get(), pathname + "/geo/geoMast", WbWriterVtkXmlBinary::getInstance());

         GbTriFaceMesh3DPtr geoDisk(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+diskFilename, "disk"));
         if (myid == 0) GbSystem3D::writeGeoObject(geoDisk.get(), pathname + "/geo/geoDisk", WbWriterVtkXmlBinary::getInstance());
         if (myid==0) UBLOG(logINFO, "load geometry end");

         //bounding box
         double g_minX1 = geoMast->getX1Centroid() - 2000.0;
         double g_minX2 = geoMast->getX2Centroid() - 1000.0;
         double g_minX3 = geoMast->getX3Minimum();

         double g_maxX1 = g_minX1 + geoLength[0];
         double g_maxX2 = g_minX2 + geoLength[1];
         double g_maxX3 = g_minX3 + geoLength[2];

         double nx1_temp = floor((g_maxX2-g_minX2) /(coarseNodeDx*(double)blocknx[1]));

         coarseNodeDx = (g_maxX2-g_minX2) /(nx1_temp*(double)blocknx[1]);

         fineNodeDx = coarseNodeDx / (double)(1<<refineLevel);

         //set grid
         grid->setDeltaX(coarseNodeDx);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
         grid->setPeriodicX1(periodicx1);
         grid->setPeriodicX2(periodicx2);
         grid->setPeriodicX3(periodicx3);


         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         //////////////////////////////////////////////////////////////////////////
         if (myid == 0)
         {
            UBLOG(logINFO, "*****************************************");
            UBLOG(logINFO, "* Parameters                            *");
            UBLOG(logINFO, "* Re                = "<<Re);
            UBLOG(logINFO, "* nu_LB             = "<<nu_LB);
            UBLOG(logINFO, "* u_LB              = "<<u_LB);
            UBLOG(logINFO, "* cdx               = "<<coarseNodeDx<<" mm");
            UBLOG(logINFO, "* fdx               = "<<fineNodeDx<<" mm");
            UBLOG(logINFO, "* nx1/2/3           = "<<grid->getNX1()<<"/"<<grid->getNX2()<<"/"<<grid->getNX3());
            UBLOG(logINFO, "* blocknx1/2/3      = "<<blocknx[0]<<"/"<<blocknx[1]<<"/"<<blocknx[2]);
            UBLOG(logINFO, "* x1Periodic        = "<<periodicx1);
            UBLOG(logINFO, "* x2Periodic        = "<<periodicx2);
            UBLOG(logINFO, "* x3Periodic        = "<<periodicx3);
            UBLOG(logINFO, "* number of levels  = "<<refineLevel+1);
            UBLOG(logINFO, "* path              = "<<pathname);

            UBLOG(logINFO, "*****************************************");
            UBLOG(logINFO, "* number of threads    = "<<numOfThreads);
            UBLOG(logINFO, "* number of processes  = "<<comm->getNumberOfProcesses());
            UBLOG(logINFO, "*****************************************");
            //UBLOGML(logINFO, "UnitConverter:"<<unitConverter->toString());
            //UBLOG(logINFO, "*****************************************");     
         }

         ////walls
         GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1-blockLengthx1, g_minX2-blockLengthx1, g_minX3-blockLengthx1, g_maxX1+blockLengthx1, g_minX2, g_maxX3+blockLengthx1));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1-blockLengthx1, g_maxX2, g_minX3-blockLengthx1, g_maxX1+blockLengthx1, g_maxX2+blockLengthx1, g_maxX3+blockLengthx1));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1 - blockLengthx1, g_minX2 - blockLengthx1, g_minX3 - blockLengthx1, g_maxX1 + blockLengthx1, g_maxX2 + blockLengthx1, g_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname + "/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1 - blockLengthx1, g_minX2 - blockLengthx1, g_maxX3, g_maxX1 + blockLengthx1, g_maxX2 + blockLengthx1, g_maxX3 + blockLengthx1));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname + "/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcNoSlip(new D3Q27NoSlipBCAdapter(bbOption));
         D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcNoSlip, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcNoSlip, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcNoSlip, Interactor3D::SOLID));
         D3Q27BoundaryConditionAdapterPtr bcSlip(new D3Q27SlipBCAdapter());
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcSlip, Interactor3D::SOLID));

         D3Q27TriFaceMeshInteractorPtr mastInt(new D3Q27TriFaceMeshInteractor(geoMast, grid, bcNoSlip, Interactor3D::SOLID, Interactor3D::POINTS));
         D3Q27TriFaceMeshInteractorPtr diskInt(new D3Q27TriFaceMeshInteractor(geoDisk, grid, bcNoSlip, Interactor3D::SOLID, Interactor3D::POINTS));

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLengthx1, g_minX2-blockLengthx1, g_minX3-blockLengthx1, g_minX1, g_maxX2+blockLengthx1, g_maxX3+blockLengthx1));
         if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr inflowInt(new D3Q27Interactor(geoInflow, grid, Interactor3D::SOLID));

         //inflow
         mu::Parser inflowProfile;
         //inflowProfile.SetExpr("u_ref*(((x3+Z_ref)/Z_ref)^a)");
         inflowProfile.SetExpr("u_ref");
         inflowProfile.DefineConst("u_ref", u_LB);
         inflowProfile.DefineConst("Z_ref", 300.0);
         inflowProfile.DefineConst("a", 0.143);

         D3Q27BoundaryConditionAdapterPtr velBCAdapter = D3Q27BoundaryConditionAdapterPtr(new D3Q27VelocityBCAdapter(true, false, false, inflowProfile, 0, D3Q27BCFunction::INFCONST));
         inflowInt->addBCAdapter(velBCAdapter);

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLengthx1, g_minX3-blockLengthx1, g_maxX1+blockLengthx1, g_maxX2+blockLengthx1, g_maxX3+blockLengthx1));
         if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rho_LB));
         D3Q27InteractorPtr outflowInt(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));

         {
            if (myid == 0) UBLOG(logINFO, "Write blocks - start"); 
            BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
            if (myid == 0)
               ppblocks->update(0);
            if (myid == 0) UBLOG(logINFO, "Write blocks - end"); 
         }
         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));	
        
         //////////////////////////////////////////////////////////////////////////
         //refinement
         double diameter = geoDisk->getLengthX2();
         GbCuboid3DPtr refineDiskBox(new GbCuboid3D(geoDisk->getX1Centroid()-0.2*diameter, geoDisk->getX2Centroid()-0.6*diameter, geoDisk->getX3Minimum()-diameter, 
            geoDisk->getX1Centroid() + 1.0*diameter, geoDisk->getX2Centroid()+0.6*diameter, geoDisk->getX3Maximum() + 0.05*diameter));
         if (myid == 0) GbSystem3D::writeGeoObject(refineDiskBox.get(), pathname + "/geo/refineDiskBox", WbWriterVtkXmlASCII::getInstance());

         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            refineHelper.addGbObject(refineDiskBox, refineLevel);
            refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         }

         {
            if (myid == 0) UBLOG(logINFO, "Write blocks - start");
            BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
            if (myid == 0)
               ppblocks->update(1);
            if (myid == 0) UBLOG(logINFO, "Write blocks - end");
         }
         ////////////////////////////////////////////
         /////delete solid blocks
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(mastInt);
         intHelper.addInteractor(diskInt);
         //intHelper.addInteractor(addWallYminInt);
         //intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(outflowInt);
         intHelper.addInteractor(inflowInt);
         intHelper.selectBlocks();
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	 
         //////////////////////////////////////

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx[0])* (blocknx[1])* (blocknx[2]);
         unsigned long nod = nob * (blocknx[0])* (blocknx[1])* (blocknx[2]);
         unsigned long nodg = nob * (blocknx[0] + gl) * (blocknx[1] + gl) * (blocknx[2] + gl);
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
               UBLOG(logINFO, "Number of blocks for level " << level << " = " << nobl);
               UBLOG(logINFO, "Number of nodes for level " << level << " = " << nobl*nodb);
            }
            UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
            UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
            UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
         }

         ////////////////////////////
         LBMKernel3DPtr kernel;
         //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));
         //with sponge layer
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLBWithSpongeLayer(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));
         kernel->setWithSpongeLayer(true);
         kernel->setSpongeLayer(spongeLayer);

         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         BoundaryConditionPtr densityBC(new NonEqDensityBoundaryCondition());
         BoundaryConditionPtr noSlipBC(new NoSlipBoundaryCondition());
         BoundaryConditionPtr velocityBC(new VelocityBoundaryCondition());
         BoundaryConditionPtr slipBC(new SlipBoundaryCondition());

         bcProc->addBC(densityBC);
         bcProc->addBC(noSlipBC);
         bcProc->addBC(velocityBC);
         bcProc->addBC(slipBC);
         kernel->setBCProcessor(bcProc);
         SetKernelBlockVisitor kernelVisitor(kernel, nu_LB, availMem, needMem);
         grid->accept(kernelVisitor);
         //////////////////////////////////
         //undef nodes
         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }
         //////////////////////////////////////////
         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept( setConnsVisitor );

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         intHelper.setBC();

         BoundaryConditionBlockVisitor bcVisitor;
         grid->accept(bcVisitor);

         //initialization of decompositions
         D3Q27ETInitDistributionsBlockVisitor initVisitor( nu_LB,rho_LB);
         //initVisitor.setVx1(inflowProfile);
         initVisitor.setVx1(u_LB);
         grid->accept(initVisitor);

         //Postprozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), 
            unitConverter, true));
         ppgeo->update(0);
         ppgeo.reset();
         geoSch.reset();

         if(myid == 0) UBLOG(logINFO,"Preprozess - end");      
      }
      else
      {
         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept( setConnsVisitor );
         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);
         //SetSpongeLayerBlockVisitor ssp(spongeLayer);
         //grid->accept(ssp);
         if(myid == 0) UBLOG(logINFO,"Restart - end"); 
      }
      UbSchedulerPtr visSch(new UbScheduler(outTime));

      double startStep = 80000;

      if (averaging)
      {
         UbSchedulerPtr resSchRMS(new UbScheduler());
         resSchRMS->addSchedule(100000, 80000, 10000000);
         UbSchedulerPtr resSchMeans(new UbScheduler());
         resSchMeans->addSchedule(100000, 80000, 10000000);
         UbSchedulerPtr stepAvSch(new UbScheduler());
         int averageInterval=100;
         stepAvSch->addSchedule(averageInterval, 0, 10000000);

         AverageValuesPostprocessor Avpp(grid, pathname, WbWriterVtkXmlBinary::getInstance(),
            visSch/*wann wird rausgeschrieben*/, stepAvSch/*wann wird gemittelt*/, resSchMeans, resSchRMS/*wann wird resettet*/, restart);
      }

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), unitConverter);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 10, 30));
      nupsSch->addSchedule(1000, 1000, 1000000000);
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      if(myid == 0)
      {
         UBLOG(logINFO,"PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      //double endTime = 80001;
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

   if (argv != NULL)
   {
      if (argv[1] != NULL)
      {
         run(string(argv[1]));
      }
      else
      {
         cout << "Configuration file is missing!" << endl;
      }
   }

   return 0;
}


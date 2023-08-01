#include <iostream>
#include <string>
#include <math.h> 

#include <vfluids.h>

using namespace std;

void run(const char *cstr1, const char *cstr2)
{
   try
   {
      string pathname; 
      string pathGeo;
      string pathLog;
      int numOfThreads = 1;
      bool logfile = false;
      stringstream logFilename;
      double availMem = 0;

      CommunicatorPtr comm = vf::parallel::MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      string machine = string(cstr1);

      if(machine == "my") 
      {
         //Sleep(30000);
         pathname = "d:/temp/town";
         pathGeo = "d:/Data/town";
         pathLog = "d:/temp/town";
         numOfThreads = 1;
         logfile = false;
         availMem = 15.0e9;
      }
      else if(machine == "Ludwig")      
      {
         pathname = "/work/koskuche/town";
         pathGeo = "/home/koskuche/data/town";
         pathLog = pathname;
         numOfThreads = 8;
         availMem = 12.0e9;///8*numOfThreads;
         logfile = true;
      }
      else if(machine == "HLRS")      
      {
         pathname = "/univ_1/ws1/ws/xrmkuchr-plate3-0";
         pathGeo = "/zhome/academic/HLRS/xrm/xrmkuchr/data/plate";
         pathLog = "/zhome/academic/HLRS/xrm/xrmkuchr/work/plate";
         numOfThreads = 16;
         availMem = 2.0e9;
         logfile = true;
      }
      else if(machine == "HLRN")      
      {
         pathname = "/gfs1/work/niivfcpu/scratch/plateEx";
         pathGeo = "/gfs1/work/niivfcpu/data/plate";
         pathLog = pathname;
         numOfThreads = 24;
         availMem = 64.0e9/24.0*numOfThreads;
         logfile = true;
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

#if defined(__unix__)
      if (myid==0) 
      {
         const char* str = pathLog.c_str();
         int status=mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      }
#endif 

      if(myid == 0 && logfile)
      {
         //UbLog::reportingLevel() = logDEBUG5;
         logFilename <<  pathLog + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(myid)+".txt";
         UbLog::output_policy::setStream(logFilename.str());
      }

      if(myid==0) UBLOG(logINFO,"Testcase town");

      //string townFilename = pathGeo + "/Manhattan.stl";
      string townFilename = pathGeo + "/town.stl"; 


      ///////////////Knotenabmessungen:
      int blocknx[3], nx[3];
      blocknx[0] = 8;
      blocknx[1] = 8;
      blocknx[2] = 8;

      nx[0] = 12;
      nx[1] = 12;
      nx[2] = 3;

      int baseLevel   = 0;
      int refineLevel = 2;

      LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter());

      //////////////////////////////////////////////////////////////////////////
      //physik
      //////////////////////////////////////////////////////////////////////////
      LBMReal uLB = 0.05;
      LBMReal rhoLB = 0.0;
      LBMReal nuLB = 1e-5;

      Grid3DPtr grid(new Grid3D(comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(1000,1000,10000000));
      RestartPostprocessor rp(grid, rSch, comm, pathname, RestartPostprocessor::BINARY);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep() == 0)
      {

         if(myid==0) UBLOG(logINFO,"Neustart..");

         //////////////////////////////////////////////////////////////////////////
         //town
         GbTriFaceMesh3DPtr town(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(townFilename, "Netz"));
         if(myid == 0) GbSystem3D::writeGeoObject( town.get(), pathname+"/geo/town", WbWriterVtkXmlBinary::getInstance() );
         //////////////////////////////////////////////////////////////////////////

         //double cdx = 0.8;
         double cdx = town->getX3Maximum() / (double)(nx[2] * blocknx[2]);
         double fdx = cdx/double(1<<refineLevel);
 
         double blockLengthx = blocknx[0]*cdx; //geowerte

         double geoLength[] = { nx[0] * blockLengthx, nx[1] * blockLengthx, nx[2] * blockLengthx };

         double originX1 = town->getX1Minimum();
         double originX2 = town->getX2Minimum();
         double originX3 = town->getX3Minimum();


         bool periodicx1 = true;
         bool periodicx2 = true;
         bool periodicx3 = false;

         //bounding box
         double g_minX1 = originX1-3.0*blockLengthx;
         double g_minX2 = originX2-3.0*blockLengthx;
         double g_minX3 = originX3;

         double g_maxX1 = originX1 + geoLength[0]+3.0*blockLengthx;
         double g_maxX2 = originX2 + geoLength[1]+1.0*blockLengthx;
         double g_maxX3 = originX3 + geoLength[2]+2.0*blockLengthx;

         //double g_maxX1 = town->getX1Maximum()+blockLengthx;
         //double g_maxX2 = town->getX2Maximum()+2.0*blockLengthx;
         //double g_maxX3 = town->getX3Maximum()+2.0*blockLengthx;


         //set grid
         grid->setDeltaX(cdx);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
         grid->setPeriodicX1(periodicx1);
         grid->setPeriodicX2(periodicx2);
         grid->setPeriodicX3(periodicx3);

         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if(myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);


         //////////////////////////////////////////////////////////////////////////
         if(myid == 0)
         {
            UBLOG(logINFO, "*****************************************");
            UBLOG(logINFO, "* Parameters                            *");
            //UBLOG(logINFO, "* Re            ="<<Re);
            UBLOG(logINFO, "* nuLB          ="<<nuLB);
            UBLOG(logINFO, "* uLB           ="<<uLB);
            UBLOG(logINFO, "* cdx           ="<<cdx);
            UBLOG(logINFO, "* fdx           ="<<fdx);
            UBLOG(logINFO, "* blocknx1/2/3  ="<<blocknx[0]<<"/"<<blocknx[1]<<"/"<<blocknx[2]);
            UBLOG(logINFO, "* x1Periodic    ="<<periodicx1);
            UBLOG(logINFO, "* x2Periodic    ="<<periodicx2);
            UBLOG(logINFO, "* x3Periodic    ="<<periodicx3);
            UBLOG(logINFO, "* number of levels  ="<<refineLevel+1);
            UBLOG(logINFO, "* path          ="<<pathname);

            UBLOG(logINFO, "*****************************************");
            UBLOG(logINFO, "* number of threads    ="<<numOfThreads);
            UBLOG(logINFO, "* number of processes  ="<<comm->getNumberOfProcesses());
            UBLOG(logINFO, "*****************************************");
            UBLOG(logINFO, "*****************************************");     
         }
         //////////////////////////////////////////////////////////////////////////


         //////////////////////////////////////////////////////////////////////////
         //refinement

         /////////////////////////////////////////////////
         ///interactor
         int bbOption1 = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcNoSlip(new D3Q27NoSlipBCAdapter(bbOption1));
         D3Q27TriFaceMeshInteractorPtr triTownInteractor(new D3Q27TriFaceMeshInteractor(town, grid, bcNoSlip, Interactor3D::SOLID));

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLengthx, g_minX2-blockLengthx, g_maxX3, g_maxX1+blockLengthx, g_maxX2+blockLengthx, g_maxX3+blockLengthx));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr velBCInteractor(new D3Q27Interactor(addWallZmax, grid, Interactor3D::SOLID));

         double raiseVelSteps = 0;
         vector<D3Q27BCFunction> velcX2BCs, dummy;

         mu::Parser inflowProfile;
         inflowProfile.SetExpr("uLB");
         inflowProfile.DefineConst("uLB", uLB);
         velcX2BCs.push_back(D3Q27BCFunction(inflowProfile, raiseVelSteps, D3Q27BCFunction::INFCONST));

         D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter(dummy, velcX2BCs, dummy));
         velBCInteractor->addBCAdapter(velBCAdapter);

         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLengthx, g_minX2-blockLengthx, g_minX3-blockLengthx, g_maxX1+blockLengthx, g_maxX2+blockLengthx, g_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcNoSlip, Interactor3D::SOLID));

         //GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLengthx, g_minX2-blockLengthx, g_maxX3, g_maxX1+blockLengthx, g_maxX2+blockLengthx, g_maxX3+blockLengthx));
         //if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());
         //D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcNoSlip, Interactor3D::SOLID));

         GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1-blockLengthx, g_minX2-blockLengthx, g_minX3-blockLengthx, g_maxX1+blockLengthx, g_minX2, g_maxX3+blockLengthx));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcNoSlip, Interactor3D::SOLID));

         GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1-blockLengthx, g_minX2-blockLengthx, g_maxX3, g_maxX1+blockLengthx, g_maxX2+blockLengthx, g_maxX3+blockLengthx));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcNoSlip, Interactor3D::SOLID));

         GbCuboid3DPtr addWallXmin(new GbCuboid3D(g_minX1-blockLengthx, g_minX2-blockLengthx, g_minX3-blockLengthx, g_minX1, g_maxX2+blockLengthx, g_maxX3+blockLengthx));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallXmin.get(), pathname+"/geo/addWallXmin", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr addWallXminInt(new D3Q27Interactor(addWallXmin, grid, bcNoSlip, Interactor3D::SOLID));

         GbCuboid3DPtr addWallXmax(new GbCuboid3D(g_maxX1, g_minX2-blockLengthx, g_minX3-blockLengthx, g_maxX1+blockLengthx, g_maxX2+blockLengthx, g_maxX3+blockLengthx));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallXmax.get(), pathname+"/geo/addWallXmax", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr addWallXmaxInt(new D3Q27Interactor(addWallXmax, grid, bcNoSlip, Interactor3D::SOLID));

         GbCuboid3DPtr refineTownBox(new GbCuboid3D(town->getX1Minimum(), town->getX2Minimum(), town->getX3Minimum(), town->getX1Maximum(), town->getX2Maximum(), town->getX3Maximum()));
         if (myid == 0) GbSystem3D::writeGeoObject(refineTownBox.get(), pathname + "/geo/refineTownBox", WbWriterVtkXmlASCII::getInstance());

         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
            //RefineAroundGbObjectHelper refineHelper(grid, refineLevel, boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triTownInteractor), 0.0, 3.0, comm);
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            //refineHelper.addGbObject(refineTownBox, refineLevel);
            refineHelper.addGbObject(town, refineLevel);
            refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         }

         //Grid3D::BlockIDMap bmap = grid->getBlockIDs();
         //bmap.clear();
         //(grid->getBlockIDs()).clear();
         //grid->deleteBlockIDs();

         //RenumberBlockVisitor renumber;
         //grid->accept(renumber);


         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));	

         ////////////////////////////////////////////
         /////delete solid blocks
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         //intHelper.addInteractor(triTownInteractor);
         intHelper.addInteractor(velBCInteractor);
         intHelper.addInteractor(addWallZminInt);
         //intHelper.addInteractor(addWallZmaxInt);
         //intHelper.addInteractor(addWallYminInt);
         //intHelper.addInteractor(addWallYmaxInt);
         //intHelper.addInteractor(addWallXminInt);
         //intHelper.addInteractor(addWallXmaxInt);
         intHelper.selectBlocks();
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	 
         //////////////////////////////////////

         //grid->accept(renumber);

         //if (myid == 0)
         {
            UBLOG(logINFO, "Write blocks - start");
            BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
            ppblocks->update(0);
            UBLOG(logINFO, "Write blocks - end");
         }

         


         //domain decomposition for threads
         if(numOfThreads > 1)
         {
            PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
            grid->accept(pqPartVisitor);
         }


         unsigned long nob = grid->getNumberOfBlocks();
         unsigned long nod = nob * blocknx[0]*blocknx[1]*blocknx[2];
         unsigned long nod_real = nob * (blocknx[0]+3)*(blocknx[1]+3)*(blocknx[2]+3);
         unsigned long nodb = (blocknx[0]) * (blocknx[1]) * (blocknx[2]);

         double needMemAll  = double(nod_real*(27*sizeof(double) + sizeof(int)));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());
         
         double nup = 0; 

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
               nup += nobl*nodb*double(1<<level); 
            }
            UBLOG(logINFO,"Hypothetically time for calculation step for 120 nodes  = " << nup/6.0e5/(120*8)  << " s");
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
            UBLOG(logINFO,"Available memory per node/8.0 = " << (availMem/8.0) << " bytes");
         }
         //////////////////////////////////////////
         //set connectors
         if(myid == 0) UBLOG(logINFO,"set connectors - start");
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept( setConnsVisitor );
         if(myid == 0) UBLOG(logINFO,"set connectors - end");

         ////////////////////////////
         LBMKernel3DPtr kernel;
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));

         //mu::Parser fctForcingX2;
         //fctForcingX2.SetExpr("Fx2*dx");
         //fctForcingX2.DefineConst("Fx2", 5e-6);

         //kernel->setForcingX2(fctForcingX2);
         //kernel->setWithForcing(true);

         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         //BCProcessorPtr bcProc(new D3Q27ETForThinWallBCProcessor());
         kernel->setBCProcessor(bcProc);
         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);
         //////////////////////////////////
         //undef nodes
         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }


         intHelper.setBC();

         //initialization of decompositions
         D3Q27ETInitDistributionsBlockVisitor initVisitor( nuLB,rhoLB);
         initVisitor.setVx2(uLB);
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
         //domain decomposition for threads
         if(numOfThreads > 1)
         {
            PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
            grid->accept(pqPartVisitor);
         }
         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept( setConnsVisitor );

         if(myid == 0) UBLOG(logINFO,"Restart - end"); 
      }
      UbSchedulerPtr visSch(new UbScheduler());
      visSch->addSchedule(1,0,3);
      //visSch->addSchedule(100,100,1000);
      //visSch->addSchedule(1000,1000,5000);
      //visSch->addSchedule(5000,5000,100000);
      //visSch->addSchedule(100000,100000,10000000);

      visSch->addSchedule(1000,1000,10000000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), unitConverter);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 10, 30));
      nupsSch->addSchedule(500,500,1e6);
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      //UbSchedulerPtr emSch(new UbScheduler(100));
      //EmergencyExitPostprocessor empr(grid, emSch, pathname, RestartPostprocessorPtr(&rp), comm);

      if(myid == 0)
      {
         UBLOG(logINFO,"PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      string lastStep = "1000000";// string(cstr2);
      double endTime = UbSystem::stringTo<double>(lastStep);
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
//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   if (argc == 1)
   {
      cout<<"Command line argument isn't specified!"<<endl;
      cout<<"plate2 <machine name>"<<endl;
      return 1;
   }
   run(argv[1], argv[2]);

   return 0;
}


#include <iostream>
#include <string>
#include "VirtualFluids.h"



//#include <thread>

using namespace std;

//////////////////////////////////////////////////////////////////////////
void run(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathname          = config.getString("pathname");
      string          pathGeo           = config.getString("pathGeo");
      int             numOfThreads      = config.getInt("numOfThreads");
      string          sampleFilename    = config.getString("sampleFilename");
      vector<int>     pmNX              = config.getVector<int>("pmNX");
      double          lthreshold        = config.getDouble("lthreshold");
      double          uthreshold        = config.getDouble("uthreshold");
      vector<float>   voxelDeltaX       = config.getVector<float>("voxelDeltaX");
      vector<int>     blocknx           = config.getVector<int>("blocknx");
      double          u_LB              = config.getDouble("u_LB");
      double          restartStep       = config.getDouble("restartStep");
      double          restartStepStart  = config.getDouble("restartStepStart");
      double          endTime           = config.getDouble("endTime");
      double          outTime           = config.getDouble("outTime");
      double          availMem          = config.getDouble("availMem");
      bool            rawFile           = config.getBool("rawFile");
      bool            logToFile         = config.getBool("logToFile");
      bool            writeSample       = config.getBool("writeSample");
      vector<double>  pmL               = config.getVector<double>("pmL");
      double          deltaXfine        = config.getDouble("deltaXfine");
      int             refineLevel       = config.getInt("refineLevel");
      bool            thinWall          = config.getBool("thinWall");
      double          Re                = config.getDouble("Re");
      double          channelHigh       = config.getDouble("channelHigh");
      double          lengthFactor      = config.getDouble("lengthFactor");
      double          forcing           = config.getDouble("forcing");
      bool            changeQs          = config.getBool("changeQs"); 
      double          timeAvStart       = config.getDouble("timeAvStart");
      double          timeAvStop        = config.getDouble("timeAvStop");


      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (logToFile)
      {
#if defined(__unix__)
         if (myid == 0)
         {
            const char* str = pathname.c_str();
            mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         }
#endif 

         if (myid == 0)
         {
            stringstream logFilename;
            logFilename << pathname + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }

      //Sleep(30000);

      if (myid == 0) UBLOG(logINFO, "Testcase porous channel");

      LBMReal rho_LB = 0.0;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      double deltaXcoarse = deltaXfine*(double)(1<<refineLevel);

      double coord[6];
      bool restart;

      vector<double> origin(3);
      origin[0] = 0;
      origin[1] = 0;
      origin[2] = 0;

      //real velocity is 49.63 m/s
      double nu_LB;

      Grid3DPtr grid(new Grid3D(comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(restartStep, restartStepStart));
      RestartCoProcessor rp(grid, rSch, comm, pathname, RestartCoProcessor::TXT);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep() == 0)
      {
         if (myid == 0) UBLOG(logINFO, "new start...");
         restart = false;

         double offsetMinX3 = pmL[2];
         
         double offsetMaxX1 = pmL[0]*lengthFactor;
         double offsetMaxX2 = pmL[1]*2.0;
         double offsetMaxX3 = pmL[2] + channelHigh; //DLR-F15  //pmL[2]*2.0;

         //bounding box
         double g_minX1 = origin[0];
         double g_minX2 = origin[1];
         double g_minX3 = origin[2];

         double g_maxX1 = origin[0] + offsetMaxX1;
         double g_maxX2 = origin[1] + offsetMaxX2;
         double g_maxX3 = origin[2] + offsetMaxX3;
//////////////////////////////////////////////////////////////////////////
         double nx1_temp = floor((g_maxX1-g_minX1) /(deltaXcoarse*(double)blocknx[0]));

         deltaXcoarse = (g_maxX1-g_minX1) /(nx1_temp*(double)blocknx[0]);

         g_maxX1 -= 0.5* deltaXcoarse;
//////////////////////////////////////////////////////////////////////////
         double blockLength = (double)blocknx[0]*deltaXcoarse;

         grid->setPeriodicX1(true);
         grid->setPeriodicX2(true);
         grid->setPeriodicX3(false);
         grid->setDeltaX(deltaXcoarse);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);
         double channel_high = channelHigh; // g_maxX3-g_minX3;
         double channel_high_LB = channel_high/deltaXcoarse;
//////////////////////////////////////////////////////////////////////////
         //nu_LB = 0.005;
         nu_LB = (u_LB*channel_high_LB)/Re;
//////////////////////////////////////////////////////////////////////////
         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "u_LB = " << u_LB);
            UBLOG(logINFO, "rho_LB = " << rho_LB);
            UBLOG(logINFO, "nu_LB = " << nu_LB);
            UBLOG(logINFO, "dx coarse = " << deltaXcoarse << " m");
            UBLOG(logINFO, "dx fine = " << grid->getDeltaX(refineLevel) << " m");
            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
            UBLOG(logINFO, "number of threads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
            UBLOG(logINFO, "Preprocess - start");
         }

         //////////////////////////////////////////////////////////////////////////
         //refinement
         double blockLengthX3Fine = grid->getDeltaX(refineLevel) * blocknx[2];

         GbCuboid3DPtr refineBoxTop(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3-blockLengthX3Fine, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxTop.get(), pathname + "/geo/refineBoxTop", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+offsetMinX3+blockLengthX3Fine));
         //GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLengthX3Fine, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+blockLengthX3Fine));
         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxBottom.get(), pathname + "/geo/refineBoxBottom", WbWriterVtkXmlASCII::getInstance());

         if (refineLevel > 0)
         {
            if (myid == 0) UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            refineHelper.addGbObject(refineBoxTop, refineLevel);
            refineHelper.addGbObject(refineBoxBottom, refineLevel);
            refineHelper.refine();
            if (myid == 0) UBLOG(logINFO, "Refinement - end");
         }
         //////////////////////////////////////////////////////////////////////////

         //walls
         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());


         //wall interactors
         int bbOption = 1;
         D3Q27BoundaryConditionAdapterPtr bcNoSlip(new D3Q27NoSlipBCAdapter(bbOption));
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcNoSlip, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcNoSlip, Interactor3D::SOLID));

		 ////////////////////////////////////////////////
		 //TEST
		 //GbObject3DPtr testCube(new GbCuboid3D(g_minX1 + 2.0 * blockLength, g_minX2 + 2.0 * blockLength, g_minX3 + 5.0 * blockLength, 
			// g_minX1 + 3.0 * blockLength, g_minX2 + 3.0 * blockLength, g_minX3 + 6.0 * blockLength));
		 //if (myid == 0) GbSystem3D::writeGeoObject(testCube.get(), pathname + "/geo/testCube", WbWriterVtkXmlBinary::getInstance());
		 //D3Q27InteractorPtr testCubeInt(new D3Q27Interactor(testCube, grid, bcNoSlip, Interactor3D::SOLID));
		 ///////////////////////////////////////////////

         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         //Zoltan
         //Grid3DVisitorPtr zoltanVisitor(new ZoltanPartitioningGridVisitor(comm, D3Q27System::BSW, 1));
         //grid->accept(zoltanVisitor);
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
		 //////////////////////////////////////////////////////////////////////////
		 //TEST
		 //intHelper.addInteractor(testCubeInt);
         //////////////////////////////////////////////////////////////////////////
		 intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
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

         LBMKernel3DPtr kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));

         mu::Parser fctForcingX1;
         fctForcingX1.SetExpr("Fx1");
         fctForcingX1.DefineConst("Fx1", 1.0e-6);

         kernel->setWithForcing(true);

         BCProcessorPtr bcProc;
         BoundaryConditionPtr noSlipBC;

         if (thinWall)
         {
            bcProc = BCProcessorPtr(new D3Q27ETForThinWallBCProcessor());
            noSlipBC = BoundaryConditionPtr(new ThinWallNoSlipBoundaryCondition());
         }
         else
         {
            bcProc = BCProcessorPtr(new D3Q27ETBCProcessor());
            noSlipBC = BoundaryConditionPtr(new NoSlipBoundaryCondition());
         }

         bcProc->addBC(noSlipBC);

         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nu_LB, availMem, needMem);
         grid->accept(kernelVisitor);

         //////////////////////////////////
         //undef nodes for refinement
         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }


         //BC
         intHelper.setBC();

         ////porous media
         {
            string samplePathname = pathGeo + "/" + sampleFilename;

            GbVoxelMatrix3DPtr sample(new GbVoxelMatrix3D(pmNX[0], pmNX[1], pmNX[2], 0, lthreshold, uthreshold));
            if (rawFile)
            {
               sample->readBufferedMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
            }
            else
            {
               sample->readMatrixFromVtiASCIIFile(samplePathname);
            }

            sample->setVoxelMatrixDelta(voxelDeltaX[0], voxelDeltaX[1], voxelDeltaX[2]);
            sample->setVoxelMatrixMininum(origin[0], origin[1], origin[2]);

            int bounceBackOption = 1;
            bool vxFile = false;
            int i = 0;
            for (int x = 0; x < lengthFactor; x+=2)
            {
               double offset = pmL[0] * (double)x;
               //sample 0
               if (myid == 0) UBLOG(logINFO, "sample # " << i);
               sample->setVoxelMatrixMininum(origin[0]+offset, origin[1], origin[2]);
               Utilities::voxelMatrixDiscretisation(sample, pathname, myid, i, grid, bounceBackOption, vxFile);
               i++;

               if (myid == 0)
               {
                  UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
                  UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
                  UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
               }

               //sample 1
               if (myid == 0) UBLOG(logINFO, "sample # " << i);
               sample->setVoxelMatrixMininum(origin[0]+pmL[0]+offset, origin[1], origin[2]);
               sample->mirrorX();
               Utilities::voxelMatrixDiscretisation(sample, pathname, myid, i, grid, bounceBackOption, vxFile);
               i++;

               if (myid == 0)
               {
                  UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
                  UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
                  UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
               }

               //sample 2
               if (myid == 0) UBLOG(logINFO, "sample # " << i);
               sample->setVoxelMatrixMininum(origin[0]+pmL[0]+offset, origin[1]+pmL[1], origin[2]);
               sample->mirrorY();
               Utilities::voxelMatrixDiscretisation(sample, pathname, myid, i, grid, bounceBackOption, vxFile);
               i++;

               if (myid == 0)
               {
                  UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
                  UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
                  UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
               }

               //sample 3
               if (myid == 0) UBLOG(logINFO, "sample # " << i);
               sample->setVoxelMatrixMininum(origin[0]+offset, origin[1]+pmL[1], origin[2]);
               sample->mirrorX();
               Utilities::voxelMatrixDiscretisation(sample, pathname, myid, i, grid, bounceBackOption, vxFile);
               sample->mirrorY();
               i++;
            }

         }
         BoundaryConditionBlockVisitor bcVisitor;
         grid->accept(bcVisitor);

         mu::Parser inflowProfile;
         inflowProfile.SetExpr("x3 < h ? 0.0 : uLB+1*x1-1*x2");
		   //inflowProfile.SetExpr("uLB+1*x1-1*x2");
         //inflowProfile.SetExpr("uLB");
         inflowProfile.DefineConst("uLB", u_LB);
         inflowProfile.DefineConst("h", pmL[2]);

         D3Q27ETInitDistributionsBlockVisitor initVisitor(nu_LB, rho_LB);
         initVisitor.setVx1(inflowProfile);
         //initVisitor.setVx1(u_LB);
         grid->accept(initVisitor);

         ////set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         //ConnectorFactoryPtr factory(new Block3DConnectorFactory());
         //ConnectorBlockVisitor setConnsVisitor(comm, nu_LB, iProcessor, factory);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         MacroscopicQuantitiesCoProcessorPtr ppgeo(
            new MacroscopicQuantitiesCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
         ppgeo->process(0);
         ppgeo.reset();

         coord[0] = g_minX1;
         coord[1] = g_minX2;
         coord[2] = g_minX3;
         coord[3] = g_maxX1;
         coord[4] = g_maxX2;
         coord[5] = g_maxX3;

         ////////////////////////////////////////////////////////
         FILE * pFile;
         string str = pathname + "/checkpoints/coord.txt";
         pFile = fopen(str.c_str(), "w");
         fprintf(pFile, "%g\n", deltaXcoarse);
         fprintf(pFile, "%g\n", nu_LB);
         fprintf(pFile, "%g\n", coord[0]);
         fprintf(pFile, "%g\n", coord[1]);
         fprintf(pFile, "%g\n", coord[2]);
         fprintf(pFile, "%g\n", coord[3]);
         fprintf(pFile, "%g\n", coord[4]);
         fprintf(pFile, "%g\n", coord[5]);
         fclose(pFile);
         ////////////////////////////////////////////////////////

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         ////////////////////////////////////////////////////////
         FILE * pFile;
         string str = pathname + "/checkpoints/coord.txt";
         pFile = fopen(str.c_str(), "r");
         fscanf(pFile, "%lg\n", &deltaXcoarse);
         fscanf(pFile, "%lg\n", &nu_LB);
         fscanf(pFile, "%lg\n", &coord[0]);
         fscanf(pFile, "%lg\n", &coord[1]);
         fscanf(pFile, "%lg\n", &coord[2]);
         fscanf(pFile, "%lg\n", &coord[3]);
         fscanf(pFile, "%lg\n", &coord[4]);
         fscanf(pFile, "%lg\n", &coord[5]);
         fclose(pFile);
         ////////////////////////////////////////////////////////

         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "u_LB = " << u_LB);
            UBLOG(logINFO, "rho_LB = " << rho_LB);
            UBLOG(logINFO, "nu_LB = " << nu_LB);
            UBLOG(logINFO, "dx coarse = " << deltaXcoarse << " m");
            UBLOG(logINFO, "dx fine = " << grid->getDeltaX(refineLevel) << " m");
            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
         }

         BoundaryConditionBlockVisitor bcVisitor;
         grid->accept(bcVisitor);

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept(setConnsVisitor);

         restart = true;

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      UbSchedulerPtr stepSch(new UbScheduler(outTime));

      MacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv);

      double startStep = grid->getTimeStep();

      //UbSchedulerPtr visSch(new UbScheduler());
      //visSch->addSchedule(40000,40000,40000000);
      //UbSchedulerPtr resSchRMS(new UbScheduler());
      //resSchRMS->addSchedule(40000, startStep, 40000000);
      //UbSchedulerPtr resSchMeans(new UbScheduler());
      //resSchMeans->addSchedule(40000, startStep, 40000000);
      //UbSchedulerPtr stepAvSch(new UbScheduler());
      //stepAvSch->addSchedule(100, 0, 10000000);
      //AverageValuesCoProcessor Avpp(grid, pathname, WbWriterVtkXmlBinary::getInstance(),
      //   stepSch/*wann wird rausgeschrieben*/, stepAvSch/*wann wird gemittelt*/, resSchMeans, resSchRMS/*wann wird resettet*/, restart);


      UbSchedulerPtr AdjForcSch(new UbScheduler());
      AdjForcSch->addSchedule(100, 0, 10000000);
      D3Q27IntegrateValuesHelperPtr intValHelp(new D3Q27IntegrateValuesHelper(grid, comm,
         coord[0], coord[1], coord[2],
         coord[3], coord[4], coord[5]));
      if (myid == 0) GbSystem3D::writeGeoObject(intValHelp->getBoundingBox().get(), pathname + "/geo/IntValHelp", WbWriterVtkXmlBinary::getInstance());

      double vxTarget=u_LB;
      AdjustForcingCoProcessor AdjForcPPPtr(grid, AdjForcSch, pathname, intValHelp, vxTarget, forcing, comm);

      //mu::Parser decrViscFunc;
      //decrViscFunc.SetExpr("nue0+c0/(t+1)/(t+1)");
      //decrViscFunc.DefineConst("nue0", nu_LB*4.0);
      //decrViscFunc.DefineConst("c0", 0.1);
      //UbSchedulerPtr DecrViscSch(new UbScheduler());
      //DecrViscSch->addSchedule(10, 0, 1000);
      //DecreaseViscosityCoProcessor decrViscPPPtr(grid, DecrViscSch, &decrViscFunc, comm);

	  //if (changeQs)
	  //{
		 // double z1 = pmL[2];
		 // D3Q27IntegrateValuesHelperPtr intValHelp2(new D3Q27IntegrateValuesHelper(grid, comm,
			//  coord[0], coord[1], z1 - deltaXfine,
			//  coord[3], coord[4], z1 + deltaXfine));
		 // if (myid == 0) GbSystem3D::writeGeoObject(intValHelp2->getBoundingBox().get(), pathname + "/geo/intValHelp2", WbWriterVtkXmlBinary::getInstance());
		 // Utilities::ChangeRandomQs(intValHelp2);
	  //}


      UbSchedulerPtr tavSch(new UbScheduler(1, timeAvStart, timeAvStop));
      TimeAveragedValuesCoProcessorPtr tav(new TimeAveragedValuesCoProcessor(grid, pathname, WbWriterVtkXmlBinary::getInstance(), tavSch,
         TimeAveragedValuesCoProcessor::Velocity | TimeAveragedValuesCoProcessor::Fluctuations | TimeAveragedValuesCoProcessor::Triplecorrelations));
      
      //UbSchedulerPtr catalystSch(new UbScheduler(1));
      //InSituCatalystCoProcessor catalyst(grid, catalystSch, "pchannel.py");

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
      calculation->setTimeAveragedValuesCoProcessor(tav);
      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
      if (myid == 0) UBLOG(logINFO, "Simulation-end");
   }
   catch (exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch (string& s)
   {
      cerr << s << endl;
   }
   catch (...)
   {
      cerr << "unknown exception" << endl;
   }

}
//////////////////////////////////////////////////////////////////////////
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

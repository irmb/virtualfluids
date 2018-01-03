#include <iostream>
#include <string>

#include <VirtualFluids.h>

using namespace std;


//void pflowForcing(string configname)
//{
//   try
//   {
//      ConfigurationFile   config;
//      config.load(configname);
//
//      string          pathname = config.getString("pathname");
//      int             numOfThreads = config.getInt("numOfThreads");
//      vector<int>     blocknx = config.getVector<int>("blocknx");
//      vector<double>  gridnx = config.getVector<double>("gridnx");
//      double          nuLB = config.getDouble("nuLB");
//      double          endTime = config.getDouble("endTime");
//      double          outTime = config.getDouble("outTime");
//      double          availMem = config.getDouble("availMem");
//      int             refineLevel = config.getInt("refineLevel");
//      bool            logToFile = config.getBool("logToFile");
//      double          restartStep = config.getDouble("restartStep");
//      double          forcing = config.getDouble("forcing");
//      bool            thinWall = config.getBool("thinWall");
//      double          deltax = config.getDouble("deltax");
//
//
//      CommunicatorPtr comm = MPICommunicator::getInstance();
//      int myid = comm->getProcessID();
//
//      if (logToFile)
//      {
//#if defined(__unix__)
//         if (myid == 0)
//         {
//            const char* str = pathname.c_str();
//            mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//         }
//#endif 
//
//         if (myid == 0)
//         {
//            stringstream logFilename;
//            logFilename << pathname + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
//            UbLog::output_policy::setStream(logFilename.str());
//         }
//      }
//
//      double dx = deltax;
//
//      const int blocknx1 = blocknx[0];
//      const int blocknx2 = blocknx[1];
//      const int blocknx3 = blocknx[2];
//
//      LBMReal rhoLB = 0.0;
//
//      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());
//
//      const int baseLevel = 0;
//
//      //bounding box
//      double g_minX1 = 0;
//      double g_minX2 = 0;
//      double g_minX3 = 0;
//
//      double g_maxX1 = gridnx[0];
//      double g_maxX2 = gridnx[1];
//      double g_maxX3 = gridnx[2];
//
//      double blockLength = blocknx1*dx;
//
//      Grid3DPtr grid(new Grid3D(comm));
//      grid->setPeriodicX1(true);
//      grid->setPeriodicX2(true);
//      grid->setPeriodicX3(false);
//      grid->setDeltaX(dx);
//      grid->setBlockNX(blocknx1, blocknx2, blocknx3);
//
//      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
//      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
//
//      //////////////////////////////////////////////////////////////////////////
//      //restart
//      UbSchedulerPtr rSch(new UbScheduler(restartStep, restartStep));
//      RestartCoProcessor rp(grid, rSch, comm, pathname, RestartCoProcessor::TXT);
//      //////////////////////////////////////////////////////////////////////////
//
//      if (grid->getTimeStep() == 0)
//      {
//         GenBlocksGridVisitor genBlocks(gridCube);
//         grid->accept(genBlocks);
//
//         if (myid == 0)
//         {
//            UBLOG(logINFO, "Parameters:");
//            UBLOG(logINFO, "forcing = " << forcing);
//            UBLOG(logINFO, "rho = " << rhoLB);
//            UBLOG(logINFO, "nu = " << nuLB);
//            UBLOG(logINFO, "dx = " << dx);
//            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
//            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
//            UBLOG(logINFO, "Preprozess - start");
//         }
//
//         //////////////////////////////////////////////////////////////////////////
//         //refinement
//         double blockLengthX3Fine = grid->getDeltaX(refineLevel) * blocknx[2];
//
//         GbCuboid3DPtr refineBoxTop(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_maxX3 - blockLengthX3Fine, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
//         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxTop.get(), pathname + "/geo/refineBoxTop", WbWriterVtkXmlASCII::getInstance());
//
//         //GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+offsetMinX3+blockLengthX3Fine));
//         GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLengthX3Fine, g_maxX1 + blockLength, g_maxX2 + blockLength, g_minX3 + blockLengthX3Fine));
//         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxBottom.get(), pathname + "/geo/refineBoxBottom", WbWriterVtkXmlASCII::getInstance());
//
//         if (refineLevel > 0)
//         {
//            if (myid == 0) UBLOG(logINFO, "Refinement - start");
//            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
//            refineHelper.addGbObject(refineBoxTop, refineLevel);
//            refineHelper.addGbObject(refineBoxBottom, refineLevel);
//            refineHelper.refine();
//            if (myid == 0) UBLOG(logINFO, "Refinement - end");
//         }
//         //////////////////////////////////////////////////////////////////////////
//
//         //walls
//         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_minX3));
//         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname + "/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());
//
//         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_maxX3, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
//         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname + "/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());
//
//         //wall interactors
//         int bbOption = 1;
//         D3Q27BoundaryConditionAdapterPtr bcNoSlip(new D3Q27NoSlipBCAdapter(bbOption));
//         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcNoSlip, Interactor3D::SOLID));
//         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcNoSlip, Interactor3D::SOLID));
//
//         ////////////////////////////////////////////
//         //METIS
//         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
//         ////////////////////////////////////////////
//         /////delete solid blocks
//         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
//         InteractorsHelper intHelper(grid, metisVisitor);
//         intHelper.addInteractor(addWallZminInt);
//         intHelper.addInteractor(addWallZmaxInt);
//         intHelper.selectBlocks();
//         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
//         //////////////////////////////////////
//
//         //set connectors
//         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
//         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
//         grid->accept(setConnsVisitor);
//
//         //domain decomposition for threads
//         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
//         grid->accept(pqPartVisitor);
//
//         WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
//         ppblocks->process(0);
//         ppblocks.reset();
//
//         unsigned long nob = grid->getNumberOfBlocks();
//         int gl = 3;
//         unsigned long nodb = (blocknx1) * (blocknx2) * (blocknx3);
//         unsigned long nod = nob * (blocknx1) * (blocknx2) * (blocknx3);
//         unsigned long nodg = nob * (blocknx1 + gl) * (blocknx2 + gl) * (blocknx3 + gl);
//         double needMemAll = double(nodg*(27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
//         double needMem = needMemAll / double(comm->getNumberOfProcesses());
//
//         if (myid == 0)
//         {
//            UBLOG(logINFO, "Number of blocks = " << nob);
//            UBLOG(logINFO, "Number of nodes  = " << nod);
//            int minInitLevel = grid->getCoarsestInitializedLevel();
//            int maxInitLevel = grid->getFinestInitializedLevel();
//            for (int level = minInitLevel; level <= maxInitLevel; level++)
//            {
//               int nobl = grid->getNumberOfBlocks(level);
//               UBLOG(logINFO, "Number of blocks for level " << level << " = " << nob);
//               UBLOG(logINFO, "Number of nodes for level " << level << " = " << nob*nodb);
//            }
//            UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
//            UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
//            UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
//         }
//
//         LBMKernel3DPtr kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));
//
//         mu::Parser fctForcingX1;
//         fctForcingX1.SetExpr("Fx1");
//         fctForcingX1.DefineConst("Fx1", forcing);
//
//         kernel->setWithForcing(true);
//         kernel->setForcingX1(fctForcingX1);
//
//         BCProcessorPtr bcProc;
//         BoundaryConditionPtr noSlipBC;
//
//         if (thinWall)
//         {
//            bcProc = BCProcessorPtr(new D3Q27ETForThinWallBCProcessor());
//            noSlipBC = BoundaryConditionPtr(new ThinWallNoSlipBoundaryCondition());
//         }
//         else
//         {
//            bcProc = BCProcessorPtr(new D3Q27ETBCProcessor());
//            noSlipBC = BoundaryConditionPtr(new NoSlipBoundaryCondition());
//         }
//
//         bcProc->addBC(noSlipBC);
//
//         kernel->setBCProcessor(bcProc);
//
//         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
//         grid->accept(kernelVisitor);
//
//         //////////////////////////////////
//         //undef nodes for refinement
//         if (refineLevel > 0)
//         {
//            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
//            grid->accept(undefNodesVisitor);
//         }
//
//         //BC
//         intHelper.setBC();
//         BoundaryConditionBlockVisitor bcVisitor;
//         grid->accept(bcVisitor);
//
//         //initialization of distributions
//         D3Q27ETInitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
//         grid->accept(initVisitor);
//
//         //Postrozess
//         UbSchedulerPtr geoSch(new UbScheduler(1));
//         MacroscopicQuantitiesCoProcessorPtr ppgeo(
//            new MacroscopicQuantitiesCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
//         ppgeo->process(0);
//         ppgeo.reset();
//
//         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
//      }
//      else
//      {
//         mu::Parser fctForcingX1;
//         mu::Parser fctForcingX2;
//         mu::Parser fctForcingX3;
//         fctForcingX1.SetExpr("Fx1");
//         fctForcingX1.DefineConst("Fx1", forcing);
//         fctForcingX2.SetExpr("0.0");
//         fctForcingX3.SetExpr("0.0");
//
//         SetForcingBlockVisitor forcingVisitor(fctForcingX1, fctForcingX2, fctForcingX3);
//         grid->accept(forcingVisitor);
//
//         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
//         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
//         grid->accept(setConnsVisitor);
//
//         //domain decomposition for threads
//         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
//         grid->accept(pqPartVisitor);
//      }
//
//      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
//      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);
//
//      UbSchedulerPtr stepSch(new UbScheduler(outTime));
//
//      MacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv);
//
//      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
//      if (myid == 0) UBLOG(logINFO, "Simulation-start");
//      calculation->calculate();
//      if (myid == 0) UBLOG(logINFO, "Simulation-end");
//   }
//   catch (std::exception& e)
//   {
//      cerr << e.what() << endl << flush;
//   }
//   catch (std::string& s)
//   {
//      cerr << s << endl;
//   }
//   catch (...)
//   {
//      cerr << "unknown exception" << endl;
//   }
//
//}
//////////////////////////////////////////////////////////////////////////
void pflowdp(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathname     = config.getString("pathname");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx      = config.getVector<int>("blocknx");
      vector<double>  boundingBox       = config.getVector<double>("boundingBox");
      double          nuLB         = config.getValue<double>("nuLB");
      double          endTime      = config.getValue<double>("endTime");
      double          outTime      = config.getValue<double>("outTime");
      double          availMem     = config.getValue<double>("availMem");
      int             refineLevel  = config.getValue<int>("refineLevel");
      bool            logToFile    = config.getValue<bool>("logToFile");
      double          restartStep  = config.getValue<double>("restartStep");
      double          dpLB         = config.getValue<double>("dpLB");
      bool            thinWall     = config.getValue<bool>("thinWall");
      double          deltax       = config.getValue<double>("deltax");
      double          cpStep       = config.getDouble("cpStep");
      double          cpStepStart  = config.getDouble("cpStepStart");
      bool            newStart     = config.getValue<bool>("newStart");

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      LBMReal rhoLB = 0.0;
      double rhoLBinflow = dpLB*3.0;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;

      //bounding box
      double g_minX1 = 0;
      double g_minX2 = 0;
      double g_minX3 = 0;
      
      double g_maxX1 = boundingBox[0];
      double g_maxX2 = boundingBox[1];
      double g_maxX3 = boundingBox[2];

      double blockLength = 3.0*deltax;

      double h = (g_maxX2) / 2.0;
      double dex = g_maxX1;
      double Umax = (1.0 / (4.0*nuLB))*(dpLB / dex)*(h*h);
      double Re = (4 * h*Umax) / (3 * nuLB);

      //bc
      LBMReal uLB = 0.01;
      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);
      BCAdapterPtr denBCAdapterInflow(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      //denBCAdapterInflow->setBcAlgorithm(BCAlgorithmPtr(new VelocityWithDensityBCAlgorithm()));
      denBCAdapterInflow->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithm()));

      //BCAdapterPtr denBCAdapterOutflow(new DensityBCAdapter(rhoLB));
      //denBCAdapterOutflow->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithm()));
      ////denBCAdapterOutflow->setBcAlgorithm(BCAlgorithmPtr(new NonEqDensityBCAlgorithm()));

      //BCAdapterPtr slipBCAdapter(new SlipBCAdapter());
      ////slipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingSlipBCAlgorithm()));
      //slipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new SlipBCAlgorithm()));
      //

      //BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
      //noSlipBCAdapter->setBcAlgorithm(NoSlipBCAlgorithmPtr(new NoSlipBCAlgorithm()));

      //BoundaryConditionsBlockVisitor bcVisitor;
      //bcVisitor.addBC(noSlipBCAdapter);
      //bcVisitor.addBC(slipBCAdapter);
      //bcVisitor.addBC(denBCAdapterInflow);
      //bcVisitor.addBC(denBCAdapterOutflow);



      BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(NoSlipBCAlgorithmPtr(new NoSlipBCAlgorithm()));

      //BCAdapterPtr denBCAdapterInflow(new DensityBCAdapter(rhoLBinflow));
      //denBCAdapterInflow->setBcAlgorithm(BCAlgorithmPtr(new NonEqDensityBCAlgorithm()));

      BCAdapterPtr denBCAdapterOutflow(new DensityBCAdapter(rhoLB));
      denBCAdapterOutflow->setBcAlgorithm(BCAlgorithmPtr(new NonEqDensityBCAlgorithm()));

      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(denBCAdapterInflow);
      bcVisitor.addBC(denBCAdapterOutflow);

      Grid3DPtr grid(new Grid3D(comm));
      grid->setPeriodicX1(false);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(false);
      grid->setDeltaX(deltax);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

      double k1 = 4;
      double k2 = 8;

      GbObject3DPtr refineCube1_1(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2 / k1 - 1.0, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube1_1.get(), pathname + "/geo/refineCube1_1", WbWriterVtkXmlBinary::getInstance());

      GbObject3DPtr refineCube1_2(new GbCuboid3D(g_minX1, g_maxX2 - g_maxX2 / k1 + 1.0, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube1_2.get(), pathname + "/geo/refineCube1_2", WbWriterVtkXmlBinary::getInstance());

      GbObject3DPtr refineCube2_1(new GbCuboid3D(g_minX1 + 2 * blockLength + 2 * deltax, g_minX2, g_minX3, g_maxX1 - 2 * blockLength - 2 * deltax, g_maxX2 / k2 - 1.0, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube2_1.get(), pathname + "/geo/refineCube2_1", WbWriterVtkXmlBinary::getInstance());

      GbObject3DPtr refineCube2_2(new GbCuboid3D(g_minX1 + 2 * blockLength + 2 * deltax, g_maxX2 - g_maxX2 / k2 + 1.0, g_minX3, g_maxX1 - 2 * blockLength - 2 * deltax, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube2_2.get(), pathname + "/geo/refineCube2_2", WbWriterVtkXmlBinary::getInstance());

      GbObject3DPtr refineCube2_3(new GbCuboid3D(g_minX1 + blockLength + 2 * deltax, g_minX2 + blockLength + 2 * deltax, g_minX3 + blockLength + 2 * deltax, g_maxX1 - blockLength - 2 * deltax, g_maxX2 - blockLength - 2 * deltax, g_maxX3 - blockLength - 2 * deltax));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube2_3.get(), pathname + "/geo/refineCube2_3", WbWriterVtkXmlBinary::getInstance());

      //////////////////////////////////////////////////////////////////////////
      //restart
      //UbSchedulerPtr rSch(new UbScheduler(restartStep));
      //RestartCoProcessor rp(grid, rSch, comm, pathname, RestartCoProcessor::TXT);

      UbSchedulerPtr rSch2(new UbScheduler(cpStep, cpStepStart));
      MPIIORestart1CoProcessor rcp(grid, rSch2, pathname, comm);

      LBMKernelPtr kernel;
      kernel = LBMKernelPtr(new IncompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], IncompressibleCumulantLBMKernel::NORMAL));

      BCProcessorPtr bcProc(new BCProcessor());
      //BCProcessorPtr bcProc = BCProcessorPtr(new ThinWallBCProcessor());
      kernel->setBCProcessor(bcProc);

      rcp.setLBMKernel(kernel);
      rcp.setBCProcessor(bcProc);
      //rcp.setChunk(1);
      //////////////////////////////////////////////////////////////////////////

      if (newStart)
      {
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "h = " << h);
            UBLOG(logINFO, "rho = " << rhoLB);
            UBLOG(logINFO, "nue = " << nuLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << deltax);
            UBLOG(logINFO, "dpLB = " << dpLB);
            UBLOG(logINFO, "Umax = " << Umax);
            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
            UBLOG(logINFO, "Preprozess - start");
         }

         //walls
         GbCuboid3DPtr addWallYmin (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_minX2, g_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax (new GbCuboid3D(g_minX1-blockLength, g_maxX2, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());


         //GbCuboid3DPtr addWallXmax(new GbCuboid3D(g_maxX1-4.0*deltax, g_maxX2, g_minX3 - 4.0*blockLength, g_maxX1 + 4.0*blockLength, g_maxX2 + 4.0*blockLength, g_maxX3 + 4.0*blockLength));
         //if (myid == 0) GbSystem3D::writeGeoObject(addWallXmax.get(), pathname+"/geo/addWallXmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_minX1, g_maxX2 + blockLength, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2+blockLength, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         //GbCuboid3DPtr geoOutflowSolid(new GbCuboid3D(g_maxX1-1.0*deltax, g_minX2 - 4.0*blockLength, g_minX3 - 4.0*blockLength, g_maxX1 + 4.0*blockLength, g_maxX2+4.0*blockLength, g_maxX3 + 4.0*blockLength));
         //if (myid == 0) GbSystem3D::writeGeoObject(geoOutflowSolid.get(), pathname + "/geo/geoOutflowSolid", WbWriterVtkXmlASCII::getInstance());

         ////inflow
         //GbCuboid3DPtr geoInflow (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_minX3-4.0*blockLength, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_minX3));
         //if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         ////outflow
         //GbCuboid3DPtr geoOutflow (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_maxX3, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_maxX3+4.0*blockLength));
         //if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         if (refineLevel > 0)
         {
            if (myid == 0) UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel,comm);
            //refineHelper.addGbObject(refineCube1_1, 1);
            //refineHelper.addGbObject(refineCube1_2, 1);
            //refineHelper.addGbObject(refineCube2_1, 2);
            //refineHelper.addGbObject(refineCube2_2, 2);
            refineHelper.addGbObject(refineCube2_3, refineLevel);
            refineHelper.refine();
            if (myid == 0) UBLOG(logINFO, "Refinement - end");
         }

         //walls
         D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, noSlipBCAdapter,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, noSlipBCAdapter,Interactor3D::SOLID));

         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));

         //D3Q27InteractorPtr addWallXmaxInt(new D3Q27Interactor(addWallXmax, grid, denBCAdapterOutflow,Interactor3D::SOLID));

         D3Q27InteractorPtr inflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoInflow, grid, denBCAdapterInflow, Interactor3D::SOLID));

         //outflow
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, denBCAdapterOutflow, Interactor3D::SOLID));
         //D3Q27InteractorPtr outflowSolidInt = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, noSlipBCAdapter, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);

         //intHelper.addInteractor(addWallYminInt);
         //intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);

         intHelper.addInteractor(inflowInt);
         
         intHelper.addInteractor(outflowInt);

         //die Geschwindigkeit Randbedingung soll Ausflüß überdecken !!!!!



         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         //set connectors
         //InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         ppblocks->process(0);
         ppblocks.reset();

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nod = nob * (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nodg = nob * (blocknx[0] + gl) * (blocknx[1] + gl) * (blocknx[1] + gl);
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

         //LBMKernelPtr kernel;
         //kernel = LBMKernelPtr(new IncompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], IncompressibleCumulantLBMKernel::NORMAL));
         //kernel = LBMKernelPtr(new CompressibleCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], CompressibleCumulantLBMKernel::NORMAL));
         //}
         //kernel->setWithForcing(true);
         //kernel->setForcingX1(2e-6);
         //BCProcessorPtr bcProc(new BCProcessor());
         //kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //walls
         intHelper.setBC();

         grid->accept(bcVisitor);
         

         //initialization of distributions
         //mu::Parser fct;
         //fct.SetExpr("-(1.0/(2.0*nu))*(dp/dx)*((x2-h)^2 - h^2)");
         //fct.DefineConst("dp", dpLB);
         //fct.DefineConst("dx", dex);
         //fct.DefineConst("h", h);
         //fct.DefineConst("nu", nuLB);

         mu::Parser fct;
         fct.SetExpr("(x1max-x1)/l*dp*3.0");
         fct.DefineConst("dp", dpLB);
         fct.DefineConst("x1max", g_maxX1);
         fct.DefineConst("l", g_maxX1-g_minX1);

         InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         //initVisitor.setVx1(fct);
         //initVisitor.setVx1(uLB);
         //initVisitor.setVx3(fct);
         //initVisitor.setRho(fct);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         WriteBoundaryConditionsCoProcessorPtr ppgeo(
            new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         ppgeo->process(0);
         ppgeo.reset();

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         //rcp.readBlocks(restartStep);
         //Grid3DVisitorPtr newMetisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         //grid->accept(newMetisVisitor);
         //rcp.readDataSet(restartStep);
         //rcp.readBoundaryConds(restartStep);

         rcp.restart((int)restartStep);

         grid->setTimeStep(restartStep);

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         grid->accept(bcVisitor);

         UbSchedulerPtr geoSch(new UbScheduler(1));
         WriteBoundaryConditionsCoProcessor ppgeo = WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm);
         ppgeo.process(1);

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      UbSchedulerPtr stepSch(new UbScheduler(outTime));

      WriteMacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm);

      UbSchedulerPtr AdjForcSch(new UbScheduler());
      AdjForcSch->addSchedule(10, 0, 10000000);
      IntegrateValuesHelperPtr intValHelp(new IntegrateValuesHelper(grid, comm,
         g_minX1, g_minX2, g_minX3,
         g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(intValHelp->getBoundingBox().get(), pathname + "/geo/IntValHelp", WbWriterVtkXmlBinary::getInstance());

      double vxTarget=uLB;
      AdjustForcingCoProcessor AdjForcPPPtr(grid, AdjForcSch, pathname, intValHelp, vxTarget, comm);

      const std::shared_ptr<ConcreteCalculatorFactory> calculatorFactory = std::make_shared<ConcreteCalculatorFactory>(stepSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, calculatorFactory, CalculatorType::HYBRID));
       //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch, CalculationManager::MPI));
      //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch, CalculationManager::PrePostBc));
      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
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
//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   if (argv != NULL)
   {
      if (argv[1] != NULL)
      {
         //pflowForcing(string(argv[1]));
         pflowdp(string(argv[1]));
      }
      else
      {
         cout << "Configuration file is missing!" << endl;
      }
   }

   return 0;
}

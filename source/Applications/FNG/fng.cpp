#include <iostream>
#include <string>

#include <boost/pointer_cast.hpp>

#include "VirtualFluids.h"

using namespace std;


void run(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathOut = config.getString("pathOut");
      string          pathGeo = config.getString("pathGeo");
      string          fngFileWhole = config.getString("fngFileWhole");
      string          fngFileTrailingEdge = config.getString("fngFileTrailingEdge");
      string          fngFileBodyPart = config.getString("fngFileBodyPart");
      int             numOfThreads = config.getInt("numOfThreads");
      vector<int>     blockNx = config.getVector<int>("blockNx");
      vector<int>     boundingBox = config.getVector<int>("boundingBox");
      double          uLB = config.getDouble("uLB");
      double          restartStep = config.getDouble("restartStep");
      double          restartStepStart = config.getDouble("restartStepStart");
      double          endTime = config.getDouble("endTime");
      double          outTime = config.getDouble("outTime");
      double          availMem = config.getDouble("availMem");
      int             refineLevel = config.getInt("refineLevel");
      bool            logToFile = config.getBool("logToFile");
      bool            porousTralingEdge = config.getBool("porousTralingEdge");
      double          deltaXfine = config.getDouble("deltaXfine");

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (logToFile)
      {
#if defined(__unix__)
         if (myid == 0)
         {
            const char* str = pathOut.c_str();
            mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         }
#endif 

         if (myid == 0)
         {
            stringstream logFilename;
            logFilename << pathOut + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }

      GbTriFaceMesh3DPtr fngMeshWhole;
      GbTriFaceMesh3DPtr fngMeshBodyPart;
      GbTriFaceMesh3DPtr fngMeshTrailingEdge;
      if (porousTralingEdge)
      {
         if (myid == 0) UBLOG(logINFO, "Read fngFileBodyPart:start");
         fngMeshBodyPart = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo + "/" + fngFileBodyPart, "fngMeshBody", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         if (myid == 0) UBLOG(logINFO, "Read fngFileBodyPart:end");
         if (myid == 0) GbSystem3D::writeGeoObject(fngMeshBodyPart.get(), pathOut + "/geo/fngMeshBody", WbWriterVtkXmlBinary::getInstance());

         if (myid == 0) UBLOG(logINFO, "Read fngFileTrailingEdge:start");
         fngMeshTrailingEdge = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo + "/" + fngFileTrailingEdge, "fngMeshTrailingEdge", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         if (myid == 0) UBLOG(logINFO, "Read fngFileTrailingEdge:end");
         if (myid == 0) GbSystem3D::writeGeoObject(fngMeshTrailingEdge.get(), pathOut + "/geo/fngMeshTrailingEdge", WbWriterVtkXmlBinary::getInstance());
      }
      else
      {
         if (myid == 0) UBLOG(logINFO, "Read fngFileWhole:start");
         fngMeshWhole = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo + "/" + fngFileWhole, "fngMeshWhole", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         if (myid == 0) UBLOG(logINFO, "Read fngFileWhole:end");
         if (myid == 0) GbSystem3D::writeGeoObject(fngMeshWhole.get(), pathOut + "/geo/fngMeshWhole", WbWriterVtkXmlBinary::getInstance());
      }
      
      double g_minX1 = boundingBox[0];
      double g_minX2 = boundingBox[2];
      double g_minX3 = boundingBox[4];

      double g_maxX1 = boundingBox[1];
      double g_maxX2 = boundingBox[3];
      double g_maxX3 = boundingBox[5];

      //////////////////////////////////////////////////////////////////////////
      double deltaXcoarse = deltaXfine*(double)(1 << refineLevel);
      double nx1_temp = floor((g_maxX1 - g_minX1) / (deltaXcoarse*(double)blockNx[0]));

      deltaXcoarse = (g_maxX1 - g_minX1) / (nx1_temp*(double)blockNx[0]);

      g_maxX1 -= 0.5* deltaXcoarse;
      //////////////////////////////////////////////////////////////////////////
      double blockLength = (double)blockNx[0] * deltaXcoarse;

      //##########################################################################
      //## physical parameters
      //##########################################################################
      double Re = 1e6;

      double rhoLB = 0.0;
      double rhoReal = 1.2041; //(kg/m3)
      double nueReal = 153.5e-7; //m^2/s

      double lReal = 3.0;//m
      double uReal = Re*nueReal / lReal;

      //##Machzahl:
      //#Ma     = uReal/csReal
      double Ma = 0.15;//Ma-Real!
      //double csReal = uReal / Ma;
      //double hLB = lReal / deltaXcoarse;

      //LBMUnitConverter unitConverter(lReal, csReal, rhoReal, hLB);

      //double u_LB = uReal   * unitConverter.getFactorVelocityWToLb();
      //double nu_LB = nueReal * unitConverter.getFactorViscosityWToLb();
      double l_LB = 0.3 / deltaXcoarse;
      double nuLB = (uLB*l_LB) / Re; //0.005;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;

      ////////////////////////////////////////////////////////////////////////
      //Grid
      //////////////////////////////////////////////////////////////////////////
      Grid3DPtr grid(new Grid3D(comm));
      grid->setDeltaX(deltaXcoarse);
      grid->setBlockNX(blockNx[0], blockNx[1], blockNx[2]);

      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      //gridCube->setCenterCoordinates(geo->getX1Centroid(), geo->getX2Centroid(), geo->getX3Centroid());
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut + "/geo/gridCube", WbWriterVtkXmlASCII::getInstance());
      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      grid->setPeriodicX2(true);
      grid->setPeriodicX3(true);

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(restartStep, restartStep));
      RestartCoProcessor rp(grid, rSch, comm, pathOut, RestartCoProcessor::TXT);
      //////////////////////////////////////////////////////////////////////////


      if (grid->getTimeStep() == 0)
      {
         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "* Re            =" << Re);
            UBLOG(logINFO, "* Ma            =" << Ma);
            UBLOG(logINFO, "* uReal         =" << uReal);
            UBLOG(logINFO, "* nuReal        =" << nueReal);
            UBLOG(logINFO, "* nuLB          =" << nuLB);
            UBLOG(logINFO, "* velocity      =" << uLB);
            UBLOG(logINFO, "* dx_base       =" << deltaXcoarse);
            UBLOG(logINFO, "* dx_refine     =" << deltaXfine);

            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads     = " << numOfThreads);
            UBLOG(logINFO, "Preprozess - start");
         }


         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_minX1, g_maxX2 + blockLength, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathOut + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr noSlipBCAdapter(new D3Q27NoSlipBCAdapter(bbOption));

         Interactor3DPtr fngIntrWhole;
         Interactor3DPtr fngIntrBodyPart;
         Interactor3DPtr fngIntrTrailingEdge;
         if (porousTralingEdge)
         {
            fngIntrBodyPart = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshBodyPart, grid, noSlipBCAdapter, Interactor3D::SOLID));
            fngIntrTrailingEdge = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshTrailingEdge, grid, noSlipBCAdapter, Interactor3D::SOLID));
         }
         else
         {
            fngIntrWhole = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshWhole, grid, noSlipBCAdapter, Interactor3D::SOLID));
         }

         if (refineLevel > 0)
         {
            if (myid == 0) UBLOG(logINFO, "Refinement - start");
            //RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            //refineHelper.addGbObject(geo, refineLevel);
            //refineHelper.refine();
            
            //RefineAroundGbObjectHelper refineHelper1(grid, refineLevel-1, boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(geoIntr1), 0.0, 10.0, comm);
            //refineHelper1.refine();
            //RefineAroundGbObjectHelper refineHelper2(grid, refineLevel, boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(geoIntr2), -1.0, 5.0, comm);
            //refineHelper2.refine();

            if (porousTralingEdge)
            {
               int rank = grid->getRank();
               grid->setRank(0);
               boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrBodyPart)->refineBlockGridToLevel(refineLevel - 1, 0.0, 10.0);
               grid->setRank(rank);
            }
            else
            {
               int rank = grid->getRank();
               grid->setRank(0);
               boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrWhole)->refineBlockGridToLevel(refineLevel - 1, 0.0, 10.0);
               grid->setRank(rank);
            }


            ////////////////////////////////////////////
            //METIS
            Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
            ////////////////////////////////////////////
            /////delete solid blocks
            if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
            InteractorsHelper intHelper(grid, metisVisitor);
            if (porousTralingEdge)
            {
               intHelper.addInteractor(fngIntrBodyPart);
            }
            else
            {
               intHelper.addInteractor(fngIntrWhole);
            }
            //////////////////////////////////////////////////////////////////////////
            intHelper.selectBlocks();
            if (porousTralingEdge)
            {
               fngIntrBodyPart->removeSolidBlocks();
               fngIntrBodyPart->removeTransBlocks();
            }
            else
            {
               fngIntrWhole->removeSolidBlocks();
               fngIntrWhole->removeTransBlocks();
            }

            if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
            //////////////////////////////////////

            if (porousTralingEdge)
            {
               //grid->setRank(0);
               //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrTrailingEdge)->refineBlockGridToLevel(refineLevel, -1.0, 10.0);
               //grid->setRank(rank);

               GbObject3DPtr trailingEdgeCube(new GbCuboid3D(fngMeshTrailingEdge->getX1Minimum()-blockLength, fngMeshTrailingEdge->getX2Minimum(), fngMeshTrailingEdge->getX3Minimum()-blockLength,
                  fngMeshTrailingEdge->getX1Maximum()+blockLength, fngMeshTrailingEdge->getX2Maximum(), fngMeshTrailingEdge->getX3Maximum()+blockLength));
               if (myid == 0) GbSystem3D::writeGeoObject(trailingEdgeCube.get(), pathOut + "/geo/trailingEdgeCube", WbWriterVtkXmlASCII::getInstance());

               RefineCrossAndInsideGbObjectBlockVisitor refVisitor(trailingEdgeCube, refineLevel);
               grid->accept(refVisitor);
            }

            RatioBlockVisitor ratioVisitor(refineLevel);
            grid->accept(ratioVisitor);

            RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
            grid->accept(ratioSmoothVisitor);

            OverlapBlockVisitor overlapVisitor(refineLevel, false);
            grid->accept(overlapVisitor);

            std::vector<int> dirs;
            for (int i = D3Q27System::E; i <= D3Q27System::TS; i++)
            {
               dirs.push_back(i);
            }
            SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
            grid->accept(interDirsVisitor);

            if (myid == 0) UBLOG(logINFO, "Refinement - end");
         }


         //geo->scale(1.0/scaleFactorX, 1.0, 1.0/scaleFactorX);
         //geo = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(geoFile,"geo"));
         //if (myid == 0) GbSystem3D::writeGeoObject(geo.get(), pathOut + "/geo/geo3", WbWriterVtkXmlASCII::getInstance());

         mu::Parser fct;
         fct.SetExpr("U");
         fct.DefineConst("U", uLB);

         //inflow
         D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter(true, false, false, fct, 0, D3Q27BCFunction::INFCONST));
         velBCAdapter->setSecondaryBcOption(2);
         D3Q27InteractorPtr inflowIntr = D3Q27InteractorPtr(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         denBCAdapter->setSecondaryBcOption(0);
         D3Q27InteractorPtr outflowIntr = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(inflowIntr);
         intHelper.addInteractor(outflowIntr);
         if (porousTralingEdge)
         {
            intHelper.addInteractor(fngIntrBodyPart);
         } 
         else
         {
            intHelper.addInteractor(fngIntrWhole);
         }
         
         //////////////////////////////////////////////////////////////////////////
         intHelper.selectBlocks();

         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm));
         ppblocks->process(0);
         ppblocks.reset();

         unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
         int ghostLayer = 3;
         unsigned long long numberOfNodesPerBlock = (unsigned long long)(blockNx[0])* (unsigned long long)(blockNx[1])* (unsigned long long)(blockNx[2]);
         unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
         unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blockNx[0] + ghostLayer) * (blockNx[1] + ghostLayer) * (blockNx[2] + ghostLayer);
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

         LBMKernel3DPtr kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blockNx[0], blockNx[1], blockNx[2], LBMKernelETD3Q27CCLB::NORMAL));

         BCProcessorPtr bcProc;
         BoundaryConditionPtr noSlipBC;
         BoundaryConditionPtr velBC = VelocityBoundaryConditionPtr(new VelocityBoundaryCondition());
         BoundaryConditionPtr denBC = NonEqDensityBoundaryConditionPtr(new NonEqDensityBoundaryCondition());

         //if (thinWall)
         //{
            bcProc = BCProcessorPtr(new D3Q27ETForThinWallBCProcessor());
            noSlipBC = BoundaryConditionPtr(new ThinWallNoSlipBoundaryCondition());
         //}
         //else
         //{
         //   bcProc = BCProcessorPtr(new D3Q27ETBCProcessor());
         //   noSlipBC = BoundaryConditionPtr(new NoSlipBoundaryCondition());
         //}

         bcProc->addBC(noSlipBC);
         bcProc->addBC(velBC);
         bcProc->addBC(denBC);

         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //BC
         intHelper.setBC();
         
         //////////////////////////////////////////////////////////////////////////
         ///porous media

         //////////////////////////////////////////////////////////////////////////

         BoundaryConditionBlockVisitor bcVisitor;
         grid->accept(bcVisitor);

         //initialization of distributions
         D3Q27ETInitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         initVisitor.setVx1(fct);
         initVisitor.setNu(nuLB);
         grid->accept(initVisitor);

         ////set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         //ConnectorFactoryPtr factory(new Block3DConnectorFactory());
         //ConnectorBlockVisitor setConnsVisitor(comm, nu_LB, iProcessor, factory);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         MacroscopicQuantitiesCoProcessorPtr ppgeo(
            new MacroscopicQuantitiesCoProcessor(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, true));
         ppgeo->process(0);
         ppgeo.reset();

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);
      }

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      UbSchedulerPtr stepSch(new UbScheduler(outTime));

      MacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv);

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
      //calculation->setTimeAveragedValuesCoProcessor(tav);
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
         cout << "Configuration file must be set!: " << argv[0] << " <config file>" << endl << std::flush;
      }
   }

   return 0;
}


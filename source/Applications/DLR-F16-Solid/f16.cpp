#include <iostream>
#include <string>

#include <PointerDefinitions.h>
#include "VirtualFluids.h"
#include <omp.h>
using namespace std;

void run(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathOut = config.getValue<string>("pathOut");
      string          pathGeo = config.getValue<string>("pathGeo");
      string          fngFileWhole = config.getValue<string>("fngFileWhole");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blockNx = config.getVector<int>("blockNx");
      vector<double>  boundingBox = config.getVector<double>("boundingBox");
      double          restartStep = config.getValue<double>("restartStep");
      double          cpStart = config.getValue<double>("cpStart");
      double          cpStep = config.getValue<double>("cpStep");
      int             endTime = config.getValue<int>("endTime");
      double          outTimeStep = config.getValue<double>("outTimeStep");
      double          outTimeStart = config.getValue<double>("outTimeStart");
      double          availMem = config.getValue<double>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      bool            logToFile = config.getValue<bool>("logToFile");
      double          deltaXfine = config.getValue<double>("deltaXfine");
      double          refineDistance = config.getValue<double>("refineDistance");
      double          startDistance = config.getValue<double>("startDistance");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");
      bool            newStart = config.getValue<bool>("newStart");
      bool            writeBlocks = config.getValue<bool>("writeBlocks");
      string          pathReInit = config.getValue<string>("pathReInit");
      int             stepReInit = config.getValue<int>("stepReInit");

      double          pcpStart = config.getValue<double>("pcpStart");
      double          pcpStop  = config.getValue<double>("pcpStop");

      double          timeAvStart       = config.getValue<double>("timeAvStart");
      double          timeAvStop        = config.getValue<double>("timeAvStop");

      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (logToFile)
      {
#if defined(__unix__)
         if (myid==0)
         {
            const char* str = pathOut.c_str();
            mkdir(str, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
         }
#endif 

         if (myid==0)
         {
            stringstream logFilename;
            logFilename<<pathOut+"/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 1");
         UBLOG(logINFO, "PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }


      //the geometry is in mm

      double g_minX1 = boundingBox[0];
      double g_minX2 = boundingBox[2];
      double g_minX3 = boundingBox[4];

      double g_maxX1 = boundingBox[1];
      double g_maxX2 = boundingBox[3];
      double g_maxX3 = boundingBox[5];

      //////////////////////////////////////////////////////////////////////////
      double deltaXcoarse = deltaXfine*(double)(1<<refineLevel);
      //////////////////////////////////////////////////////////////////////////
      double blockLength = (double)blockNx[0]*deltaXcoarse;

      //##########################################################################
      //## physical parameters
      //##########################################################################
      double Re = 1e6;

      double rhoLB = 0.0;
      double rhoReal = 1.2041; //(kg/m3)
      //double nueReal = 153.5e-7; //m^2/s
      double uReal = 55; //m/s
      double lReal = 0.3;//m
      //double uReal = Re*nueReal / lReal;
      double nuReal = (uReal*lReal)/Re; //m^2/s

      //##Machzahl:
      //#Ma     = uReal/csReal
      double Ma = 0.15;//Ma-Real!
      double csReal = uReal / Ma;
      double hLB = lReal / deltaXcoarse;

      LBMUnitConverter unitConverter(lReal, csReal, rhoReal, hLB);

      double uLB = uReal   * unitConverter.getFactorVelocityWToLb();
      double nuLB = nuReal * unitConverter.getFactorViscosityWToLb();
      double lLB = lReal*1000.0/deltaXcoarse;
      //double nuLB = (uLB*lLB)/Re; //0.005;
      //double nuLB = 0.005;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      ////////////////////////////////////////////////////////////////////////
      //Grid
      //////////////////////////////////////////////////////////////////////////
      SPtr<Grid3D> grid(new Grid3D(comm));

      //BC adapters
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      SPtr<BCAdapter> slipBCAdapter(new SlipBCAdapter());
      slipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SlipBCAlgorithm()));

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);
      SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

      //fct.SetExpr("U");
      //fct.DefineConst("U", 0.01);
      //SPtr<BCAdapter> velBCAdapterOut(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      //velBCAdapterOut->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));

      SPtr<BCAdapter> outflowBCAdapter(new DensityBCAdapter(rhoLB));
      outflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor; 
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(velBCAdapter);
      bcVisitor.addBC(outflowBCAdapter);

      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
      //dynamicPointerCast<CompressibleCumulantLBMKernel>(kernel)->setRelaxationParameter(CompressibleCumulantLBMKernel::NORMAL);
      
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      //dynamicPointerCast<CompressibleCumulant4thOrderViscosityLBMKernel>(kernel)->setBulkViscosity(nuLB*2.0e3);


      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new BCProcessor());
      kernel->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIORestartCoProcessor> restartCoProcessor(new MPIIORestartCoProcessor(grid, rSch, pathOut, comm));
      restartCoProcessor->setLBMKernel(kernel);
      restartCoProcessor->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 2");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }


      if (newStart)
      {
         ////////////////////////////////////////////////////////////////////////
         //define grid
         //////////////////////////////////////////////////////////////////////////
         grid->setDeltaX(deltaXcoarse);
         grid->setBlockNX(blockNx[0], blockNx[1], blockNx[2]);

         SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid==0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         grid->setPeriodicX1(false);
         grid->setPeriodicX2(true);
         grid->setPeriodicX3(false);

         if (myid==0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "* Re                  = "<<Re);
            UBLOG(logINFO, "* Ma                  = "<<Ma);
            UBLOG(logINFO, "* velocity (uReal)    = "<<uReal<<" m/s");
            UBLOG(logINFO, "* viscosity (nuReal)  = "<<nuReal<<" m^2/s");
            UBLOG(logINFO, "* chord length (lReal)= "<<lReal<<" m");
            UBLOG(logINFO, "* velocity LB (uLB)   = "<<uLB);
            UBLOG(logINFO, "* viscosity LB (nuLB) = "<<nuLB);
            UBLOG(logINFO, "* chord length (l_LB) = "<<lLB<<" dx_base");
            UBLOG(logINFO, "* dx_base             = "<<deltaXcoarse<<" m");
            UBLOG(logINFO, "* dx_refine           = "<<deltaXfine<<" m");
            UBLOG(logINFO, "* blocknx             = "<<blockNx[0]<<"x"<<blockNx[1]<<"x"<<blockNx[2]);
            UBLOG(logINFO, "* refineDistance      = "<<refineDistance);
            UBLOG(logINFO, "* number of levels    = "<<refineLevel+1);
            UBLOG(logINFO, "* number of threads   = "<<numOfThreads);
            UBLOG(logINFO, "* number of processes = "<<comm->getNumberOfProcesses());
            UBLOG(logINFO, "* path = "<<pathOut);
            UBLOG(logINFO, "Preprozess - start");
         }


         //SPtr<GbObject3D> fngMeshWhole(new GbCylinder3D(15.0, 0.0, 0.0, 15.0, 100.0, 0.0, 25.0));
         //GbSystem3D::writeGeoObject(fngMeshWhole.get(), pathOut + "/geo/fngMeshWholeCylinder", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbTriFaceMesh3D> fngMeshWhole;
         if (myid==0) UBLOG(logINFO, "Read fngFileWhole:start");
         fngMeshWhole = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileWhole, "fngMeshWhole", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         if (myid==0) UBLOG(logINFO, "Read fngFileWhole:end");
         fngMeshWhole->rotate(0.0, 0.5, 0.0);
         //fngMeshWhole->scale(1e-3,1e-3,1e-3);
         //fngMeshWhole->translate(-150.0,-50.0,-1.28);
         if (myid==0) GbSystem3D::writeGeoObject(fngMeshWhole.get(), pathOut+"/geo/fngMeshWhole", WbWriterVtkXmlBinary::getInstance());

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 3");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }
         //////////////////////////////////////////////////////////////////////////
         SPtr<Interactor3D> fngIntrWhole;
         //fngIntrWhole = SPtr<D3Q27Interactor>(new D3Q27Interactor(fngMeshWhole, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::POINTS));
         fngIntrWhole = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshWhole, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::POINTS));

         if (refineLevel>0 && myid==0 && writeBlocks)
         {
            if (myid==0) UBLOG(logINFO, "Refinement - start");
            int rank = grid->getRank();
            grid->setRank(0);

            int level;

            //level 1
            level = 1;
            SPtr<GbObject3D> refCylinderL2(new GbCylinder3D(0.015, 0.0, 0.0, 0.015, 0.1, 0.0, 0.040));
            GbSystem3D::writeGeoObject(refCylinderL2.get(), pathOut + "/geo/refCylinderL2", WbWriterVtkXmlBinary::getInstance());
            RefineCrossAndInsideGbObjectBlockVisitor refVisitorCylinderL2(refCylinderL2, level);
            grid->accept(refVisitorCylinderL2);

            //SPtr<GbObject3D> refBoxL2(new GbCuboid3D(0.015, 0.0, -0.04, 1.100, 0.1, 0.04));
            SPtr<GbObject3D> refBoxL2(new GbCuboid3D(0.015, 0.0, -0.04, 0.5, 0.1, 0.04));
            if (myid==0) GbSystem3D::writeGeoObject(refBoxL2.get(), pathOut+"/geo/refBoxL2", WbWriterVtkXmlASCII::getInstance());
            RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL2(refBoxL2, level);
            grid->accept(refVisitorBoxL2);

            ////level 2
            //level = 2;
            //SPtr<GbObject3D> refCylinderL3(new GbCylinder3D(15.0, 0.0, 0.0, 15.0, 100.0, 0.0, 30.0));
            //GbSystem3D::writeGeoObject(refCylinderL3.get(), pathOut + "/geo/refCylinderL3", WbWriterVtkXmlBinary::getInstance());
            //RefineCrossAndInsideGbObjectBlockVisitor refVisitorCylinderL3(refCylinderL3, level);
            //grid->accept(refVisitorCylinderL3);

            //SPtr<GbObject3D> refBoxL3(new GbCuboid3D(15.0, 0.0, -30.0, 700.0, 100.0, 30.0));
            //if (myid==0) GbSystem3D::writeGeoObject(refBoxL3.get(), pathOut+"/geo/refBoxL3", WbWriterVtkXmlASCII::getInstance());
            //RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL3(refBoxL3, level);
            //grid->accept(refVisitorBoxL3);

            ////level 3
            //level = 3;
            //SPtr<GbObject3D> refCylinderL4(new GbCylinder3D(15.0, 0.0, 0.0, 15.0, 100.0, 0.0, 25.0));
            //GbSystem3D::writeGeoObject(refCylinderL4.get(), pathOut + "/geo/refCylinderL4", WbWriterVtkXmlBinary::getInstance());
            //RefineCrossAndInsideGbObjectBlockVisitor refVisitorCylinderL4(refCylinderL4, level);
            //grid->accept(refVisitorCylinderL4);

            //SPtr<GbObject3D> refBoxL4(new GbCuboid3D(15.0, 0.0, -25.0, 400.0, 100.0, 25.0));
            //if (myid==0) GbSystem3D::writeGeoObject(refBoxL4.get(), pathOut+"/geo/refBoxL4", WbWriterVtkXmlASCII::getInstance());
            //RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL4(refBoxL4, level);
            //grid->accept(refVisitorBoxL4);

            ////level 4
            //level = 4;
            //SPtr<GbObject3D> refBoxL5(new GbCuboid3D(120.0, 0.0, -9.0, 320.0, 100.0, 18.0));
            //if (myid==0) GbSystem3D::writeGeoObject(refBoxL5.get(), pathOut+"/geo/refBoxL5", WbWriterVtkXmlASCII::getInstance());
            //RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL5(refBoxL5, level);
            //grid->accept(refVisitorBoxL5);

            //last level
            //dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole)->refineBlockGridToLevel(refineLevel, startDistance, refineDistance);


            grid->setRank(rank);

            {
               WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(0);
            }

            /////delete solid blocks
            if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");

            SetSolidBlockVisitor v(fngIntrWhole, BlockType::SOLID);
            grid->accept(v);
            std::vector<SPtr<Block3D>>& sb = fngIntrWhole->getSolidBlockSet();
            for(SPtr<Block3D> block : sb)
            {
               grid->deleteBlock(block);
            }
            fngIntrWhole->removeSolidBlocks();
            fngIntrWhole->removeBcBlocks();

            if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - end");
            //////////////////////////////////////


            RatioBlockVisitor ratioVisitor(refineLevel);
            CheckRatioBlockVisitor checkRatio(refineLevel);
            int count = 0;

            do {
               grid->accept(ratioVisitor);
               checkRatio.resetState();
               grid->accept(checkRatio);
               if (myid==0) UBLOG(logINFO, "count = "<<count++<<" state = "<<checkRatio.getState());
            } while (!checkRatio.getState());


            {
               WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(1);
            }

            OverlapBlockVisitor overlapVisitor(refineLevel, false);
            grid->accept(overlapVisitor);

            if (myid==0) UBLOG(logINFO, "Refinement - end");
         }
         grid->updateDistributedBlocks(comm);

         std::vector<int> dirs;
         for (int i = D3Q27System::E; i<=D3Q27System::TS; i++)
         {
            dirs.push_back(i);
         }
         SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
         grid->accept(interDirsVisitor);

         //walls
         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
         if (myid==0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathOut+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathOut+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //wall interactors
         SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, velBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, velBCAdapter, Interactor3D::SOLID));

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathOut+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         //inflow
         SPtr<D3Q27Interactor> inflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //outflow
         SPtr<D3Q27Interactor> outflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowBCAdapter, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(inflowIntr);
         intHelper.addInteractor(outflowIntr);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(fngIntrWhole);
         intHelper.selectBlocks();

         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 4");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }
         //////////////////////////////////////

         {
            WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
            ppblocks.process(2);
         }

         unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
         int ghostLayer = 3;
         unsigned long long numberOfNodesPerBlock = (unsigned long long)(blockNx[0])* (unsigned long long)(blockNx[1])* (unsigned long long)(blockNx[2]);
         unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
         unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blockNx[0]+ghostLayer) * (blockNx[1]+ghostLayer) * (blockNx[2]+ghostLayer);
         double needMemAll = double(numberOfNodesPerBlockWithGhostLayer*(27*sizeof(double)+sizeof(int)+sizeof(float)*4));
         double needMem = needMemAll/double(comm->getNumberOfProcesses());

         if (myid==0)
         {
            UBLOG(logINFO, "Number of blocks = "<<numberOfBlocks);
            UBLOG(logINFO, "Number of nodes  = "<<numberOfNodes);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for (int level = minInitLevel; level<=maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO, "Number of blocks for level "<<level<<" = "<<nobl);
               UBLOG(logINFO, "Number of nodes for level "<<level<<" = "<<nobl*numberOfNodesPerBlock);
            }
            UBLOG(logINFO, "Necessary memory  = "<<needMemAll<<" bytes");
            UBLOG(logINFO, "Necessary memory per process = "<<needMem<<" bytes");
            UBLOG(logINFO, "Available memory per process = "<<availMem<<" bytes");
         }


         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (myid==0) UBLOG(logINFO, "SetKernelBlockVisitor:end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 5");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         if (refineLevel>0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         if (myid==0) UBLOG(logINFO, "SetUndefinedNodesBlockVisitor:end");

         //BC
         intHelper.setBC();
         if (myid==0) UBLOG(logINFO, "intHelper.setBC():end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 6");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor1;
         initVisitor1.setVx1(fct);
         grid->accept(initVisitor1);

         ////set connectors
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         //PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         //grid->accept(pqPartVisitor);

         //bcVisitor should be accept after initialization!!!!
         grid->accept(bcVisitor);
         if (myid==0) UBLOG(logINFO, "grid->accept(bcVisitor):end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 7");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         //Post process
         {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm);
            ppgeo.process(0);
         }

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 8");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         ////sponge layer
         ////////////////////////////////////////////////////////////////////////////

         GbCuboid3DPtr spongeLayerX1max(new GbCuboid3D(g_maxX1-8.0*blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1max.get(), pathOut+"/geo/spongeLayerX1max", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX1max;
         slVisitorX1max.setBoundingBox(spongeLayerX1max);
         SPtr<LBMKernel> spKernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
         dynamicPointerCast<CompressibleCumulantLBMKernel>(spKernel)->setRelaxationParameter(CompressibleCumulantLBMKernel::NORMAL);
         spKernel->setBCProcessor(bcProc);
         slVisitorX1max.setKernel(spKernel);
         grid->accept(slVisitorX1max);

         //GbCuboid3DPtr spongeLayerX1min(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1+75, g_maxX2+blockLength, g_maxX3+blockLength));
         //if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1min.get(), pathOut+"/geo/spongeLayerX1min", WbWriterVtkXmlASCII::getInstance());
         //SpongeLayerBlockVisitor slVisitorX1min(spongeLayerX1min);
         //grid->accept(slVisitorX1min);

         //GbCuboid3DPtr spongeLayerX3min(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+75));
         //if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX3min.get(), pathOut+"/geo/spongeLayerX3min", WbWriterVtkXmlASCII::getInstance());
         //SpongeLayerBlockVisitor slVisitorX3min(spongeLayerX3min);
         //grid->accept(slVisitorX3min);

         //GbCuboid3DPtr spongeLayerX3max(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3-75, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         //if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX3max.get(), pathOut+"/geo/spongeLayerX3max", WbWriterVtkXmlASCII::getInstance());
         //SpongeLayerBlockVisitor slVisitorX3max(spongeLayerX3max);
         //grid->accept(slVisitorX3max);

         /////////////////////////////////////////////////////////////////////////////
         if (myid==0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         restartCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);
         ////////////////////////////////////////////////////////////////////////////
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         grid->accept(bcVisitor);

         ////sponge layer
         ////////////////////////////////////////////////////////////////////////////

         GbCuboid3DPtr spongeLayerX1max(new GbCuboid3D(g_maxX1-8.0*blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1max.get(), pathOut+"/geo/spongeLayerX1max", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX1max;
         slVisitorX1max.setBoundingBox(spongeLayerX1max);
         SPtr<LBMKernel> spKernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
         dynamicPointerCast<CompressibleCumulantLBMKernel>(spKernel)->setRelaxationParameter(CompressibleCumulantLBMKernel::NORMAL);
         spKernel->setBCProcessor(bcProc);
         slVisitorX1max.setKernel(spKernel);
         grid->accept(slVisitorX1max);
      }

      SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      std::shared_ptr<NUPSCounterCoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> stepSch(new UbScheduler(outTimeStep, outTimeStart));

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 9");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }

      SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      SPtr<UbScheduler> tavSch(new UbScheduler(1, timeAvStart, timeAvStop));
      SPtr<TimeAveragedValuesCoProcessor> tav(new TimeAveragedValuesCoProcessor(grid, pathOut, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,
         TimeAveragedValuesCoProcessor::Density | TimeAveragedValuesCoProcessor::Velocity | TimeAveragedValuesCoProcessor::Fluctuations));
      tav->setWithGhostLayer(true);

      SPtr<IntegrateValuesHelper> mic1(new IntegrateValuesHelper(grid, comm, 300-deltaXcoarse, 35, -600-deltaXcoarse, 300, 65, -600));
      if (myid==0) GbSystem3D::writeGeoObject(mic1->getBoundingBox().get(), pathOut+"/geo/mic1", WbWriterVtkXmlBinary::getInstance());
      SPtr<UbScheduler> stepMV(new UbScheduler(1));
      //TimeseriesCoProcessor tsp1(grid, stepMV, mic1, pathOut+"/mic/mic1", comm);

      omp_set_num_threads(numOfThreads);
      SPtr<Calculator> calculator(new BasicCalculator(grid, tavSch, endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(restartCoProcessor);
      calculator->addCoProcessor(writeMQCoProcessor);
      calculator->addCoProcessor(tav);
      

      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      if (myid==0) UBLOG(logINFO, "Simulation-end");

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 10");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }
   }
   catch (std::exception& e)
   {
      cerr<<e.what()<<endl<<flush;
   }
   catch (std::string& s)
   {
      cerr<<s<<endl;
   }
   catch (...)
   {
      cerr<<"unknown exception"<<endl;
   }

}

void test_run()
{
   try
   {

      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 1");
         UBLOG(logINFO, "PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }

      double g_minX1 = 0;
      double g_minX2 = 0;
      double g_minX3 = 0;

      double g_maxX1 = 5;
      double g_maxX2 = 5;
      double g_maxX3 = 5;

      int blockNx[3] ={ 5, 5, 5 };

      string pathOut = "d:/temp/DLR-F16-Solid-test";

      double deltaX = 1;
      double rhoLB = 0.0;
      double uLB = 0.0866025;
      double nuLB = 0.001; //4.33013e-06;
      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());
      ////////////////////////////////////////////////////////////////////////
      //Grid
      //////////////////////////////////////////////////////////////////////////
      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setDeltaX(deltaX);
      grid->setBlockNX(blockNx[0], blockNx[1], blockNx[2]);
      grid->setPeriodicX1(false);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(false);

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid==0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());
      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
      ppblocks.process(0);

      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      kernel->setNX(std::array<int,3>{{blockNx[0], blockNx[1], blockNx[2]}});
      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new BCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nuLB, 1e9, 12);
      grid->accept(kernelVisitor);

      //initialization of distributions
      InitDistributionsBlockVisitor initVisitor1;
      initVisitor1.setVx1(0.001);
      grid->accept(initVisitor1);

      SPtr<UbScheduler> stepSch(new UbScheduler(1));
      SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      //omp_set_num_threads(numOfThreads);
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepSch, 2));
      calculator->addCoProcessor(writeMQCoProcessor);
   

      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      if (myid==0) UBLOG(logINFO, "Simulation-end");
   }
   catch (std::exception& e)
   {
      cerr<<e.what()<<endl<<flush;
   }
   catch (std::string& s)
   {
      cerr<<s<<endl;
   }
   catch (...)
   {
      cerr<<"unknown exception"<<endl;
   }
}

int main(int argc, char* argv[])
{

   if (argv!=NULL)
   {
      if (argv[1]!=NULL)
      {
         run(string(argv[1]));
      }
      else
      {
         cout<<"Configuration file must be set!: "<<argv[0]<<" <config file>"<<endl<<std::flush;
      }
   }

   //test_run();

   return 0;
}


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
      string          fngFileWhole1 = config.getValue<string>("fngFileWhole1");
      string          fngFileWhole2 = config.getValue<string>("fngFileWhole2");
      //string          tapeFile = config.getValue<string>("tapeFile");
      int             accuracy = config.getValue<int>("accuracy");
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
      //string          pathReInit = config.getValue<string>("pathReInit");
      //int             stepReInit = config.getValue<int>("stepReInit");
      //bool            reinit = config.getValue<bool>("reinit");

      //double          pcpStart = config.getValue<double>("pcpStart");
      //double          pcpStop  = config.getValue<double>("pcpStop");

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

      double g_minX1 = boundingBox[0];//*1000.0;
      double g_minX2 = boundingBox[2];//*1000.0;
      double g_minX3 = boundingBox[4];//*1000.0;
      double g_maxX1 = boundingBox[1];//*1000.0;
      double g_maxX2 = boundingBox[3];//*1000.0;
      double g_maxX3 = boundingBox[5];//*1000.0;
      //deltaXfine *=1000.0;

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
      double uReal = 50; //m/s
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
      double lLB = lReal/deltaXcoarse;
      //double nuLB = (uLB*lLB)/Re; //0.005;
      //double nuLB = 0.005;

      if (myid==0) UBLOG(logINFO, unitConverter.toString());
     

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      //SPtr<GbObject3D> mic6(new GbCuboid3D(0.3, 0.015, -0.46+4.25*deltaXcoarse, 0.3+deltaXcoarse, 0.015+deltaXcoarse, -0.46+5.25*deltaXcoarse));
      //if (myid==0) GbSystem3D::writeGeoObject(mic6.get(), pathOut+"/geo/mic6", WbWriterVtkXmlBinary::getInstance());

      //GbCuboid3DPtr spongeLayerX1max(new GbCuboid3D(g_maxX1-0.35, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      //if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1max.get(), pathOut+"/geo/spongeLayerX1max", WbWriterVtkXmlASCII::getInstance());
      //return;


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

      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new BCProcessor());

      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      //t = 21.8, P = 1.0145 atm, Relative Humidity = 45.8, Second Coefficient of Viscosity = 3120      //Ash, R. L., Zuckerwar, A. J., & Zheng, Z. (1991). Second coefficient of viscosity in air.
      double bulckViscosity = 3120 * nuLB;
      dynamicPointerCast<CompressibleCumulant4thOrderViscosityLBMKernel>(kernel)->setBulkViscosity(bulckViscosity);

      kernel->setBCProcessor(bcProc);

      SPtr<LBMKernel> spKernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
      spKernel->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////
      //restart
      //SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
      //SPtr<MPIIORestartCoProcessor> restartCoProcessor(new MPIIORestartCoProcessor(grid, rSch, pathOut, comm));
      //restartCoProcessor->setLBMKernel(kernel);
      //restartCoProcessor->setBCProcessor(bcProc);

      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIOMigrationCoProcessor> migCoProcessor(new MPIIOMigrationCoProcessor(grid, mSch, pathOut+"/mig", comm));
      migCoProcessor->setLBMKernel(kernel);
      migCoProcessor->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 2");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }

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


         //GbCuboid3DPtr spongeLayerX1max(new GbCuboid3D(g_maxX1-0.4, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         //if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1max.get(), pathOut+"/geo/spongeLayerX1max", WbWriterVtkXmlASCII::getInstance());

         if (myid==0)
         {
            UBLOG(logINFO, "Preprocessing - start");
         }

         {
            WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
            ppblocks.process(0);
         }
         
         //SPtr<GbObject3D> fngMeshWhole(new GbCylinder3D(15.0, 0.0, 0.0, 15.0, 100.0, 0.0, 25.0));
         //GbSystem3D::writeGeoObject(fngMeshWhole.get(), pathOut + "/geo/fngMeshWholeCylinder", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbTriFaceMesh3D> fngMeshWhole1;
         if (myid==0) UBLOG(logINFO, "Read fngFileWhole1:start");
         fngMeshWhole1 = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileWhole1, "fngMeshWhole1", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         if (myid==0) UBLOG(logINFO, "Read fngFileWhole1:end");
         fngMeshWhole1->rotate(0.0, 0.5, 0.0);
         //fngMeshWhole->scale(1e3,1e3,1e3);
         //fngMeshWhole->translate(1.932008e-5-149.867,-0.03-49.95,-0.0172298-1.32814);
         if (myid==0) GbSystem3D::writeGeoObject(fngMeshWhole1.get(), pathOut+"/geo/fngMeshWhole1", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbTriFaceMesh3D> fngMeshWhole2;
         if (myid==0) UBLOG(logINFO, "Read fngFileWhole2:start");
         fngMeshWhole2 = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileWhole2, "fngMeshWhole2", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         if (myid==0) UBLOG(logINFO, "Read fngFileWhole2:end");
         fngMeshWhole2->rotate(0.0, 0.5, 0.0);
         //fngMeshWhole->scale(1e3,1e3,1e3);
         //fngMeshWhole->translate(1.932008e-5-149.867,-0.03-49.95,-0.0172298-1.32814);
         if (myid==0) GbSystem3D::writeGeoObject(fngMeshWhole2.get(), pathOut+"/geo/fngMeshWhole2", WbWriterVtkXmlBinary::getInstance());

         //SPtr<GbTriFaceMesh3D> tapeMesh;
         //if (myid==0) UBLOG(logINFO, "Read fngFileWhole:start");
         //tapeMesh = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+tapeFile, "tapeMesh", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         //if (myid==0) UBLOG(logINFO, "Read fngFileWhole:end");
         //tapeMesh->rotate(0.0, 0.5, 0.0);
         ////fngMeshWhole->scale(1e3,1e3,1e3);
         //tapeMesh->translate(0.0,0.0,-0.001085);
         //if (myid==0) GbSystem3D::writeGeoObject(tapeMesh.get(), pathOut+"/geo/tapeMesh", WbWriterVtkXmlBinary::getInstance());

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 3");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }
         //////////////////////////////////////////////////////////////////////////
         SPtr<Interactor3D> fngIntrWhole1 = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshWhole1, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));//, Interactor3D::POINTS));
         SPtr<Interactor3D> fngIntrWhole2 = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshWhole2, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));

         if (refineLevel>0 && myid==0 && writeBlocks)
         {
            if (myid==0) UBLOG(logINFO, "Refinement - start");
            int rank = grid->getRank();
            grid->setRank(0);


            int level;

            level = 1;
            if (refineLevel - level >= 0)
            {
               dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole1)->refineBlockGridToLevel(level, startDistance, refineDistance);
               //SPtr<GbObject3D> refCylinderL1(new GbCylinder3D(0.3, -0.03, 0.001, 0.3, 0.06, 0.001, 0.1));
               //RefineCrossAndInsideGbObjectBlockVisitor refVisitorCylinderL1(refCylinderL1, level);
               //grid->accept(refVisitorCylinderL1);

               //SPtr<GbObject3D> refBoxL1(new GbCuboid3D(0.15, -0.03, -0.035, 0.45, 0.06, 0.035));
               //if (myid==0) GbSystem3D::writeGeoObject(refBoxL1.get(), pathOut+"/geo/refBoxL1", WbWriterVtkXmlASCII::getInstance());
               //RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL1(refBoxL1, level);
               //grid->accept(refVisitorBoxL1);
            }
            
            level = 2;
            if (refineLevel - level >= 0)
            {
               dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole1)->refineBlockGridToLevel(level, startDistance, refineDistance);
            }

            level = 3;
            if (refineLevel - level >= 0)
            {
               //SPtr<GbObject3D> refCylinderL1(new GbCylinder3D(0.015, -0.03, 0.0, 0.015, 0.06, 0.0, 0.03));
               //GbSystem3D::writeGeoObject(refCylinderL1.get(), pathOut + "/geo/refCylinderL1", WbWriterVtkXmlBinary::getInstance());
               //RefineCrossAndInsideGbObjectBlockVisitor refVisitorCylinderL1(refCylinderL1, level);
               //grid->accept(refVisitorCylinderL1);

               ////SPtr<GbObject3D> refBoxL1(new GbCuboid3D(0.015, -0.03, -0.03, 1.100, 0.06, 0.03));
               ////SPtr<GbObject3D> refBoxL1(new GbCuboid3D(0.12, -0.02625, -0.03, 1.0, 0.06, 0.03));
               //SPtr<GbObject3D> refBoxL1(new GbCuboid3D(0.12, -0.02625, -0.03, 0.34, 0.06, 0.03));
               //if (myid==0) GbSystem3D::writeGeoObject(refBoxL1.get(), pathOut+"/geo/refBoxL1", WbWriterVtkXmlASCII::getInstance());
               //RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL1(refBoxL1, level);
               //grid->accept(refVisitorBoxL1);

               //SPtr<GbObject3D> refCylinderL1(new GbCylinder3D(0.3, -0.03, 0.001, 0.3, 0.06, 0.001, 0.03));
               //GbSystem3D::writeGeoObject(refCylinderL1.get(), pathOut + "/geo/refCylinderL1", WbWriterVtkXmlBinary::getInstance());
               //RefineCrossAndInsideGbObjectBlockVisitor refVisitorCylinderL1(refCylinderL1, level);
               //grid->accept(refVisitorCylinderL1);

               SPtr<GbObject3D> refBoxL2(new GbCuboid3D(0.15, -0.03, -0.015, 0.42, 0.06, 0.015));
               if (myid==0) GbSystem3D::writeGeoObject(refBoxL2.get(), pathOut+"/geo/refBoxL2", WbWriterVtkXmlASCII::getInstance());
               RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL2(refBoxL2, level);
               grid->accept(refVisitorBoxL2);

               //dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole1)->refineBlockGridToLevel(level, startDistance, 24.0*refineDistance);
               dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole1)->refineBlockGridToLevel(level, startDistance, 12.0*refineDistance);
            }

            //level = 4;
            //if (refineLevel - level >= 0)
            //{
            //   //SPtr<GbObject3D> refCylinderL2(new GbCylinder3D(0.015, -0.03, 0.0, 0.015, 0.06, 0.0, 0.03));
            //   //GbSystem3D::writeGeoObject(refCylinderL2.get(), pathOut + "/geo/refCylinderL2", WbWriterVtkXmlBinary::getInstance());
            //   //RefineCrossAndInsideGbObjectBlockVisitor refVisitorCylinderL2(refCylinderL2, level);
            //   //grid->accept(refVisitorCylinderL2);

            //   //SPtr<GbObject3D> refBoxL2(new GbCuboid3D(0.15, -0.03, -0.015, 0.7, 0.06, 0.015));
            //   SPtr<GbObject3D> refBoxL2(new GbCuboid3D(0.15, -0.03, -0.015, 0.42, 0.06, 0.015));
            //   if (myid==0) GbSystem3D::writeGeoObject(refBoxL2.get(), pathOut+"/geo/refBoxL2", WbWriterVtkXmlASCII::getInstance());
            //   RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL2(refBoxL2, level);
            //   grid->accept(refVisitorBoxL2);

            //   //SPtr<GbObject3D> refCylinderL1(new GbCylinder3D(0.3, -0.03, 0.001, 0.3, 0.06, 0.001, 0.03));
            //   //GbSystem3D::writeGeoObject(refCylinderL1.get(), pathOut + "/geo/refCylinderL1", WbWriterVtkXmlBinary::getInstance());
            //   //RefineCrossAndInsideGbObjectBlockVisitor refVisitorCylinderL1(refCylinderL1, level);
            //   //grid->accept(refVisitorCylinderL1);

            //   dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole1)->refineBlockGridToLevel(level, startDistance, 6.0*refineDistance);
            //}

            //level = 5;
            //if (refineLevel - level >= 0)
            //{
            //   //SPtr<GbObject3D> refCylinderL3(new GbCylinder3D(0.015, -0.03, 0.0, 0.015, 0.06, 0.0, 0.025));
            //   //GbSystem3D::writeGeoObject(refCylinderL3.get(), pathOut + "/geo/refCylinderL3", WbWriterVtkXmlBinary::getInstance());
            //   //RefineCrossAndInsideGbObjectBlockVisitor refVisitorCylinderL3(refCylinderL3, level);
            //   //grid->accept(refVisitorCylinderL3);

            //   //SPtr<GbObject3D> refBoxL3(new GbCuboid3D(0.15, -0.03, -0.010, 0.32, 0.06, 0.012));
            //   //if (myid==0) GbSystem3D::writeGeoObject(refBoxL3.get(), pathOut+"/geo/refBoxL3", WbWriterVtkXmlASCII::getInstance());
            //   //RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL3(refBoxL3, level);
            //   //grid->accept(refVisitorBoxL3);

            //   dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole1)->refineBlockGridToLevel(level, startDistance, 6.0*refineDistance);
            //}

            //level = 6;
            //if (refineLevel - level >= 0)
            //{
            //   //SPtr<GbObject3D> refBoxL4(new GbCuboid3D(0.15, -0.03, -0.005, 0.31, 0.06, 0.01));
            //   //if (myid==0) GbSystem3D::writeGeoObject(refBoxL4.get(), pathOut+"/geo/refBoxL4", WbWriterVtkXmlASCII::getInstance());
            //   //RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL4(refBoxL4, level);
            //   //grid->accept(refVisitorBoxL4);

            //   dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole1)->refineBlockGridToLevel(level, startDistance, 3.0*refineDistance);
            //}

            //level = 7;
            //if (refineLevel - level >= 0)
            //{
            //   dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole1)->refineBlockGridToLevel(level, startDistance, refineDistance);
            //}

            //dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrWhole1)->refineBlockGridToLevel(refineLevel, startDistance, refineDistance);
            //

            /////delete solid blocks
            if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");

            SetSolidBlocksBlockVisitor v(fngIntrWhole1);
            grid->accept(v);
            std::vector<SPtr<Block3D>>& sb = fngIntrWhole1->getSolidBlockSet();
            for (SPtr<Block3D> block : sb)
            {
               grid->deleteBlock(block);
            }
            fngIntrWhole1->removeSolidBlocks();
            fngIntrWhole1->removeBcBlocks();

            //SPtr<GbObject3D> delBox(new GbCuboid3D(0.03, -0.03, -0.010, 0.2, 0.06, 0.012));
            //if (myid==0) GbSystem3D::writeGeoObject(delBox.get(), pathOut+"/geo/delBox", WbWriterVtkXmlASCII::getInstance());
            //SPtr<D3Q27Interactor> delBoxInter(new D3Q27Interactor(delBox, grid, noSlipBCAdapter, Interactor3D::SOLID));
            //SetSolidBlockVisitor v(delBoxInter, BlockType::SOLID);
            //grid->accept(v);
            //std::vector<SPtr<Block3D>>& sb = delBoxInter->getSolidBlockSet();
            //for (SPtr<Block3D> block : sb)
            //{
            //   grid->deleteBlock(block);
            //}
            //delBoxInter->removeSolidBlocks();
            //delBoxInter->removeBcBlocks();

            if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - end");
            ////////////////////////////////////////

            {
               WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(1);
            }

            grid->setRank(rank);

            RatioBlockVisitor ratioVisitor(refineLevel);
            CheckRatioBlockVisitor checkRatio(refineLevel);
            int count = 0;

            do {
               if (myid==0) UBLOG(logINFO, "ratioVisitor - start");
               grid->accept(ratioVisitor);
               if (myid==0) UBLOG(logINFO, "ratioVisitor - end");
               if (myid==0) UBLOG(logINFO, "checkRatio - start");
               checkRatio.resetState();
               grid->accept(checkRatio);
               if (myid==0) UBLOG(logINFO, "checkRatio - end");
               if (myid==0) UBLOG(logINFO, "count = "<<count++<<" state = "<<checkRatio.getState());
            } while (!checkRatio.getState());


            {
               WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(2);
            }

            OverlapBlockVisitor overlapVisitor(refineLevel, false);
            grid->accept(overlapVisitor);

            if (myid==0) UBLOG(logINFO, "Refinement - end");
         }
         //else
         //{
         //   migCoProcessor->readBlocks(0);
         //}
         grid->updateDistributedBlocks(comm);

         //if (writeBlocks)
         //{
         //   migCoProcessor->writeBlocks(0);
         //}

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
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, MetisPartitioner::RECURSIVE));
         //std::dynamic_pointer_cast<MetisPartitioningGridVisitor>(metisVisitor)->setNumberOfProcesses(4000);
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(inflowIntr);
         intHelper.addInteractor(outflowIntr);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(fngIntrWhole2);
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
            ppblocks.process(3);
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
         InitDistributionsBlockVisitor initVisitor;
         initVisitor.setVx1(fct);
         grid->accept(initVisitor);

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 7");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         //Post process
         {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
            ppgeo.process(0);
         }

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 8");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         if (myid==0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         //restartCoProcessor->restart((int)restartStep);
         migCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);
         ////////////////////////////////////////////////////////////////////////////
         WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
         ppblocks.process(0);
         ////////////////////////////////////////////////////////////////////////////
      }
      
      ////set connectors
      SPtr<InterpolationProcessor> iProcessor(new CompressibleOffsetMomentsInterpolationProcessor());
      dynamicPointerCast<CompressibleOffsetMomentsInterpolationProcessor>(iProcessor)->setBulkViscosity(nuLB, bulckViscosity);
      SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
      grid->accept(setConnsVisitor);

      //bcVisitor should be accept after initialization!!!!
      grid->accept(bcVisitor);
      if (myid == 0) UBLOG(logINFO, "grid->accept(bcVisitor):end");

      ////sponge layer
      GbCuboid3DPtr spongeLayerX1max(new GbCuboid3D(g_maxX1-0.35, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1max.get(), pathOut+"/geo/spongeLayerX1max", WbWriterVtkXmlASCII::getInstance());
      SpongeLayerBlockVisitor slVisitorX1max(spongeLayerX1max, spKernel, nuLB, D3Q27System::E);
      grid->accept(slVisitorX1max);

      SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      std::shared_ptr<CoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> stepSch(new UbScheduler(outTimeStep, outTimeStart));

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 9");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }

      SPtr<CoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      SPtr<GbObject3D> bbBox(new GbCuboid3D(g_minX1-blockLength, (g_maxX2-g_minX2)/2.0, g_minX3-blockLength, g_maxX1+blockLength, (g_maxX2-g_minX2)/2.0+deltaXcoarse, g_maxX3+blockLength));
      if (myid==0) GbSystem3D::writeGeoObject(bbBox.get(), pathOut+"/geo/bbBox", WbWriterVtkXmlASCII::getInstance());
      SPtr<WriteMQFromSelectionCoProcessor> writeMQSelectCoProcessor(new WriteMQFromSelectionCoProcessor(grid, stepSch, bbBox, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      SPtr<UbScheduler> tavSch(new UbScheduler(1, timeAvStart, timeAvStop));
      SPtr<TimeAveragedValuesCoProcessor> tav(new TimeAveragedValuesCoProcessor(grid, pathOut, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,
         TimeAveragedValuesCoProcessor::Density | TimeAveragedValuesCoProcessor::Velocity | TimeAveragedValuesCoProcessor::Fluctuations));
      tav->setWithGhostLayer(true);

      SPtr<UbScheduler> stepMV(new UbScheduler(1, 0, 1000000));
      SPtr<MicrophoneArrayCoProcessor> micCoProcessor(new MicrophoneArrayCoProcessor(grid, stepSch, pathOut, comm) );
      double offsetX1 = 0.017;
      double offsetZ1 = 0.11375;
      std::vector<UbTupleFloat3> nodes;
      //for (int i = 0; i <= 10; i++)
      for (int i = 0; i <= 6; i++)
      {
         micCoProcessor->addMicrophone(Vector3D(0.3+deltaXcoarse+offsetX1*double(i), 0.015, 0.0-offsetZ1*double(i)));
         nodes.push_back(UbTupleFloat3(float(0.3+deltaXcoarse+offsetX1*float(i)), float(0.015), float(0.0-offsetZ1*float(i))));
      }
      //double offsetX2 = 0.1;
      //for (int i = 0; i <= 6; i++)
      //{
      //   micCoProcessor->addMicrophone(Vector3D(0.17+offsetX2*double(i), 0.015, -1.1375));
      //   nodes.push_back(UbTupleFloat3(float(0.17+offsetX2*float(i)), float(0.015), float(-1.1375)));
      //}
      
      if (myid==0) WbWriterVtkXmlBinary::getInstance()->writeNodes(pathOut+"/geo/mic", nodes);

      //omp_set_num_threads(numOfThreads);
      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(micCoProcessor);
      //calculator->addCoProcessor(restartCoProcessor);
      calculator->addCoProcessor(migCoProcessor);
      //calculator->addCoProcessor(writeMQSelectCoProcessor);
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


int main(int argc, char* argv[])
{
   //Sleep(30000);

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
   
   //SuperMUC
   //MPI_Finalize();

   return 0;
}


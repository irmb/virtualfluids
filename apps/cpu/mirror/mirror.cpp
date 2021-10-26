#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;

void run(string configname)
{
   try
   {
      vf::basics::ConfigurationFile   config;
      config.load(configname);

      string          pathOut = config.getValue<string>("pathOut");
      string          pathGeo = config.getValue<string>("pathGeo");
      string          pathMesh = config.getValue<string>("pathMesh");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blockNx = config.getVector<int>("blockNx");
      double          restartStep = config.getValue<double>("restartStep");
      double          restartStepStart = config.getValue<double>("restartStepStart");
      double          endTime = config.getValue<double>("endTime");
      double          outTime = config.getValue<double>("outTime");
      double          availMem = config.getValue<double>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      bool            logToFile = config.getValue<bool>("logToFile");
      double          deltaXcoarse = config.getValue<double>("deltaXcoarse");
      double          deltaXfine = config.getValue<double>("deltaXfine");
      double          refineDistance = config.getValue<double>("refineDistance");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");

      vector<double>  WTUNNEL1 = config.getVector<double>("WTUNNEL1");
      vector<double>  VRES0100 = config.getVector<double>("VRES0100");
      vector<double>  VRES0200 = config.getVector<double>("VRES0200");
      vector<double>  VRES0300 = config.getVector<double>("VRES0300");
      vector<double>  VRES0400 = config.getVector<double>("VRES0400");
      vector<double>  VRES0500 = config.getVector<double>("VRES0500");
      vector<double>  VRES0700 = config.getVector<double>("VRES0700");
      vector<double>  VRES0900 = config.getVector<double>("VRES0900");

      string          SAE = config.getValue<string>("SAE");
      string          VRES0600_chopped = config.getValue<string>("VRES0600_chopped");
      string          VRES0700_chopped = config.getValue<string>("VRES0700_chopped");
      string          VRES0800_Fahrzeug = config.getValue<string>("VRES0800_Fahrzeug");
      //string          VRES0900 = config.getValue<string>("VRES0900");
      string          VRES1000_ASaeule = config.getValue<string>("VRES1000_ASaeule");
      string          VRES1000_Scheibe = config.getValue<string>("VRES1000_Scheibe");
      string          VRES1000_Spiegel = config.getValue<string>("VRES1000_Spiegel");
      string          VRES1100_Spiegel_fein = config.getValue<string>("VRES1100_Spiegel_fein");


      SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
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


      double g_minX1 = WTUNNEL1[0];
      double g_minX2 = WTUNNEL1[2];
      double g_minX3 = WTUNNEL1[4];

      double g_maxX1 = WTUNNEL1[1];
      double g_maxX2 = WTUNNEL1[3];
      double g_maxX3 = WTUNNEL1[5];

      double blockLength = (double)blockNx[0]*deltaXcoarse;

      //##########################################################################
      //## physical parameters
      //##########################################################################


      double rhoLB = 0.0;
      double rhoReal = 1.2041; //(kg/m3)
      double nueReal = 153.5e-7; //m^2/s

      double lReal = 2.048;//m
      double uReal = 140.0/3.6;

      double Re = uReal*lReal/nueReal;

      //##Machzahl:
      //#Ma     = uReal/csReal
      double Ma = 140.0/1236.0;//Ma-Real!

      double uLB = Ma*sqrt(1.0/3.0);
      double nuLB = (uLB*1.0)/Re;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      ////////////////////////////////////////////////////////////////////////
      //Grid
      //////////////////////////////////////////////////////////////////////////
      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setDeltaX(deltaXcoarse);
      grid->setBlockNX(blockNx[0], blockNx[1], blockNx[2]);

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid==0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());
      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      grid->setPeriodicX1(false);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(false);

      //BC adapters
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      SPtr<BCAdapter> slipBCAdapter(new SlipBCAdapter());
      slipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SlipBCAlgorithm()));

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);
      SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapter(new DensityBCAdapter(rhoLB));
      denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(slipBCAdapter);
      bcVisitor.addBC(velBCAdapter);
      bcVisitor.addBC(denBCAdapter);

      //////////////////////////////////////////////////////////////////////////
      //restart
    
      //////////////////////////////////////////////////////////////////////////


      if (grid->getTimeStep()==0)
      {
         if (myid==0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "* Re                  = "<<Re);
            UBLOG(logINFO, "* Ma                  = "<<Ma);
            UBLOG(logINFO, "* velocity (uReal)    = "<<uReal<<" m/s");
            UBLOG(logINFO, "* viscosity (nuReal)  = "<<nueReal<<" m^2/s");
            UBLOG(logINFO, "* velocity LB (uLB)   = "<<uLB);
            UBLOG(logINFO, "* viscosity LB (nuLB) = "<<nuLB);
            UBLOG(logINFO, "* dx_base             = "<<deltaXcoarse<<" m");
            UBLOG(logINFO, "* dx_refine           = "<<deltaXfine<<" m");
            UBLOG(logINFO, "* number of levels    = "<<refineLevel+1);
            UBLOG(logINFO, "* number of threads   = "<<numOfThreads);
            UBLOG(logINFO, "* number of processes = "<<comm->getNumberOfProcesses());
            UBLOG(logINFO, "Preprozess - start");
         }

         GbCuboid3DPtr geoVRES0100(new GbCuboid3D(VRES0100[0], VRES0100[2], VRES0100[4], VRES0100[1], VRES0100[3], VRES0100[5]));
         if (myid==0) GbSystem3D::writeGeoObject(geoVRES0100.get(), pathOut+"/geo/geoVRES0100", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr geoVRES0200(new GbCuboid3D(VRES0200[0], VRES0200[2], VRES0200[4], VRES0200[1], VRES0200[3], VRES0200[5]));
         if (myid==0) GbSystem3D::writeGeoObject(geoVRES0200.get(), pathOut+"/geo/geoVRES0200", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr geoVRES0300(new GbCuboid3D(VRES0300[0], VRES0300[2], VRES0300[4], VRES0300[1], VRES0300[3], VRES0300[5]));
         if (myid==0) GbSystem3D::writeGeoObject(geoVRES0300.get(), pathOut+"/geo/geoVRES0300", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr geoVRES0400(new GbCuboid3D(VRES0400[0], VRES0400[2], VRES0400[4], VRES0400[1], VRES0400[3], VRES0400[5]));
         if (myid==0) GbSystem3D::writeGeoObject(geoVRES0400.get(), pathOut+"/geo/geoVRES0400", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr geoVRES0500(new GbCuboid3D(VRES0500[0], VRES0500[2], VRES0500[4], VRES0500[1], VRES0500[3], VRES0500[5]));
         if (myid==0) GbSystem3D::writeGeoObject(geoVRES0500.get(), pathOut+"/geo/geoVRES0500", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr geoVRES0700(new GbCuboid3D(VRES0700[0], VRES0700[2], VRES0700[4], VRES0700[1], VRES0700[3], VRES0700[5]));
         if (myid==0) GbSystem3D::writeGeoObject(geoVRES0700.get(), pathOut+"/geo/geoVRES0700", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr geoVRES0900(new GbCuboid3D(VRES0900[0], VRES0900[2], VRES0900[4], VRES0900[1], VRES0900[3], VRES0900[5]));
         if (myid==0) GbSystem3D::writeGeoObject(geoVRES0900.get(), pathOut+"/geo/geoVRES0900", WbWriterVtkXmlASCII::getInstance());

         SPtr<D3Q27Interactor> geoVRES0700Int(new D3Q27Interactor(geoVRES0700, grid, noSlipBCAdapter, Interactor3D::SOLID));

         //GEO
         if (myid==0) UBLOG(logINFO, "Read geoSAE:start");
         SPtr<GbTriFaceMesh3D> geoSAE = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+SAE, "meshSAE", GbTriFaceMesh3D::KDTREE_SAHPLIT, true));
         if (myid==0) UBLOG(logINFO, "Read meshSAE:end");
         if (myid==0) GbSystem3D::writeGeoObject(geoSAE.get(), pathOut+"/geo/meshSAE", WbWriterVtkXmlBinary::getInstance());

         SPtr<D3Q27TriFaceMeshInteractor> geoSAEInteractor(new D3Q27TriFaceMeshInteractor(geoSAE, grid, noSlipBCAdapter, Interactor3D::SOLID));



         if (myid==0)
         {
            //////////////////////////////////////////
            //meshes
            if (myid==0) UBLOG(logINFO, "Read meshVRES0600:start");
            SPtr<GbTriFaceMesh3D> meshVRES0600 = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathMesh+"/"+VRES0600_chopped, "meshVRES0600", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read meshVRES0600:end");
            if (myid==0) GbSystem3D::writeGeoObject(meshVRES0600.get(), pathOut+"/geo/meshVRES0600", WbWriterVtkXmlBinary::getInstance());
            SPtr<D3Q27TriFaceMeshInteractor> meshVRES0600Interactor(new D3Q27TriFaceMeshInteractor(meshVRES0600, grid, noSlipBCAdapter, Interactor3D::SOLID));

            if (myid==0) UBLOG(logINFO, "Read meshVRES0700:start");
            SPtr<GbTriFaceMesh3D> meshVRES0700 = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathMesh+"/"+VRES0700_chopped, "meshVRES0700", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read meshVRES0700:end");
            if (myid==0) GbSystem3D::writeGeoObject(meshVRES0700.get(), pathOut+"/geo/meshVRES0700", WbWriterVtkXmlBinary::getInstance());
            SPtr<D3Q27TriFaceMeshInteractor> meshVRES0700Interactor(new D3Q27TriFaceMeshInteractor(meshVRES0700, grid, noSlipBCAdapter, Interactor3D::SOLID));

            if (myid==0) UBLOG(logINFO, "Read meshVRES0800:start");
            SPtr<GbTriFaceMesh3D> meshVRES0800 = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathMesh+"/"+VRES0800_Fahrzeug, "meshVRES0800", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read meshVRES0800:end");
            if (myid==0) GbSystem3D::writeGeoObject(meshVRES0800.get(), pathOut+"/geo/meshVRES0800", WbWriterVtkXmlBinary::getInstance());
            SPtr<D3Q27TriFaceMeshInteractor> meshVRES0800Interactor(new D3Q27TriFaceMeshInteractor(meshVRES0800, grid, noSlipBCAdapter, Interactor3D::SOLID));

            //if (myid==0) UBLOG(logINFO, "Read meshVRES0900:start");
            //SPtr<GbTriFaceMesh3D> meshVRES0900 = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathMesh+"/"+VRES0900, "meshVRES0900", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            //if (myid==0) UBLOG(logINFO, "Read meshVRES0900:end");
            //if (myid==0) GbSystem3D::writeGeoObject(meshVRES0900.get(), pathOut+"/geo/meshVRES0900", WbWriterVtkXmlBinary::getInstance());
            //SPtr<D3Q27TriFaceMeshInteractor> meshVRES0900Interactor(new D3Q27TriFaceMeshInteractor(meshVRES0900, grid, noSlipBCAdapter, Interactor3D::SOLID));

            if (myid==0) UBLOG(logINFO, "Read meshVRES1000ASaeule:start");
            SPtr<GbTriFaceMesh3D> meshVRES1000ASaeule = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathMesh+"/"+VRES1000_ASaeule, "meshVRES1000ASaeule", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read meshVRES1000ASaeule:end");
            if (myid==0) GbSystem3D::writeGeoObject(meshVRES1000ASaeule.get(), pathOut+"/geo/meshVRES1000ASaeule", WbWriterVtkXmlBinary::getInstance());
            SPtr<D3Q27TriFaceMeshInteractor> meshVRES1000ASaeuleInteractor(new D3Q27TriFaceMeshInteractor(meshVRES1000ASaeule, grid, noSlipBCAdapter, Interactor3D::SOLID));

            if (myid==0) UBLOG(logINFO, "Read meshVRES1000Scheibe:start");
            SPtr<GbTriFaceMesh3D> meshVRES1000Scheibe = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathMesh+"/"+VRES1000_Scheibe, "meshVRES1000Scheibe", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read meshVRES1000Scheibe:end");
            if (myid==0) GbSystem3D::writeGeoObject(meshVRES1000Scheibe.get(), pathOut+"/geo/meshVRES1000Scheibe", WbWriterVtkXmlBinary::getInstance());
            SPtr<D3Q27TriFaceMeshInteractor> meshVRES1000ScheibeInteractor(new D3Q27TriFaceMeshInteractor(meshVRES1000Scheibe, grid, noSlipBCAdapter, Interactor3D::SOLID));

            if (myid==0) UBLOG(logINFO, "Read meshVRES1000Spiegel:start");
            SPtr<GbTriFaceMesh3D> meshVRES1000Spiegel = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathMesh+"/"+VRES1000_Spiegel, "meshSpiegel", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read meshVRES1000Spiegel:end");
            if (myid==0) GbSystem3D::writeGeoObject(meshVRES1000Spiegel.get(), pathOut+"/geo/meshVRES1000Spiegel", WbWriterVtkXmlBinary::getInstance());
            SPtr<D3Q27TriFaceMeshInteractor> meshVRES1000SpiegelInteractor(new D3Q27TriFaceMeshInteractor(meshVRES1000Spiegel, grid, noSlipBCAdapter, Interactor3D::SOLID));

            if (myid==0) UBLOG(logINFO, "Read meshVRES1100SpiegelFine:start");
            SPtr<GbTriFaceMesh3D> meshVRES1100SpiegelFine = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathMesh+"/"+VRES1100_Spiegel_fein, "meshSpiegelFine", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read meshVRES1100SpiegelFine:end");
            if (myid==0) GbSystem3D::writeGeoObject(meshVRES1100SpiegelFine.get(), pathOut+"/geo/meshVRES1100SpiegelFine", WbWriterVtkXmlBinary::getInstance());
            SPtr<D3Q27TriFaceMeshInteractor> meshVRES1100SpiegelFineInteractor(new D3Q27TriFaceMeshInteractor(meshVRES1100SpiegelFine, grid, noSlipBCAdapter, Interactor3D::SOLID));

            UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
            ////refineHelper.addGbObject(geoVRES0100, refineLevel-4);
            //refineHelper.addGbObject(geoVRES0200, 1);
            //refineHelper.addGbObject(geoVRES0300, 2);
            //refineHelper.addGbObject(geoVRES0400, 3);
            //refineHelper.addGbObject(geoVRES0500, 4);
            //refineHelper.addGbObject(geoVRES0700, 7);
            //refineHelper.addGbObject(geoVRES0900, 9);
            //refineHelper.refine();

            RefineCrossAndInsideGbObjectBlockVisitor geoVRES0200RefVisitor(geoVRES0200, 1);
            grid->accept(geoVRES0200RefVisitor);
            RefineCrossAndInsideGbObjectBlockVisitor geoVRES0300RefVisitor(geoVRES0300, 2);
            grid->accept(geoVRES0300RefVisitor);
            RefineCrossAndInsideGbObjectBlockVisitor geoVRES0400RefVisitor(geoVRES0400, 3);
            grid->accept(geoVRES0400RefVisitor);
            RefineCrossAndInsideGbObjectBlockVisitor geoVRES0500RefVisitor(geoVRES0500, 4);
            grid->accept(geoVRES0500RefVisitor);


            int rank = grid->getRank();
            grid->setRank(0);
            meshVRES0600Interactor->refineBlockGridToLevel(5, 0.0, 0.0);
            meshVRES0700Interactor->refineBlockGridToLevel(6, -0.6, 0.0);

            UBLOG(logINFO, "Refinement - geoVRES0700");
            RefineCrossAndInsideGbObjectBlockVisitor geoVRES0700RefVisitor(geoVRES0700, 7);
            grid->accept(geoVRES0700RefVisitor);

            UBLOG(logINFO, "Refinement - geoSAEInteractor");
            meshVRES0800Interactor->refineBlockGridToLevel(8, -0.5, 0.0);
            //geoSAEInteractor->refineBlockGridToLevel(8, 0.0, 0.1);

            //SetSolidOrTransBlockVisitor v(geoSAEInteractor, SetSolidOrTransBlockVisitor::SOLID);
            //grid->accept(v);
            //std::vector<SPtr<Block3D>>& sb = geoSAEInteractor->getSolidBlockSet();
            //BOOST_FOREACH(SPtr<Block3D> block, sb)
            //{
            //   grid->deleteBlock(block);
            //}
            //geoSAEInteractor->removeSolidBlocks();
            //geoSAEInteractor->removeTransBlocks();

            UBLOG(logINFO, "Refinement - geoVRES0900RefVisitor");
            //meshVRES0900Interactor->refineBlockGridToLevel(9, 0.0, 0.0);
            RefineCrossAndInsideGbObjectBlockVisitor geoVRES0900RefVisitor(geoVRES0900, 9);
            grid->accept(geoVRES0900RefVisitor);

            UBLOG(logINFO, "Refinement - meshVRES1000ASaeuleInteractor");
            meshVRES1000ASaeuleInteractor->refineBlockGridToLevel(10, -0.1, 0.0);

            UBLOG(logINFO, "Refinement - meshVRES1000ScheibeInteractor");
            meshVRES1000ScheibeInteractor->refineBlockGridToLevel(10, -0.1, 0.0);

            UBLOG(logINFO, "Refinement - meshVRES1000SpiegelInteractor");
            meshVRES1000SpiegelInteractor->refineBlockGridToLevel(10, -0.12, 0.0);

            UBLOG(logINFO, "Refinement - meshVRES1100SpiegelFineInteractor");
            meshVRES1100SpiegelFineInteractor->refineBlockGridToLevel(11, -0.12, 0.0);
            grid->setRank(rank);

            ///////////////////////////////////////////////////////////
            ///BOX
            //GbCuboid3DPtr geoBox1(new GbCuboid3D(-0.495, -0.8, 0.545, -0.045, -0.7, 0.795));
            //if (myid==0) GbSystem3D::writeGeoObject(geoBox1.get(), pathOut+"/geo/geoBox1", WbWriterVtkXmlASCII::getInstance());
            //CoarsenCrossAndInsideGbObjectBlockVisitor geoBox1Visitor(geoBox1, 11, 11);
            //grid->accept(geoBox1Visitor);
            //////////////////////////////////////////////////////////////////////////


            if (myid==0)
            {
               WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(0);
            }

            RatioBlockVisitor ratioVisitor(refineLevel);
            CheckRatioBlockVisitor checkRatio(refineLevel);
            int count = 0;

            do {
               UBLOG(logINFO, "Refinement - RatioBlockVisitor");
               grid->accept(ratioVisitor);
               checkRatio.resetState();
               UBLOG(logINFO, "Refinement - CheckRatioBlockVisitor");
               grid->accept(checkRatio);
               if (myid==0) UBLOG(logINFO, "count ="<<count++<<" state="<<checkRatio.getState());
            } while (!checkRatio.getState());

            UBLOG(logINFO, "Refinement - OverlapBlockVisitor");
            OverlapBlockVisitor overlapVisitor(refineLevel, false);
            grid->accept(overlapVisitor);

            if (myid==0) UBLOG(logINFO, "Refinement - end");

            if (myid==0)
            {
               WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(1);
            }
         }

         grid->updateDistributedBlocks(comm);

         if (myid == 0) UBLOG(logINFO, "SetInterpolationDirsBlockVisitor");
         std::vector<int> dirs;
         for (int i = D3Q27System::E; i<=D3Q27System::TS; i++)
         {
            dirs.push_back(i);
         }
         SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
         grid->accept(interDirsVisitor);

         //////////////////////////////////////////////////////////////////////////


         //walls
         GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_minX2, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathOut+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1-blockLength, g_maxX2, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathOut+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         //wall interactors
         SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, slipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, slipBCAdapter, Interactor3D::SOLID));

         //walls
         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
         if (myid==0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathOut+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathOut+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //wall interactors
         SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, slipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, slipBCAdapter, Interactor3D::SOLID));

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathOut+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         //inflow
         SPtr<D3Q27Interactor> inflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //outflow
         SPtr<D3Q27Interactor> outflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(inflowIntr);
         intHelper.addInteractor(outflowIntr);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(geoVRES0700Int);
         intHelper.addInteractor(geoSAEInteractor);
         //////////////////////////////////////////////////////////////////////////
         intHelper.selectBlocks();

         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         if (myid==0)
         {
            SPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm));
            ppblocks->process(2);
            ppblocks.reset();
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

         SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel(blockNx[0], blockNx[1], blockNx[2], CompressibleCumulantLBMKernel::NORMAL));

         SPtr<BCProcessor> bcProc;

         bcProc = SPtr<BCProcessor>(new BCProcessor());

         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel>0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //BC
         intHelper.setBC();

         grid->accept(bcVisitor);

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         initVisitor.setVx1(fct);
         initVisitor.setNu(nuLB);
         grid->accept(initVisitor);

         ////set connectors
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //Postrozess
         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         SPtr<CoProcessor> ppgeo(new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         ppgeo->process(0);
         ppgeo.reset();

         if (myid==0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         grid->accept(bcVisitor);
      }

      SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      SPtr<UbScheduler> stepSch(new UbScheduler(outTime));

      WriteMacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm);

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe());
      }

      const SPtr<ConcreteCalculatorFactory> calculatorFactory = std::make_shared<ConcreteCalculatorFactory>(stepSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, calculatorFactory, CalculatorType::HYBRID));
      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
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

   return 0;
}


#include <iostream>
#include <string>

#include <boost/pointer_cast.hpp>

#include "VirtualFluids.h"

using namespace std;

double rangeRandom1()
{
   return (2.0*rand())/RAND_MAX-1.0;
}

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
      string          zigZagTape = config.getString("zigZagTape");
      int             numOfThreads = config.getInt("numOfThreads");
      vector<int>     blockNx = config.getVector<int>("blockNx");
      vector<double>  boundingBox = config.getVector<double>("boundingBox");
      double          uLB = config.getDouble("uLB");
      double          restartStep = config.getDouble("restartStep");
      double          restartStepStart = config.getDouble("restartStepStart");
      double          endTime = config.getDouble("endTime");
      double          outTime = config.getDouble("outTime");
      double          availMem = config.getDouble("availMem");
      int             refineLevel = config.getInt("refineLevel");
      bool            logToFile = config.getBool("logToFile");
      bool            porousTralingEdge = config.getBool("porousTralingEdge");
      double          deltaXfine = config.getDouble("deltaXfine")*1000.0;
      bool            thinWall = config.getBool("thinWall");
      double          refineDistance = config.getDouble("refineDistance");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");

      SPtr<Communicator> comm = MPICommunicator::getInstance();
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

      
      double g_minX1 = boundingBox[0]*1000.0;
      double g_minX2 = boundingBox[2]*1000.0;
      double g_minX3 = boundingBox[4]*1000.0;

      double g_maxX1 = boundingBox[1]*1000.0;
      double g_maxX2 = boundingBox[3]*1000.0;
      double g_maxX3 = boundingBox[5]*1000.0;
       
      //////////////////////////////////////////////////////////////////////////
      double deltaXcoarse = deltaXfine*(double)(1 << refineLevel);
      //double nx2_temp = floor((g_maxX2 - g_minX2) / (deltaXcoarse*(double)blockNx[0]));

      //deltaXcoarse = (g_maxX2 - g_minX2) / (nx2_temp*(double)blockNx[0]);
      //UBLOG(logINFO, "nx2_temp:"<<nx2_temp);
      //g_maxX2 -= 0.5* deltaXcoarse;
      //////////////////////////////////////////////////////////////////////////
      double blockLength = (double)blockNx[0] * deltaXcoarse;

      //##########################################################################
      //## physical parameters
      //##########################################################################
      double Re = 1;//e6;

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
      double l_LB = 300 / deltaXcoarse;
      double nuLB = (uLB*l_LB) / Re; //0.005;
      //double nuLB = 0.005;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      ////////////////////////////////////////////////////////////////////////
      //Grid
      //////////////////////////////////////////////////////////////////////////
      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setDeltaX(deltaXcoarse);
      grid->setBlockNX(blockNx[0], blockNx[1], blockNx[2]);

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      //gridCube->setCenterCoordinates(geo->getX1Centroid(), geo->getX2Centroid(), geo->getX3Centroid());
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut + "/geo/gridCube", WbWriterVtkXmlASCII::getInstance());
      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      grid->setPeriodicX1(false);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(false);

      //BC adapters
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      if (thinWall)
      {
         noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new ThinWallNoSlipBCAlgorithm()));
      }
      else
      {
         noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
      }

      SPtr<BCAdapter> slipBCAdapter(new SlipBCAdapter());
      slipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SlipBCAlgorithm()));

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);
      SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapter(new DensityBCAdapter(rhoLB));
      denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(slipBCAdapter);
      bcVisitor.addBC(velBCAdapter);
      bcVisitor.addBC(denBCAdapter);

      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> rSch(new UbScheduler(restartStep, restartStep));
      RestartCoProcessor rp(grid, rSch, comm, pathOut, RestartCoProcessor::TXT);
      //////////////////////////////////////////////////////////////////////////


      if (grid->getTimeStep() == 0)
      {
         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "* Re                  = "<<Re);
            UBLOG(logINFO, "* Ma                  = "<<Ma);
            UBLOG(logINFO, "* velocity (uReal)    = "<<uReal<<" m/s");
            UBLOG(logINFO, "* viscosity (nuReal)  = "<<nueReal<<" m^2/s");
            UBLOG(logINFO, "* velocity LB (uLB)   = "<<uLB);
            UBLOG(logINFO, "* viscosity LB (nuLB) = "<<nuLB);
            UBLOG(logINFO, "* dx_base             = "<<deltaXcoarse/1000.0<<" m");
            UBLOG(logINFO, "* dx_refine           = "<<deltaXfine/1000.0<<" m");
            UBLOG(logINFO, "* number of levels    = " << refineLevel + 1);
            UBLOG(logINFO, "* number of threads   = " << numOfThreads);
            UBLOG(logINFO, "* number of processes = " << comm->getNumberOfProcesses());
            UBLOG(logINFO, "Preprozess - start");
         }

         SPtr<GbTriFaceMesh3D> fngMeshWhole;
         SPtr<GbTriFaceMesh3D> fngMeshBodyPart;
         SPtr<GbTriFaceMesh3D> fngMeshTrailingEdge;
         if (porousTralingEdge)
         {
            if (myid==0) UBLOG(logINFO, "Read fngFileBodyPart:start");
            fngMeshBodyPart = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileBodyPart, "fngMeshBody", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read fngFileBodyPart:end");
            fngMeshBodyPart->rotate(0.0, 0.5, 0.0);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshBodyPart.get(), pathOut+"/geo/fngMeshBody", WbWriterVtkXmlBinary::getInstance());

            if (myid==0) UBLOG(logINFO, "Read fngFileTrailingEdge:start");
            fngMeshTrailingEdge = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileTrailingEdge, "fngMeshTrailingEdge", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read fngFileTrailingEdge:end");
            fngMeshTrailingEdge->rotate(0.0, 0.5, 0.0);
            fngMeshTrailingEdge->translate(0,0,1.3);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshTrailingEdge.get(), pathOut+"/geo/fngMeshTrailingEdge", WbWriterVtkXmlBinary::getInstance());
         }
         else
         {
            if (myid==0) UBLOG(logINFO, "Read fngFileWhole:start");
            fngMeshWhole = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileWhole, "fngMeshWhole", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read fngFileWhole:end");
            fngMeshWhole->rotate(0.0, 0.5, 0.0);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshWhole.get(), pathOut+"/geo/fngMeshWhole", WbWriterVtkXmlBinary::getInstance());
         }

         //////////////////////////////////////////////////////////////////////////
         // Zackenband
         //////////////////////////////////////////////////////////////////////////
         //top
         //////////////////////////////////////////////////////////////////////////
         //if (myid==0) UBLOG(logINFO, "Read zigZagTape:start");
         //string ZckbndFilename = pathGeo+"/"+zigZagTape;
         //SPtr<GbTriFaceMesh3D> meshBand1(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape1"));
         //meshBand1->rotate(0.0, 5, 0.0);
         //meshBand1->translate(15, 0, -12.850);
         //if (myid==0) GbSystem3D::writeGeoObject(meshBand1.get(), pathOut+"/geo/zigZagTape1", WbWriterVtkXmlASCII::getInstance());
         //// Zackenband2
         //SPtr<GbTriFaceMesh3D> meshBand2(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape2"));
         //meshBand2->rotate(0.0, 5, 0.0);
         //meshBand2->translate(15, 5, -12.850);
         //if (myid==0) GbSystem3D::writeGeoObject(meshBand2.get(), pathOut+"/geo/zigZagTape2", WbWriterVtkXmlASCII::getInstance());
         ////// Zackenband3
         ////SPtr<GbTriFaceMesh3D> meshBand3(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape3"));
         ////meshBand3->rotate(0.0, 5, 0.0);
         ////meshBand3->translate(15, 0, -12.35);
         ////if (myid==0) GbSystem3D::writeGeoObject(meshBand3.get(), pathOut+"/geo/zigZagTape3", WbWriterVtkXmlASCII::getInstance());
         ////// Zackenband4
         ////SPtr<GbTriFaceMesh3D> meshBand4(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape4"));
         ////meshBand4->rotate(0.0, 5, 0.0);
         ////meshBand4->translate(15, 5, -12.35);
         ////if (myid==0) GbSystem3D::writeGeoObject(meshBand4.get(), pathOut+"/geo/zigZagTape4", WbWriterVtkXmlASCII::getInstance());
        
         ////bottom
         //SPtr<GbTriFaceMesh3D> meshBand5(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape5"));
         //meshBand5->rotate(0.0, -1, 0.0);
         //meshBand5->rotate(0.0, 0.0,180.0);
         //meshBand5->translate(30, 0, -37.3);
         //if (myid==0) GbSystem3D::writeGeoObject(meshBand5.get(), pathOut+"/geo/zigZagTape5", WbWriterVtkXmlASCII::getInstance());
         //// Zackenband6
         //SPtr<GbTriFaceMesh3D> meshBand6(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape6"));
         //meshBand6->rotate(0.0, -1, 0.0);
         //meshBand6->rotate(0.0, 0.0, 180.0);
         //meshBand6->translate(30, 5, -37.3);
         //if (myid==0) GbSystem3D::writeGeoObject(meshBand6.get(), pathOut+"/geo/zigZagTape6", WbWriterVtkXmlASCII::getInstance());
         ////// Zackenband7
         ////SPtr<GbTriFaceMesh3D> meshBand7(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape7"));
         ////meshBand7->rotate(0.0, 5, 0.0);
         ////meshBand7->translate(15, 0, -12.35);
         ////if (myid==0) GbSystem3D::writeGeoObject(meshBand7.get(), pathOut+"/geo/zigZagTape7", WbWriterVtkXmlASCII::getInstance());
         ////// Zackenband8
         ////SPtr<GbTriFaceMesh3D> meshBan8(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape8"));
         ////meshBan8->rotate(0.0, 5, 0.0);
         ////meshBan8->translate(15, 5, -12.35);
         ////if (myid==0) GbSystem3D::writeGeoObject(meshBan8.get(), pathOut+"/geo/zigZagTape8", WbWriterVtkXmlASCII::getInstance());
         //if (myid==0) UBLOG(logINFO, "Read zigZagTape:end");

         //////////////////////////////////////////////////////////////////////////

         SPtr<Interactor3D> fngIntrWhole;
         SPtr<Interactor3D> fngIntrBodyPart;
         SPtr<Interactor3D> fngIntrTrailingEdge;
         if (porousTralingEdge)
         {
            fngIntrBodyPart = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshBodyPart, grid, noSlipBCAdapter, Interactor3D::SOLID));
            fngIntrTrailingEdge = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshTrailingEdge, grid, noSlipBCAdapter, Interactor3D::SOLID));
         }
         else
         {
            fngIntrWhole = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshWhole, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));
         }

         //SPtr<D3Q27TriFaceMeshInteractor> triBand1Interactor(new D3Q27TriFaceMeshInteractor(meshBand1, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));
         //SPtr<D3Q27TriFaceMeshInteractor> triBand2Interactor(new D3Q27TriFaceMeshInteractor(meshBand2, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));
         //SPtr<D3Q27TriFaceMeshInteractor> triBand3Interactor(new D3Q27TriFaceMeshInteractor(meshBand3, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));
         //SPtr<D3Q27TriFaceMeshInteractor> triBand4Interactor(new D3Q27TriFaceMeshInteractor(meshBand4, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));



         if (refineLevel > 0 && myid == 0)
         {
            if (myid == 0) UBLOG(logINFO, "Refinement - start");
            //RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            //refineHelper.addGbObject(geo, refineLevel);
            //refineHelper.refine();
            
            //RefineAroundGbObjectHelper refineHelper1(grid, refineLevel-1, boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(geoIntr1), 0.0, 10.0, comm);
            //refineHelper1.refine();
            //RefineAroundGbObjectHelper refineHelper2(grid, refineLevel, boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(geoIntr2), -1.0, 5.0, comm);
            //refineHelper2.refine();
            

            int rank = grid->getRank();
            grid->setRank(0);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand1Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand2Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand3Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand4Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            grid->setRank(rank);

            if (porousTralingEdge)
            {
               int rank = grid->getRank();
               grid->setRank(0);
               boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrBodyPart)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
               grid->setRank(rank);
            }
            else
            {
               int rank = grid->getRank();
               grid->setRank(0);
               boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrWhole)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
               grid->setRank(rank);
            }



            ////////////////////////////////////////////
            //METIS
            //SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
            ////////////////////////////////////////////
            /////delete solid blocks
            if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
            //InteractorsHelper intHelper(grid, metisVisitor);
            //if (porousTralingEdge)
            //{
            //   intHelper.addInteractor(fngIntrBodyPart);
            //}
            //else
            //{
            //   intHelper.addInteractor(fngIntrWhole);
            //}
            //////////////////////////////////////////////////////////////////////////
            
            //intHelper.selectBlocks();
            
            if (porousTralingEdge)
            {
               SetSolidBlockVisitor v(fngIntrBodyPart, BlockType::SOLID);
               grid->accept(v);
               std::vector<SPtr<Block3D>>& sb = fngIntrBodyPart->getSolidBlockSet();
               for(SPtr<Block3D> block : sb)
               {
                  grid->deleteBlock(block);
               }
               fngIntrBodyPart->removeSolidBlocks();
               fngIntrBodyPart->removeBcBlocks();
            }
            else
            {
               SetSolidBlockVisitor v(fngIntrWhole, BlockType::SOLID);
               grid->accept(v);
               std::vector<SPtr<Block3D>>& sb = fngIntrWhole->getSolidBlockSet();
               for(SPtr<Block3D> block : sb)
               {
                  grid->deleteBlock(block);
               }
               fngIntrWhole->removeSolidBlocks();
               fngIntrWhole->removeBcBlocks();
            }

            if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
            //////////////////////////////////////

            if (porousTralingEdge)
            {
               grid->setRank(0);
               boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrTrailingEdge)->refineBlockGridToLevel(refineLevel, -2.0, refineDistance);
               grid->setRank(rank);

               //SPtr<GbObject3D> trailingEdgeCube(new GbCuboid3D(fngMeshTrailingEdge->getX1Minimum()-blockLength, fngMeshTrailingEdge->getX2Minimum(), fngMeshTrailingEdge->getX3Minimum()-blockLength/2.0,
               //   fngMeshTrailingEdge->getX1Maximum()+blockLength, fngMeshTrailingEdge->getX2Maximum(), fngMeshTrailingEdge->getX3Maximum()+blockLength/2.0));
               //if (myid == 0) GbSystem3D::writeGeoObject(trailingEdgeCube.get(), pathOut + "/geo/trailingEdgeCube", WbWriterVtkXmlASCII::getInstance());

               //RefineCrossAndInsideGbObjectBlockVisitor refVisitor(trailingEdgeCube, refineLevel);
               //grid->accept(refVisitor);
            }

            RatioBlockVisitor ratioVisitor(refineLevel);
            CheckRatioBlockVisitor checkRatio(refineLevel);
            int count = 0;
            
            do {
               grid->accept(ratioVisitor);
               checkRatio.resetState();
               grid->accept(checkRatio);
               if (myid == 0) UBLOG(logINFO, "count ="<<count++<<" state="<<checkRatio.getState());
            } while (!checkRatio.getState());

            //RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
            //grid->accept(ratioSmoothVisitor);

            {
               WriteBlocksSPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm));
               ppblocks->process(0);
               ppblocks.reset();
            }

            OverlapBlockVisitor overlapVisitor(refineLevel, false);
            grid->accept(overlapVisitor);

            //std::vector<int> dirs;
            //for (int i = D3Q27System::E; i <= D3Q27System::TS; i++)
            //{
            //   dirs.push_back(i);
            //}
            //SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
            //grid->accept(interDirsVisitor);

            if (myid == 0) UBLOG(logINFO, "Refinement - end");
         }

         grid->updateDistributedBlocks(comm);


         //return;

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
         SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, slipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, slipBCAdapter, Interactor3D::SOLID));
         //SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         //SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));

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
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(inflowIntr);
         intHelper.addInteractor(outflowIntr);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         //intHelper.addInteractor(triBand1Interactor);
         //intHelper.addInteractor(triBand2Interactor);
         //intHelper.addInteractor(triBand3Interactor);
         //intHelper.addInteractor(triBand4Interactor);
         
         if (porousTralingEdge)
         {
            intHelper.addInteractor(fngIntrBodyPart);
            //intHelper.addInteractor(fngIntrTrailingEdge);
         } 
         else
         {
            intHelper.addInteractor(fngIntrWhole);
         }
         
         //////////////////////////////////////////////////////////////////////////
         intHelper.selectBlocks();

         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         WriteBlocksSPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm));
         ppblocks->process(1);
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

         SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel(blockNx[0], blockNx[1], blockNx[2], CompressibleCumulantLBMKernel::NORMAL));
         //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel(blockNx[0], blockNx[1], blockNx[2], IncompressibleCumulantLBMKernel::NORMAL));

         SPtr<BCProcessor> bcProc;

         if (thinWall)
         {
            bcProc = SPtr<BCProcessor>(new ThinWallBCProcessor());
         }
         else
         {
            bcProc = SPtr<BCProcessor>(new BCProcessor());
         }

         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //BC
         intHelper.setBC();
         
         grid->accept(bcVisitor);

         //initialization of distributions
         mu::Parser inflowProfileVx1, inflowProfileVx2, inflowProfileVx3;
         inflowProfileVx1.SetExpr("U*rangeRandom1()");
         inflowProfileVx1.DefineConst("U", uLB);
         inflowProfileVx1.DefineFun("rangeRandom1", rangeRandom1);
         inflowProfileVx2.SetExpr("0.1*U*rangeRandom1()");
         inflowProfileVx2.DefineConst("U", uLB);
         inflowProfileVx2.DefineFun("rangeRandom1", rangeRandom1);
         inflowProfileVx3.SetExpr("0.1*U*rangeRandom1()");
         inflowProfileVx3.DefineConst("U", uLB);
         inflowProfileVx3.DefineFun("rangeRandom1", rangeRandom1);
         
         InitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         //initVisitor.setVx1(fct);
         //initVisitor.setVx1(inflowProfileVx1);
         //initVisitor.setVx2(inflowProfileVx2);
         //initVisitor.setVx3(inflowProfileVx3);
         //initVisitor.setNu(nuLB);
         grid->accept(initVisitor);

         ////set connectors
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         //InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //Postrozess
         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         WriteBoundaryConditionsSPtr<CoProcessor> ppgeo(
            new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         ppgeo->process(0);
         ppgeo.reset();

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
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

      WriteMacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv,comm);

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      const SPtr<ConcreteCalculatorFactory> calculatorFactory = std::make_shared<ConcreteCalculatorFactory>(stepSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, calculatorFactory, CalculatorType::HYBRID));
      //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch, CalculationManager::PrePostBc));
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


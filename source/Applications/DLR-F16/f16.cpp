#include <iostream>
#include <string>

#include <boost/pointer_cast.hpp>

#include "VirtualFluids.h"
#include <omp.h>
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
      double          startDistance = config.getDouble("startDistance");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");
      bool            newStart = config.getBool("newStart");

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

      //the geometry is in mm

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
      //double csReal = uReal / Ma;
      //double hLB = lReal / deltaXcoarse;

      //LBMUnitConverter unitConverter(lReal, csReal, rhoReal, hLB);

      //double u_LB = uReal   * unitConverter.getFactorVelocityWToLb();
      //double nu_LB = nueReal * unitConverter.getFactorViscosityWToLb();
      double lLB = lReal*1000.0 / deltaXcoarse;
      double nuLB = (uLB*lLB)/Re; //0.005;
      //double nuLB = 0.005;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;

      ////////////////////////////////////////////////////////////////////////
      //Grid
      //////////////////////////////////////////////////////////////////////////
      Grid3DPtr grid(new Grid3D(comm));

      //BC adapters
      BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
      if (thinWall)
      {
         noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new ThinWallNoSlipBCAlgorithm()));
      }
      else
      {
         noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithm()));
      }

      BCAdapterPtr slipBCAdapter(new SlipBCAdapter());
      slipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new SlipBCAlgorithm()));

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);
      BCAdapterPtr velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingVelocityBCAlgorithm()));

      BCAdapterPtr denBCAdapter(new DensityBCAdapter(rhoLB));
      denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingDensityBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(slipBCAdapter);
      bcVisitor.addBC(velBCAdapter);
      bcVisitor.addBC(denBCAdapter);

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(restartStep, restartStepStart));
      //RestartCoProcessor rp(grid, rSch, comm, pathOut, RestartCoProcessor::TXT);
      MPIIORestartCoProcessor rcp(grid, rSch, pathOut, comm);
      //////////////////////////////////////////////////////////////////////////


      //if (grid->getTimeStep() == 0)
      if(newStart)
      {
         ////////////////////////////////////////////////////////////////////////
         //define grid
         //////////////////////////////////////////////////////////////////////////
         grid->setDeltaX(deltaXcoarse);
         grid->setBlockNX(blockNx[0], blockNx[1], blockNx[2]);

         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid==0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         grid->setPeriodicX1(false);
         grid->setPeriodicX2(true);
         grid->setPeriodicX3(false);

         if (myid == 0)
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
            UBLOG(logINFO, "* dx_base             = "<<deltaXcoarse/1000<<" m");
            UBLOG(logINFO, "* dx_refine           = "<<deltaXfine/1000<<" m");
            UBLOG(logINFO, "* blocknx             = "<<blockNx[0]<<"x"<<blockNx[1]<<"x"<<blockNx[2] );
            UBLOG(logINFO, "* refineDistance      = "<<refineDistance);
            UBLOG(logINFO, "* number of levels    = "<<refineLevel + 1);
            UBLOG(logINFO, "* number of threads   = "<<numOfThreads);
            UBLOG(logINFO, "* number of processes = "<<comm->getNumberOfProcesses());
            UBLOG(logINFO, "Preprozess - start");
         }

         GbTriFaceMesh3DPtr fngMeshWhole;
         GbTriFaceMesh3DPtr fngMeshBodyPart;
         GbTriFaceMesh3DPtr fngMeshTrailingEdge;
         if (porousTralingEdge)
         {
            if (myid==0) UBLOG(logINFO, "Read fngFileBodyPart:start");
            fngMeshBodyPart = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileBodyPart, "fngMeshBody", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read fngFileBodyPart:end");
            fngMeshBodyPart->rotate(0.0, 0.5, 0.0);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshBodyPart.get(), pathOut+"/geo/fngMeshBody", WbWriterVtkXmlBinary::getInstance());

            if (myid==0) UBLOG(logINFO, "Read fngFileTrailingEdge:start");
            fngMeshTrailingEdge = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileTrailingEdge, "fngMeshTrailingEdge", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read fngFileTrailingEdge:end");
            fngMeshTrailingEdge->rotate(0.0, 0.5, 0.0);
            fngMeshTrailingEdge->translate(0,0,1.3);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshTrailingEdge.get(), pathOut+"/geo/fngMeshTrailingEdge", WbWriterVtkXmlBinary::getInstance());
         }
         else
         {
            if (myid==0) UBLOG(logINFO, "Read fngFileWhole:start");
            fngMeshWhole = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileWhole, "fngMeshWhole", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read fngFileWhole:end");
            fngMeshWhole->rotate(0.0, 0.5, 0.0);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshWhole.get(), pathOut+"/geo/fngMeshWhole", WbWriterVtkXmlBinary::getInstance());
         }

         //////////////////////////////////////////////////////////////////////////
         // Zackenband
         //////////////////////////////////////////////////////////////////////////
         //top
         //////////////////////////////////////////////////////////////////////////
         if (myid==0) UBLOG(logINFO, "Read zigZagTape:start");
         string ZckbndFilename = pathGeo+"/"+zigZagTape;
         GbTriFaceMesh3DPtr meshBand1(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape1"));
         meshBand1->rotate(0.0, 5, 0.0);
         meshBand1->translate(15, 0, -12.850);
         if (myid==0) GbSystem3D::writeGeoObject(meshBand1.get(), pathOut+"/geo/zigZagTape1", WbWriterVtkXmlASCII::getInstance());
         // Zackenband2
         GbTriFaceMesh3DPtr meshBand2(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape2"));
         meshBand2->rotate(0.0, 5, 0.0);
         meshBand2->translate(15, 5, -12.850);
         if (myid==0) GbSystem3D::writeGeoObject(meshBand2.get(), pathOut+"/geo/zigZagTape2", WbWriterVtkXmlASCII::getInstance());
         //// Zackenband3
         //GbTriFaceMesh3DPtr meshBand3(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape3"));
         //meshBand3->rotate(0.0, 5, 0.0);
         //meshBand3->translate(15, 0, -12.35);
         //if (myid==0) GbSystem3D::writeGeoObject(meshBand3.get(), pathOut+"/geo/zigZagTape3", WbWriterVtkXmlASCII::getInstance());
         //// Zackenband4
         //GbTriFaceMesh3DPtr meshBand4(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape4"));
         //meshBand4->rotate(0.0, 5, 0.0);
         //meshBand4->translate(15, 5, -12.35);
         //if (myid==0) GbSystem3D::writeGeoObject(meshBand4.get(), pathOut+"/geo/zigZagTape4", WbWriterVtkXmlASCII::getInstance());
        
         //bottom
         GbTriFaceMesh3DPtr meshBand5(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape5"));
         meshBand5->rotate(0.0, -1, 0.0);
         meshBand5->rotate(0.0, 0.0,180.0);
         //meshBand5->translate(30, 0, -37.3);
         meshBand5->translate(30, 0, -37.2);
         if (myid==0) GbSystem3D::writeGeoObject(meshBand5.get(), pathOut+"/geo/zigZagTape5", WbWriterVtkXmlASCII::getInstance());
         // Zackenband6
         GbTriFaceMesh3DPtr meshBand6(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape6"));
         meshBand6->rotate(0.0, -1, 0.0);
         meshBand6->rotate(0.0, 0.0, 180.0);
         //meshBand6->translate(30, 5, -37.3);
         meshBand6->translate(30, 5, -37.2);
         if (myid==0) GbSystem3D::writeGeoObject(meshBand6.get(), pathOut+"/geo/zigZagTape6", WbWriterVtkXmlASCII::getInstance());
         //// Zackenband7
         //GbTriFaceMesh3DPtr meshBand7(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape7"));
         //meshBand7->rotate(0.0, 5, 0.0);
         //meshBand7->translate(15, 0, -12.35);
         //if (myid==0) GbSystem3D::writeGeoObject(meshBand7.get(), pathOut+"/geo/zigZagTape7", WbWriterVtkXmlASCII::getInstance());
         //// Zackenband8
         //GbTriFaceMesh3DPtr meshBan8(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape8"));
         //meshBan8->rotate(0.0, 5, 0.0);
         //meshBan8->translate(15, 5, -12.35);
         //if (myid==0) GbSystem3D::writeGeoObject(meshBan8.get(), pathOut+"/geo/zigZagTape8", WbWriterVtkXmlASCII::getInstance());
         if (myid==0) UBLOG(logINFO, "Read zigZagTape:end");

         

         //////////////////////////////////////////////////////////////////////////

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
            fngIntrWhole = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshWhole, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));
         }

         D3Q27TriFaceMeshInteractorPtr triBand1Interactor(new D3Q27TriFaceMeshInteractor(meshBand1, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));
         D3Q27TriFaceMeshInteractorPtr triBand2Interactor(new D3Q27TriFaceMeshInteractor(meshBand2, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));
         D3Q27TriFaceMeshInteractorPtr triBand3Interactor(new D3Q27TriFaceMeshInteractor(meshBand5, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));
         D3Q27TriFaceMeshInteractorPtr triBand4Interactor(new D3Q27TriFaceMeshInteractor(meshBand6, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::EDGES));

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

            if (porousTralingEdge)
            {
               boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrBodyPart)->refineBlockGridToLevel(refineLevel, startDistance, refineDistance);
            }
            //else
            //{
            //   boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrWhole)->refineBlockGridToLevel(refineLevel, startDistance, refineDistance);
            //}

            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand1Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand2Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand3Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand4Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);


            GbObject3DPtr fngBox(new GbCuboid3D(fngMeshWhole->getX1Minimum(), fngMeshWhole->getX2Minimum(), fngMeshWhole->getX3Minimum(),
                                                fngMeshWhole->getX1Maximum(), fngMeshWhole->getX2Maximum(), fngMeshWhole->getX3Maximum()));
            if (myid==0) GbSystem3D::writeGeoObject(fngBox.get(), pathOut+"/geo/fngBox", WbWriterVtkXmlASCII::getInstance());

            RefineCrossAndInsideGbObjectBlockVisitor refVisitor0(fngBox, refineLevel);
            grid->accept(refVisitor0);

            
            GbObject3DPtr bandTopBox(new GbCuboid3D(meshBand1->getX1Minimum(), meshBand1->getX2Minimum(), meshBand1->getX3Minimum(), 
                                                 meshBand1->getX1Maximum(), meshBand1->getX2Maximum(), meshBand1->getX3Maximum()));
            if (myid==0) GbSystem3D::writeGeoObject(bandTopBox.get(), pathOut+"/geo/bandTopBox", WbWriterVtkXmlASCII::getInstance());

            RefineCrossAndInsideGbObjectBlockVisitor refVisitor1(bandTopBox, refineLevel);
            grid->accept(refVisitor1);

            GbObject3DPtr bandBottomBox(new GbCuboid3D(meshBand5->getX1Minimum(), meshBand5->getX2Minimum(), meshBand5->getX3Minimum(), 
                                                    meshBand5->getX1Maximum(), meshBand5->getX2Maximum(), meshBand5->getX3Maximum()));
            if (myid==0) GbSystem3D::writeGeoObject(bandBottomBox.get(), pathOut+"/geo/bandBottomBox", WbWriterVtkXmlASCII::getInstance());

            RefineCrossAndInsideGbObjectBlockVisitor refVisitor2(bandBottomBox, refineLevel);
            grid->accept(refVisitor2);

            grid->setRank(rank);

            {
               WriteBlocksCoProcessor ppblocks(grid, UbSchedulerPtr(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(0);
            }

            ////////////////////////////////////////////
            //METIS
            //Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
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
               SetSolidOrTransBlockVisitor v(fngIntrBodyPart, SetSolidOrTransBlockVisitor::SOLID);
               grid->accept(v);
               std::vector<Block3DPtr>& sb = fngIntrBodyPart->getSolidBlockSet();
               BOOST_FOREACH(Block3DPtr block, sb)
               {
                  grid->deleteBlock(block);
               }
               fngIntrBodyPart->removeSolidBlocks();
               fngIntrBodyPart->removeTransBlocks();
            }
            else
            {
               SetSolidOrTransBlockVisitor v(fngIntrWhole, SetSolidOrTransBlockVisitor::SOLID);
               grid->accept(v);
               std::vector<Block3DPtr>& sb = fngIntrWhole->getSolidBlockSet();
               BOOST_FOREACH(Block3DPtr block, sb)
               {
                  grid->deleteBlock(block);
               }
               fngIntrWhole->removeSolidBlocks();
               fngIntrWhole->removeTransBlocks();
            }

            if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
            //////////////////////////////////////

            if (porousTralingEdge)
            {
               grid->setRank(0);
               boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrTrailingEdge)->refineBlockGridToLevel(refineLevel, -2.0, refineDistance);
               grid->setRank(rank);

               //GbObject3DPtr trailingEdgeCube(new GbCuboid3D(fngMeshTrailingEdge->getX1Minimum()-blockLength, fngMeshTrailingEdge->getX2Minimum(), fngMeshTrailingEdge->getX3Minimum()-blockLength/2.0,
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
               WriteBlocksCoProcessor ppblocks(grid, UbSchedulerPtr(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(1);
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
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, slipBCAdapter, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, slipBCAdapter, Interactor3D::SOLID));

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathOut+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         //inflow
         D3Q27InteractorPtr inflowIntr = D3Q27InteractorPtr(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //outflow
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
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         //intHelper.addInteractor(triBand1Interactor);
         //intHelper.addInteractor(triBand2Interactor);
         //intHelper.addInteractor(triBand3Interactor);
         //intHelper.addInteractor(triBand4Interactor);
         //
         //if (porousTralingEdge)
         //{
         //   intHelper.addInteractor(fngIntrBodyPart);
         //   //intHelper.addInteractor(fngIntrTrailingEdge);
         //} 
         //else
         //{
         //   intHelper.addInteractor(fngIntrWhole);
         //}
         
         //////////////////////////////////////////////////////////////////////////
         intHelper.selectBlocks();

         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         {
            WriteBlocksCoProcessor ppblocks(grid, UbSchedulerPtr(new UbScheduler(1)), pathOut, WbWriterVtkXmlASCII::getInstance(), comm);
            ppblocks.process(2);
         }

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

         LBMKernelPtr kernel = LBMKernelPtr(new CompressibleCumulantLBMKernel(blockNx[0], blockNx[1], blockNx[2], CompressibleCumulantLBMKernel::NORMAL));
         //LBMKernelPtr kernel = LBMKernelPtr(new IncompressibleCumulantLBMKernel(blockNx[0], blockNx[1], blockNx[2], IncompressibleCumulantLBMKernel::NORMAL));

         BCProcessorPtr bcProc;

         if (thinWall)
         {
            bcProc = BCProcessorPtr(new ThinWallBCProcessor());
         }
         else
         {
            bcProc = BCProcessorPtr(new BCProcessor());
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
         initVisitor.setVx1(fct);
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
         UbSchedulerPtr geoSch(new UbScheduler(1));
         WriteBoundaryConditionsCoProcessorPtr ppgeo(
            new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         ppgeo->process(0);
         ppgeo.reset();

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         rcp.restart();
         grid->setTimeStep(restartStepStart);

         {
            WriteBlocksCoProcessor ppblocks(grid, UbSchedulerPtr(new UbScheduler(1)), pathOut, WbWriterVtkXmlASCII::getInstance(), comm);
            ppblocks.process(3);
         }

         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         grid->accept(bcVisitor);
      }

      UbSchedulerPtr nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      UbSchedulerPtr stepSch(new UbScheduler(outTime));

      WriteMacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv,comm);

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
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


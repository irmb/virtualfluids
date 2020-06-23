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

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      string machine = string(cstr1);

      if(machine == "my") 
      {
         pathname = "d:/temp/bone";
         pathGeo = "d:/Data/Bone/BigBone";
         pathLog = "d:/temp/bone";
         numOfThreads = 4;
         logfile = false;
         availMem = 15.0e9;
      }
      else if(machine == "Ludwig")      
      {
         pathname = "/work/koskuche/SFB880/plate2Con";
         pathGeo = "/home/koskuche/data/plate";
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

      if(myid==0) UBLOG(logINFO,"Testcase bone");

      string boneFilename = pathGeo + "/bone.raw";
      
      int pmNX1=1800;  //abmessung einzelbild in x-richtung
      int pmNX2=972; //abmessung einzelbild in y richtung
      int pmNX3=1164; //anzahl der bilder
      float lthreshold = 27756.0;
      float uthreshold = 65535.0;

      GbVoxelMatrix3DPtr pmMesh(new GbVoxelMatrix3D(pmNX1,pmNX2,pmNX3,0,lthreshold,uthreshold));
      pmMesh->readMatrixFromRawFile<unsigned short>(boneFilename, GbVoxelMatrix3D::BigEndian);

      double scaleFactor = 0.001;
      double delta = 11.0*scaleFactor;
      pmMesh->setVoxelMatrixDelta(delta, delta, delta);

      pmMesh->setVoxelMatrixMininum(0.0, 0.0, 0.0);
      if(myid == 0) pmMesh->writeToLegacyVTKBinary(pathname+"/geo/bone");

      ///////////////////////////////////////////////////////
      return;


      /////////////////Knotenabmessungen:
      //int nx[3], blocknx[3];
      //nx[0]      = 90;//240;//120;//60;//86;//43;//65;//50;  //länge
      //nx[1]      = 2;//2;//6;///1;//5;// //breite
      //nx[2]      = 30;//64;//32;//18;//5;//15;//15; //höhe gebiet
      //blocknx[0] = 16;//10;//6;
      //blocknx[1] = 16;//10;//6;
      //blocknx[2] = 16;//10;//6;

      //int baseLevel   = 0;
      //int refineLevel = 4;

      //double H = 600.0; // Kanalhöhe [mm]
      //double cdx = H/(double)(nx[2]*blocknx[2]);
      //double fdx = cdx/double(1<<refineLevel);

      ////double h = 200.0; // gewünschte Plattenhöhe in Gitterpunkten
      ////double fdx = plate->getLengthX3()/h;
      ////double cdx = fdx*double(1<<refineLevel);

      //LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter());

      ////////////////////////////////////////////////////////////////////////////
      ////physik
      ////////////////////////////////////////////////////////////////////////////
      //double Re            = 1133333.3333333335; 
      //double rhoLB         = 0.0;
      //double uLB           = 0.1; 
      //double lReal         = 1000; //Plattenlänge in mm
      //double nuLB          = (uLB*(lReal/cdx))/Re;

      //int sizeSP=4;
      //mu::Parser spongeLayer;
      //spongeLayer.SetExpr("x1>=(sizeX-sizeSP)/dx ? (sizeX-(x1+1))/sizeSP/2.0 + 0.5 : 1.0");
      //spongeLayer.DefineConst("sizeX", nx[0]*blocknx[0]);
      //spongeLayer.DefineConst("sizeSP", sizeSP*blocknx[0]);

      //Grid3DPtr grid(new Grid3D(comm));

      ////////////////////////////////////////////////////////////////////////////
      ////restart
      //UbSchedulerPtr rSch(new UbScheduler(1000,1000,10000000));
      //RestartPostprocessor rp(grid, rSch, comm, pathname, RestartPostprocessor::BINARY);
      ////////////////////////////////////////////////////////////////////////////

      //if (grid->getTimeStep() == 0)
      //{

      //   if(myid==0) UBLOG(logINFO,"Neustart..");

      //   //////////////////////////////////////////////////////////////////////////
      //   //Platte
      //   GbTriFaceMesh3DPtr plate (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(PlatteFilename,"Netz"));
      //   plate->rotate(90.0,0.0,0.0);  //TriFacMesh-KO-System anders als LB-KO-System
      //   if(myid == 0) GbSystem3D::writeGeoObject( plate.get(), pathname+"/geo/platte", WbWriterVtkXmlBinary::getInstance() );
      //   //////////////////////////////////////////////////////////////////////////
      //   // Zackenband
      //   //////////////////////////////////////////////////////////////////////////
      //   GbTriFaceMesh3DPtr meshBand1 (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "NetzBand"));
      //   meshBand1->translate(-495, -700, -19.94);
      //   if(myid == 0) GbSystem3D::writeGeoObject( meshBand1.get(), pathname+"/geo/Band1", WbWriterVtkXmlASCII::getInstance() );
      //   // Zackenband2
      //   GbTriFaceMesh3DPtr meshBand2(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "NetzBand2"));
      //   meshBand2->translate(-495, -705, -19.94); 
      //   if(myid == 0) GbSystem3D::writeGeoObject( meshBand2.get(), pathname+"/geo/Band2", WbWriterVtkXmlASCII::getInstance() );
      //   // Zackenband3
      //   GbTriFaceMesh3DPtr meshBand3(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "NetzBand3"));
      //   meshBand3->translate(-495, -700, -19.64); 
      //   if(myid == 0) GbSystem3D::writeGeoObject( meshBand3.get(), pathname+"/geo/Band3", WbWriterVtkXmlASCII::getInstance() );
      //   // Zackenband4
      //   GbTriFaceMesh3DPtr meshBand4(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "NetzBand4"));
      //   meshBand4->translate(-495, -705, -19.64); 
      //   if(myid == 0) GbSystem3D::writeGeoObject( meshBand4.get(), pathname+"/geo/Band4", WbWriterVtkXmlASCII::getInstance() );
      //   //////////////////////////////////////////////////////////////////////////

      //   double blockLengthx1 = blocknx[0]*cdx; //geowerte
      //   double blockLengthx2 = blockLengthx1;
      //   double blockLengthx3 = blockLengthx1;

      //   double geoLength[]   = {  nx[0]*blockLengthx1, nx[1]*blockLengthx2, nx[2]*blockLengthx3}; 

      //   double originX1 = plate->getX1Minimum()-plate->getLengthX1()/4.0;
      //   double originX2 = plate->getX2Minimum();
      //   double originX3 = plate->getX3Minimum()-299.5;


      //   bool periodicx1 = false;
      //   bool periodicx2 = true;
      //   bool periodicx3 = true;

      //   //bounding box
      //   double g_minX1 = originX1;
      //   double g_minX2 = originX2;
      //   double g_minX3 = originX3;

      //   double g_maxX1 = originX1 + geoLength[0];
      //   double g_maxX2 = originX2 + geoLength[1];
      //   double g_maxX3 = originX3 + geoLength[2];;


      //   //set grid
      //   grid->setDeltaX(cdx);
      //   grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
      //   grid->setPeriodicX1(periodicx1);
      //   grid->setPeriodicX2(periodicx2);
      //   grid->setPeriodicX3(periodicx3);

      //   GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      //   gridCube->setCenterCoordinates(gridCube->getX1Centroid(),meshBand1->getX2Centroid(),gridCube->getX3Centroid());
      //   if(myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());

      //   originX2 = gridCube->getX2Minimum();
      //   g_minX2 = originX2;
      //   g_maxX2 = originX2 + geoLength[1];

      //   GenBlocksGridVisitor genBlocks(gridCube);
      //   grid->accept(genBlocks);

      //   //////////////////////////////////////////////////////////////////////////
      //   if(myid == 0)
      //   {
      //      UBLOG(logINFO, "*****************************************");
      //      UBLOG(logINFO, "* Parameters                            *");
      //      UBLOG(logINFO, "* Re            ="<<Re);
      //      UBLOG(logINFO, "* nuLB          ="<<nuLB);
      //      UBLOG(logINFO, "* uLB           ="<<uLB);
      //      UBLOG(logINFO, "* cdx           ="<<cdx);
      //      UBLOG(logINFO, "* fdx           ="<<fdx);
      //      double Hzb = 0.6/fdx;
      //      UBLOG(logINFO, "* Height of Zackenband ="<<Hzb);
      //      UBLOG(logINFO, "* Re on Zackenband ="<<(uLB*Hzb)/(nuLB*double(1<<refineLevel)));
      //      UBLOG(logINFO, "* nx1/2/3       ="<<nx[0]<<"/"<<nx[1]<<"/"<<nx[2]);
      //      UBLOG(logINFO, "* blocknx1/2/3  ="<<blocknx[0]<<"/"<<blocknx[1]<<"/"<<blocknx[2]);
      //      UBLOG(logINFO, "* x1Periodic    ="<<periodicx1);
      //      UBLOG(logINFO, "* x2Periodic    ="<<periodicx2);
      //      UBLOG(logINFO, "* x3Periodic    ="<<periodicx3);
      //      UBLOG(logINFO, "* number of levels  ="<<refineLevel+1);
      //      UBLOG(logINFO, "* path          ="<<pathname);

      //      UBLOG(logINFO, "*****************************************");
      //      UBLOG(logINFO, "* number of threads    ="<<numOfThreads);
      //      UBLOG(logINFO, "* number of processes  ="<<comm->getNumberOfProcesses());
      //      UBLOG(logINFO, "*****************************************");
      //      UBLOG(logINFO, "*****************************************");     
      //   }
      //   //////////////////////////////////////////////////////////////////////////


      //   //////////////////////////////////////////////////////////////////////////
      //   //refinement
      //   GbCuboid3DPtr refinePlatteBox(new GbCuboid3D(plate->getX1Minimum(), plate->getX2Minimum(), plate->getX3Minimum()+(plate->getX3Maximum()-plate->getX3Minimum())/2.0, 
      //      plate->getX1Maximum()+40.0, plate->getX2Maximum(), plate->getX3Maximum()));
      //   if(myid == 0) GbSystem3D::writeGeoObject( refinePlatteBox.get(), pathname+"/geo/refinePlatteBox", WbWriterVtkXmlASCII::getInstance() );

      //   if (refineLevel > 0)
      //   {
      //      if(myid == 0) UBLOG(logINFO,"Refinement - start");	
      //      RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
      //      refineHelper.addGbObject(refinePlatteBox, refineLevel);
      //      refineHelper.refine();
      //      if(myid == 0) UBLOG(logINFO,"Refinement - end");	
      //   }

      //   /////////////////////////////////////////////////
      //   ///interactoren
      //   int bbOption1 = 1; //0=simple Bounce Back, 1=quadr. BB
      //   D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption1));
      //   D3Q27TriFaceMeshInteractorPtr triPlateInteractor( new D3Q27TriFaceMeshInteractor(plate, grid, bcObst,Interactor3D::SOLID));
      //   D3Q27TriFaceMeshInteractorPtr triBand1Interactor( new D3Q27TriFaceMeshInteractor( meshBand1, grid, bcObst,Interactor3D::SOLID, Interactor3D::EDGES) );
      //   D3Q27TriFaceMeshInteractorPtr triBand2Interactor( new D3Q27TriFaceMeshInteractor( meshBand2, grid, bcObst,Interactor3D::SOLID, Interactor3D::EDGES) );
      //   D3Q27TriFaceMeshInteractorPtr triBand3Interactor( new D3Q27TriFaceMeshInteractor( meshBand3, grid, bcObst,Interactor3D::SOLID, Interactor3D::EDGES) );
      //   D3Q27TriFaceMeshInteractorPtr triBand4Interactor( new D3Q27TriFaceMeshInteractor( meshBand4, grid, bcObst,Interactor3D::SOLID, Interactor3D::EDGES) );

      //   //inflow
      //   GbCuboid3DPtr velBCCuboid(new GbCuboid3D(originX1-blockLengthx1, originX2-blockLengthx1, originX3-blockLengthx1, 
      //      originX1, originX2+geoLength[1]+blockLengthx1, originX3+geoLength[2]+blockLengthx1));
      //   if(myid == 0) GbSystem3D::writeGeoObject(velBCCuboid.get(), pathname+"/geo/velBCCuboid", WbWriterVtkXmlASCII::getInstance());
      //   D3Q27InteractorPtr velBCInteractor(new D3Q27Interactor(velBCCuboid,grid,Interactor3D::SOLID)); 

      //   //inflow
      //   double raiseVelSteps = 0;
      //   vector<D3Q27BCFunction> velcX1BCs,dummy;

      //   mu::Parser inflowProfile;
      //   inflowProfile.SetExpr("uLB"); 
      //   inflowProfile.DefineConst("uLB",uLB);
      //   velcX1BCs.push_back(D3Q27BCFunction(inflowProfile,raiseVelSteps,D3Q27BCFunction::INFCONST));

      //   D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (velcX1BCs,dummy,dummy));
      //   velBCInteractor->addBCAdapter(velBCAdapter);

      //   //outflow
      //   GbCuboid3DPtr densCuboid(new GbCuboid3D(originX1+geoLength[0], originX2-blockLengthx1, originX3-blockLengthx1, 
      //      originX1+geoLength[0]+blockLengthx1, originX2+geoLength[1]+blockLengthx1, originX3+geoLength[2]+blockLengthx1 ));
      //   if(myid == 0) GbSystem3D::writeGeoObject(densCuboid.get(), pathname+"/geo/densCuboid", WbWriterVtkXmlASCII::getInstance());
      //   D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
      //   D3Q27InteractorPtr densInteractor( new D3Q27Interactor(densCuboid,grid,denBCAdapter,Interactor3D::SOLID) );

      //   ////////////////////////////////////////////
      //   //METIS
      //   Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));	

      //   ////////////////////////////////////////////
      //   /////delete solid blocks
      //   if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
      //   InteractorsHelper intHelper(grid, metisVisitor);
      //   intHelper.addInteractor(triPlateInteractor);
      //   intHelper.addInteractor(triBand1Interactor);
      //   intHelper.addInteractor(triBand2Interactor);
      //   intHelper.addInteractor(triBand3Interactor);
      //   intHelper.addInteractor(triBand4Interactor);
      //   intHelper.addInteractor(densInteractor);
      //   intHelper.addInteractor(velBCInteractor);
      //   intHelper.selectBlocks();
      //   if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	 
      //   //////////////////////////////////////

      //   //domain decomposition for threads
      //   if(numOfThreads > 1)
      //   {
      //      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      //      grid->accept(pqPartVisitor);
      //   }

      //   if(myid == 0)
      //   {
      //      UBLOG(logINFO,"Write blocks - start");
      //      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
      //      ppblocks->update(0);
      //      UBLOG(logINFO,"Write blocks - end");
      //   }

      //   unsigned long nob = grid->getNumberOfBlocks();
      //   unsigned long nod = nob * blocknx[0]*blocknx[1]*blocknx[2];
      //   unsigned long nod_real = nob * (blocknx[0]+3)*(blocknx[1]+3)*(blocknx[2]+3);
      //   unsigned long nodb = (blocknx[0]) * (blocknx[1]) * (blocknx[2]);

      //   double needMemAll  = double(nod_real*(27*sizeof(double) + sizeof(int)));
      //   double needMem  = needMemAll / double(comm->getNumberOfProcesses());
      //   
      //   double nup = 0; 

      //   if(myid == 0)
      //   {
      //      UBLOG(logINFO,"Number of blocks = " << nob);
      //      UBLOG(logINFO,"Number of nodes  = " << nod);
      //      int minInitLevel = grid->getCoarsestInitializedLevel();
      //      int maxInitLevel = grid->getFinestInitializedLevel();
      //      for(int level = minInitLevel; level<=maxInitLevel; level++)
      //      {
      //         int nobl = grid->getNumberOfBlocks(level);
      //         UBLOG(logINFO,"Number of blocks for level " << level <<" = " << nobl);
      //         UBLOG(logINFO,"Number of nodes for level " << level <<" = " << nobl*nodb);
      //         nup += nobl*nodb*double(1<<level); 
      //      }
      //      UBLOG(logINFO,"Hypothetically time for calculation step for 120 nodes  = " << nup/6.0e5/(120*8)  << " s");
      //      UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
      //      UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
      //      UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
      //      UBLOG(logINFO,"Available memory per node/8.0 = " << (availMem/8.0) << " bytes");
      //   }

      //   //////////////////////////////////////////
      //   //set connectors
      //   if(myid == 0) UBLOG(logINFO,"set connectors - start");
      //   D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
      //   D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
      //   grid->accept( setConnsVisitor );
      //   if(myid == 0) UBLOG(logINFO,"set connectors - end");

      //   ////////////////////////////
      //   LBMKernel3DPtr kernel;
      //   //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));

      //   //with sponge layer
      //   kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLBWithSpongeLayer(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));
      //   kernel->setWithSpongeLayer(true);
      //   kernel->setSpongeLayer(spongeLayer);

      //   BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      //   kernel->setBCProcessor(bcProc);
      //   SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
      //   grid->accept(kernelVisitor);
      //   //////////////////////////////////
      //   //undef nodes
      //   if (refineLevel > 0)
      //   {
      //      D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
      //      grid->accept(undefNodesVisitor);
      //   }


      //   intHelper.setBC();

      //   //initialization of decompositions
      //   D3Q27ETInitDistributionsBlockVisitor initVisitor( nuLB,rhoLB);
      //   initVisitor.setVx1(uLB);
      //   grid->accept(initVisitor);

      //   //Postprozess
      //   UbSchedulerPtr geoSch(new UbScheduler(1));
      //   D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
      //      new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), 
      //      unitConverter, true));
      //   ppgeo->update(0);
      //   ppgeo.reset();
      //   geoSch.reset();

      //   if(myid == 0) UBLOG(logINFO,"Preprozess - end");      
      //}
      //else
      //{
      //   //domain decomposition for threads
      //   if(numOfThreads > 1)
      //   {
      //      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      //      grid->accept(pqPartVisitor);
      //   }
      //   //set connectors
      //   //grid->setPeriodicX3(false);
      //   D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
      //   D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
      //   grid->accept( setConnsVisitor );
      //   SetSpongeLayerBlockVisitor ssp(spongeLayer);
      //   grid->accept(ssp);

      //   //////////////////////////////////////////////////////////////////////////
      //   //////////////////////////////////////////////////////////////////////////
      //   //Platte
      //   GbTriFaceMesh3DPtr plate (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(PlatteFilename,"Netz"));
      //   plate->rotate(90.0,0.0,0.0);  //TriFacMesh-KO-System anders als LB-KO-System
      //   if(myid == 0) GbSystem3D::writeGeoObject( plate.get(), pathname+"/geo/platte", WbWriterVtkXmlBinary::getInstance() );
      //   //////////////////////////////////////////////////////////////////////////
      //   //////////////////////////////////////////////////////////////////////////
      //   // Zackenband
      //   //////////////////////////////////////////////////////////////////////////
      //   GbTriFaceMesh3DPtr meshBand1 (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "NetzBand"));
      //   meshBand1->translate(-495, -700, -19.94);
      //   if(myid == 0) GbSystem3D::writeGeoObject( meshBand1.get(), pathname+"/geo/Band1", WbWriterVtkXmlASCII::getInstance() );

      //   double blockLengthx1 = blocknx[0]*cdx; //geowerte
      //   double blockLengthx2 = blockLengthx1;
      //   double blockLengthx3 = blockLengthx1;

      //   double geoLength[]   = {  nx[0]*blockLengthx1, nx[1]*blockLengthx2, nx[2]*blockLengthx3}; 

      //   double originX1 = plate->getX1Minimum()-plate->getLengthX1()/4.0;
      //   double originX2 = plate->getX2Minimum();
      //   double originX3 = plate->getX3Minimum()-299.5;

      //   //bounding box
      //   double g_minX1 = originX1;
      //   double g_minX2 = originX2;
      //   double g_minX3 = originX3;

      //   double g_maxX1 = originX1 + geoLength[0];
      //   double g_maxX2 = originX2 + geoLength[1];
      //   double g_maxX3 = originX3 + geoLength[2];;

      //   GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      //   gridCube->setCenterCoordinates(gridCube->getX1Centroid(),meshBand1->getX2Centroid(),gridCube->getX3Centroid());
      //   if(myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());

      //   originX2 = gridCube->getX2Minimum();
      //   g_minX2 = originX2;
      //   g_maxX2 = originX2 + geoLength[1];
      //   //walls
      //   GbCuboid3DPtr addWallZmin (new GbCuboid3D(g_minX1-blockLengthx1, g_minX2-blockLengthx1, g_minX3-blockLengthx1, g_maxX1+blockLengthx1, g_maxX2+blockLengthx1, g_minX3));
      //   if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

      //   GbCuboid3DPtr addWallZmax (new GbCuboid3D(g_minX1-blockLengthx1, g_minX2-blockLengthx1, g_maxX3, g_maxX1+blockLengthx1, g_maxX2+blockLengthx1, g_maxX3+blockLengthx1));
      //   if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

      //   //walls
      //   int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB        
      //   D3Q27BoundaryConditionAdapterPtr slip(new D3Q27SlipBCAdapter(bbOption));
      //   D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, slip,Interactor3D::SOLID));
      //   D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, slip,Interactor3D::SOLID));

      //   SetSolidOrTransBlockVisitor v1(addWallZminInt, SetSolidOrTransBlockVisitor::TRANS);
      //   grid->accept(v1);
      //   addWallZminInt->initInteractor();
      //   SetSolidOrTransBlockVisitor v2(addWallZmaxInt, SetSolidOrTransBlockVisitor::TRANS);
      //   grid->accept(v2);        
      //   addWallZmaxInt->initInteractor();

      //   UbSchedulerPtr geoSch(new UbScheduler(1));
      //   D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
      //      new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), 
      //      unitConverter, true));
      //   ppgeo->update(0);
      //   ppgeo.reset();
      //   geoSch.reset();
      //   //////////////////////////////////////////////////////////////////////////

      //   if(myid == 0) UBLOG(logINFO,"Restart - end"); 
      //}
      //UbSchedulerPtr visSch(new UbScheduler());
      ////visSch->addSchedule(1,0,3);
      ////visSch->addSchedule(100,100,1000);
      ////visSch->addSchedule(1000,1000,5000);
      ////visSch->addSchedule(5000,5000,100000);
      ////visSch->addSchedule(100000,100000,10000000);

      //visSch->addSchedule(1000,1000,10000000);

      //D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), unitConverter);

      //UbSchedulerPtr resSchRMS(new UbScheduler());
      //resSchRMS->addSchedule(1000000,93000,10000000);
      //UbSchedulerPtr resSchMeans(new UbScheduler());
      //resSchMeans->addSchedule(1000000,93000,10000000);
      //UbSchedulerPtr stepAvSch(new UbScheduler());
      //int averageInterval=100;
      //stepAvSch->addSchedule(averageInterval,0,10000000);
      //AverageValuesPostprocessor Avpp(grid, pathname, WbWriterVtkXmlBinary::getInstance(), visSch/*wann wird rausgeschrieben*/, stepAvSch/*wann wird gemittelt*/, resSchMeans,resSchRMS/*wann wird resettet*/);

      //UbSchedulerPtr nupsSch(new UbScheduler(10, 10, 30));
      //nupsSch->addSchedule(500,500,1e6);
      //NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      //UbSchedulerPtr emSch(new UbScheduler(10));
      //EmergencyExitPostprocessor empr(grid, emSch, pathname, RestartPostprocessorPtr(&rp), comm);

      //if(myid == 0)
      //{
      //   UBLOG(logINFO,"PID = " << myid << " Total Physical Memory (RAM): " << MemoryUtil::getTotalPhysMem());
      //   UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used: " << MemoryUtil::getPhysMemUsed());
      //   UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used by current process: " << MemoryUtil::getPhysMemUsedByMe());
      //}

      //string lastStep = string(cstr2);
      //double endTime = UbSystem::stringTo<double>(lastStep);
      //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      //if(myid == 0) UBLOG(logINFO,"Simulation-start");
      //calculation->calculate();
      //if(myid == 0) UBLOG(logINFO,"Simulation-end");
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




#include <iostream>
#include <string>
#include <math.h> 

#include <vfluids.h>

using namespace std;


void run(const char *cstr, double endTime)
{
   try
   {
      string pathname; 
      string pathGeo;
      string pathLog;
      string PlatteFilename;
      string ZckbndFilename;
      int numOfThreads = 1;
      bool logfile = false;
      stringstream logFilename;
      double availMem = 0;

      //UbLog::reportingLevel() = logDEBUG5;

      CommunicatorPtr comm = vf::mpi::MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      string machine = string(cstr);

      if(machine == "my") 
      {
         pathname = "d:/temp/plate";
         pathGeo = "d:/Data/plate";
         pathLog = "d:/temp/plate";
         numOfThreads = 6;
         logfile = false;
         availMem = 15.0e9;
      }
      else if(machine == "Ludwig")      
      {
         pathname = "/work/koskuche/SFB880/plateR1e06";
         pathGeo = "/home/koskuche/data/plate";
         pathLog = "/work/koskuche/SFB880/plateR1e06";
         numOfThreads = 1;
         availMem = 1.0e9;
         logfile = true;
      }
      else if(machine == "Hermit")      
      {
         //Hermit
         pathname = "/univ_1/ws1/ws/xrmkuchr-plate3-0";
         pathGeo = "/zhome/academic/HLRS/xrm/xrmkuchr/data/plate";
         pathLog = "/zhome/academic/HLRS/xrm/xrmkuchr/work/plate";
         numOfThreads = 16;
         availMem = 2.0e9;
         logfile = true;
      }
      else if(machine == "HLRN")      
      {
         //Hermit
         pathname = "/gfs1/work/niivfcpu/scratch/plate";
         pathGeo = "/gfs1/work/niivfcpu/data/plate";
         pathLog = "/gfs1/work/niivfcpu/scratch/plate";
         numOfThreads = 24;
         availMem = 12.0e9;
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
         logFilename <<  pathLog + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(myid)+".txt";
      }

      if(myid ==0 && logfile)
      {
         UbLog::output_policy::setStream(logFilename.str());
      }

      if(myid==0) UBLOG(logINFO,"Testcase plate");

      PlatteFilename = pathGeo + "/platte_raw.stl"; 
      ZckbndFilename= pathGeo + "/2zackenbaender0.stl";

      int baseLevel, refineLevel,nx[3],blocknx[3];
      double Re,velocity,rhoInit,vx1Init;

      //////////////////////////////////////////////////////////////////////////
      //physik
      //////////////////////////////////////////////////////////////////////////
      Re            = 1e6; //11900;// 13286;//13286;//gemessen 18.98 m/s...*spaeter koorestur michael moessner 17m/s
      velocity      = 0.1;  
      vx1Init       = 0.1;  
      rhoInit       = 0.0;

      //int H=200;//200;//392;
      ///////////////Knotenabmessungen:
      nx[0]      = 50;//240;//120;//60;//86;//43;//65;//50;  //l�nge
      nx[1]      = 1;//2;//6;///1;//5;// //breite
      nx[2]      = 16;//64;//32;//18;//5;//15;//15; //h�he gebiet
      blocknx[0] = 25;//10;//6;
      blocknx[1] = 25;//10;//6;
      blocknx[2] = 25;//10;//6;

      baseLevel   = 0;
      refineLevel = 4;

      ///////////////Weltabmessungen:
      double kanalhoeheSI  = 60.0/100.0;//60.0/100.0;//cm, Kanalh�he
      double kanalbreiteSI = kanalhoeheSI*((double)nx[1])/((double)nx[2]);//=kanalh�he*nx1/nx2//1.65/100.0;//13.2/100.0;////40.0/100.0; //cm, Kanalbreite //13.2 zeilbreite
      double kanallaengeSI = kanalhoeheSI*((double)nx[0])/((double)nx[2]);//80.0/100.0;//cm, Kanall�nge, ist nicht angegeben

      // double refinewidth1=kanalhoeheSI/10.0;


      double fineNodeDx   = (kanalhoeheSI) / (double)( blocknx[2]*nx[2]*(1<<refineLevel)+1 ); //+1--> gitter liegt jeweils 0.5dx innerhalb
      //double fineNodeDx   = hReal/100.0;
      double coarseNodeDx = fineNodeDx * (double)(1<<refineLevel);//geowerte

      double blockLengthx1 = blocknx[0]*coarseNodeDx; //geowerte
      double blockLengthx2 = blockLengthx1;
      double blockLengthx3 = blockLengthx1;

      double originX1 = 0.0;//-50.0*propellerDurchmesser;  //geowerte
      double originX2 = 0.0;//-0.5*blockLengthx2*nx2;
      double originX3 = 0.0;// minX3 + 0.5*fineNodeDx;

      double geoLength[]   = {  nx[0]*blockLengthx1, nx[1]*blockLengthx2, nx[2]*blockLengthx3}; 

      //position vorderkante cube
      double originBridgeX1 = 20.0/100.0; //cm, geraten
      double originBridgeX2 = 0.0;//0.5*params.nx[1]*blockLengthx2-0.5*H-fineNodeDx;
      double originBridgeX3 = kanalhoeheSI*0.5;//H*0.0-fineNodeDx; //boden

      bool periodicx1 = false;
      bool periodicx2 = true;
      bool periodicx3 = true;

      //##########################################################################
      //## physical parameters
      //##########################################################################
      double rhoLB         = rhoInit;
      double rhoReal       = 1.0;
      double nuReal  = 0.000015;//0.015;

      double hReal         = 0.0105;//<-m     1.05;//Plattendicke in cm(! cm nicht m !)
      double uReal         = 15;//m/s   //Re*nueReal/hReal;
      double lReal         = 1; //m Plattenl�nge

      //##Machzahl:
      //#Ma     = uReal/csReal
      double Ma      = 0.05;//0.0553;//Ma-Real!
      double csReal  = 343; //uReal/Ma;
      double hLB     = hReal/coarseNodeDx;

      //LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter(hReal, csReal, rhoReal, hLB));
      //LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter(hReal, LBMUnitConverter::AIR_20C, hLB));
      

      double uLB           = 0.1; //uReal   * unitConverter->getFactorVelocityWToLb();
      //double nuLB         = nueReal * unitConverter->getFactorViscosityWToLb();
      double nuLB         = (uLB*(lReal/coarseNodeDx))/Re;
      //double timestep      = unitConverter->getFactorTimeLbToW(coarseNodeDx);

      
      //LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter(0, uReal, uLB, nuReal, nuLB));
      LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter());
      
      velocity = uLB;
      double viscosity = nuLB;

      Grid3DPtr grid(new Grid3D(comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(10000,10000,10000000));
      RestartPostprocessor rp(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY);
      //////////////////////////////////////////////////////////////////////////

      int sizeSP=4;
      mu::Parser spongeLayer;
      //spongeLayer.SetExpr("x1>=(sizeX-sizeSP)/dx ? (sizeX-(x1+1))/sizeSP/2.0 + 0.5 : 1.0");
      spongeLayer.SetExpr("x1>=(sizeX-sizeSP)/dx ? (sizeX-x1)/sizeSP/2.0 + 0.5 : 1.0");
      spongeLayer.DefineConst("sizeX", nx[0]*blocknx[0]);
      spongeLayer.DefineConst("sizeSP", sizeSP*blocknx[0]);

      if (grid->getTimeStep() == 0)
      {
         if(myid==0) UBLOG(logINFO,"Neustart..");
         //bounding box
         double g_minX1 = originX1;
         double g_minX2 = originX2;
         double g_minX3 = originX3;

         double g_maxX1 = originX1 + geoLength[0];
         double g_maxX2 = originX2 + geoLength[1];
         double g_maxX3 = originX3 + geoLength[2];

         //set grid
         grid->setDeltaX(coarseNodeDx);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
         grid->setPeriodicX1(periodicx1);
         grid->setPeriodicX2(periodicx2);
         grid->setPeriodicX3(periodicx3);


         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if(myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         /////////////////////////////////////////////////
         //interactoren definieren
         double geoOverlap = 3.0*coarseNodeDx;

         //inflow
         GbCuboid3DPtr velBCCuboid(new GbCuboid3D(originX1-geoOverlap, originX2-geoOverlap, originX3-geoOverlap, 
            originX1/*+coarseNodeDx*/, originX2+geoLength[1]+geoOverlap, originX3+geoLength[2]+geoOverlap));
         if(myid == 0) GbSystem3D::writeGeoObject(velBCCuboid.get(), pathname+"/geo/velBCCuboid", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr velBCInteractor(new D3Q27Interactor(velBCCuboid,grid,Interactor3D::SOLID)); 

         //inflow
         double uLB2=uLB*0.96*1.02;//*0.5;
         double raiseVelSteps = 0;
         vector<D3Q27BCFunction> velcX1BCs,dummy;

         mu::Parser inflowProfile;
         inflowProfile.SetExpr("uLB"); 

         inflowProfile.DefineConst("uLB",uLB);
         velcX1BCs.push_back(D3Q27BCFunction(inflowProfile,raiseVelSteps,D3Q27BCFunction::INFCONST));

         D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (velcX1BCs,dummy,dummy));
         velBCInteractor->addBCAdapter(velBCAdapter);

         //outflow
         GbCuboid3DPtr densCuboid(new GbCuboid3D(originX1+geoLength[0], originX2-geoOverlap, originX3-geoOverlap, 
            originX1+geoLength[0]+geoOverlap, originX2+geoLength[1]+geoOverlap, originX3+geoLength[2]+geoOverlap ));
         if(myid == 0) GbSystem3D::writeGeoObject(densCuboid.get(), pathname+"/geo/densCuboid", WbWriterVtkXmlASCII::getInstance());
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoInit));
         D3Q27InteractorPtr densInteractor( new D3Q27Interactor(densCuboid,grid,denBCAdapter,Interactor3D::SOLID) );

         //////////////////////////////////////////////////////////////////////////
         if(myid == 0)
         {
            UBLOG(logINFO, "*****************************************");
            UBLOG(logINFO, "* Parameters                            *");
            UBLOG(logINFO, "* Re            ="<<Re);
            UBLOG(logINFO, "* Ma            ="<<Ma);
            UBLOG(logINFO, "* uReal         ="<<uReal);
            UBLOG(logINFO, "* nueReal       ="<<nuReal);
            UBLOG(logINFO, "* nueLB         ="<<nuLB);
            UBLOG(logINFO, "* uLB           ="<<uLB);
            UBLOG(logINFO, "* LX3 (world/LB)="<<kanalhoeheSI<<"/"<<kanalhoeheSI/coarseNodeDx);
            UBLOG(logINFO, "* cdx           ="<<coarseNodeDx);
            UBLOG(logINFO, "* fdx           ="<<fineNodeDx);
            UBLOG(logINFO, "* dx_base       ="<<coarseNodeDx<<" == "<<coarseNodeDx);
            UBLOG(logINFO, "* dx_refine     ="<<fineNodeDx<<" == "<<fineNodeDx );
            UBLOG(logINFO, "* nx1/2/3       ="<<nx[0]<<"/"<<nx[1]<<"/"<<nx[2]);
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
            //UBLOGML(logINFO, "UnitConverter:"<<unitConverter->toString());
            UBLOG(logINFO, "*****************************************");     
         }
         //////////////////////////////////////////////////////////////////////////
         //Platte
         GbTriFaceMesh3DPtr mesh (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(PlatteFilename,"Netz"));

         double x1minMesh = mesh->getX1Minimum(); double x1maxMesh = mesh->getX1Maximum();
         double x2minMesh = mesh->getX2Minimum(); double x2maxMesh = mesh->getX2Maximum();
         double x3minMesh = mesh->getX3Minimum(); double x3maxMesh = mesh->getX3Maximum();

         double drehpunktX=x1minMesh+(x1maxMesh-x1minMesh)*0.5;//triFaceMeshS->getX1Centroid();
         double drehpunktZ=x3minMesh+(x3maxMesh-x3minMesh)*0.5;//triFaceMeshS->getX3Centroid();
         double drehpunktY=x2minMesh+(x2maxMesh-x2minMesh)*0.5;// seedX2-0.5*nodeDelta;//+nx2*deltaX2+0.5*deltaX2;

         mesh->rotate(90.0,0.0,0.0);  //TriFacMesh-KO-System anders als LB-KO-System

         x1minMesh = mesh->getX1Minimum();  x1maxMesh = mesh->getX1Maximum();
         x2minMesh = mesh->getX2Minimum();  x2maxMesh = mesh->getX2Maximum();
         x3minMesh = mesh->getX3Minimum();  x3maxMesh = mesh->getX3Maximum();

         drehpunktX=x1minMesh+(x1maxMesh-x1minMesh)*0.5;//triFaceMeshS->getX1Centroid();
         drehpunktZ=x3minMesh+(x3maxMesh-x3minMesh)*0.5;//triFaceMeshS->getX3Centroid();
         drehpunktY=x2minMesh+(x2maxMesh-x2minMesh)*0.5;// seedX2-0.5*nodeDelta;//+nx2*deltaX2+0.5*deltaX2;

         double H3=1.05/100.0;//cm, Plattendicke
         double scaleB=H3/(x3maxMesh-x3minMesh);
         double scaleX2=(geoLength[2]+2.0*coarseNodeDx)/(x2minMesh-x2maxMesh);

         mesh->scale(scaleB,scaleB,scaleB);
         x1minMesh = mesh->getX1Minimum(); x1maxMesh = mesh->getX1Maximum();
         x2minMesh = mesh->getX2Minimum(); x2maxMesh = mesh->getX2Maximum();
         x3minMesh = mesh->getX3Minimum(); x3maxMesh = mesh->getX3Maximum();
         double offsetXBridge=originBridgeX1;//originBridgeX1;
         double offsetYBridge=originBridgeX2;//originBridgeX2;
         double offsetZBridge=originBridgeX3;//originBridgeX3;//-0.5*(x3minMesh-x3maxMesh);
         //mesh->translate(-x1minMesh+offsetXBridge, -x2minMesh-0.5*offsetYBridge-coarseNodeDx, -x3minMesh+offsetZBridge); 
         mesh->translate(-x1minMesh+offsetXBridge, -x2minMesh+offsetYBridge-coarseNodeDx, -x3minMesh+offsetZBridge-(x3maxMesh-x3minMesh)*0.5/*-hReal*2.0*/); 

         x1minMesh = mesh->getX1Minimum(); x1maxMesh = mesh->getX1Maximum();
         x2minMesh = mesh->getX2Minimum(); x2maxMesh = mesh->getX2Maximum();
         x3minMesh = mesh->getX3Minimum(); x3maxMesh = mesh->getX3Maximum();

         if(myid == 0) GbSystem3D::writeGeoObject( mesh.get(), pathname+"/geo/platte", WbWriterVtkXmlBinary::getInstance() );

         //////////////////////////////////////////////////////////////////////////
         // Zackenband
         //////////////////////////////////////////////////////////////////////////
         GbTriFaceMesh3DPtr meshBand (GbTriFaceMesh3DCreator::readMeshFromFile(ZckbndFilename, "NetzBand"));
         meshBand->deleteRedundantNodes();

         double x1minMeshB = meshBand->getX1Minimum(); double x1maxMeshB = meshBand->getX1Maximum();
         double x2minMeshB = meshBand->getX2Minimum(); double x2maxMeshB = meshBand->getX2Maximum();
         double x3minMeshB = meshBand->getX3Minimum(); double x3maxMeshB = meshBand->getX3Maximum();

         x1minMeshB = meshBand->getX1Minimum();  x1maxMeshB = meshBand->getX1Maximum();
         x2minMeshB = meshBand->getX2Minimum();  x2maxMeshB = meshBand->getX2Maximum();
         x3minMeshB = meshBand->getX3Minimum();  x3maxMeshB = meshBand->getX3Maximum();

         double H1B=1.05/100.0;//*2.0;//0.05;//cm, Banddicke..nachschauen!!!
         double scaleBand=H1B/(x1maxMeshB-x1minMeshB);//H3B/(x3maxMeshB-x3minMeshB);

         meshBand->scale(scaleBand,scaleBand,scaleBand);
         x1minMeshB = meshBand->getX1Minimum(); x1maxMeshB = meshBand->getX1Maximum();
         x2minMeshB = meshBand->getX2Minimum(); x2maxMeshB = meshBand->getX2Maximum();
         x3minMeshB = meshBand->getX3Minimum(); x3maxMeshB = meshBand->getX3Maximum();
         double dBandX=0.5/100.0;//1.29; //15mm-2.1mm Absand von Bandvorderkante
         double dBandY=0.0/100.0;
         double dBandZ=0.223/100.0;//0.344;//....
         double offsetXBridgeB=x1minMesh+dBandX;//originBridgeX1+dBandX;//originBridgeX1;
         double offsetYBridgeB=originBridgeX2+dBandY;//originBridgeX2;
         double offsetZBridgeB=originBridgeX3+dBandZ;//originBridgeX3;//-0.5*(x3minMesh-x3maxMesh);
         meshBand->translate(-x1minMeshB+offsetXBridgeB, -x2minMeshB+offsetYBridgeB-coarseNodeDx, -x3minMeshB+offsetZBridgeB);//-(x3maxMeshB-x3minMeshB)*0.5); 

         x1minMeshB = meshBand->getX1Minimum(); x1maxMeshB = meshBand->getX1Maximum();
         x2minMeshB = meshBand->getX2Minimum(); x2maxMeshB = meshBand->getX2Maximum();
         x3minMeshB = meshBand->getX3Minimum(); x3maxMeshB = meshBand->getX3Maximum();

         GbSystem3D::writeGeoObject( meshBand.get(), pathname+"/geo/Band", WbWriterVtkXmlASCII::getInstance() );

         /////////////////Band2
         GbTriFaceMesh3DPtr meshBand2(GbTriFaceMesh3DCreator::readMeshFromFile(ZckbndFilename, "NetzBand2"));
         meshBand->deleteRedundantNodes();

         double x1minMeshB2 = meshBand2->getX1Minimum(); double x1maxMeshB2 = meshBand2->getX1Maximum();
         double x2minMeshB2 = meshBand2->getX2Minimum(); double x2maxMeshB2 = meshBand2->getX2Maximum();
         double x3minMeshB2 = meshBand2->getX3Minimum(); double x3maxMeshB2 = meshBand2->getX3Maximum();

         x1minMeshB2 = meshBand2->getX1Minimum();  x1maxMeshB2 = meshBand2->getX1Maximum();
         x2minMeshB2 = meshBand2->getX2Minimum();  x2maxMeshB2 = meshBand2->getX2Maximum();
         x3minMeshB2 = meshBand2->getX3Minimum();  x3maxMeshB2 = meshBand2->getX3Maximum();

         double H1B2=1.05/100.0;//0.05;//cm, Banddicke..nachschauen!!!
         double scaleBand2=H1B2/(x1maxMeshB2-x1minMeshB2);//*3.0;//H3B/(x3maxMeshB-x3minMeshB);

         meshBand2->scale(scaleBand2,scaleBand2,scaleBand2);
         x1minMeshB2 = meshBand2->getX1Minimum(); x1maxMeshB2 = meshBand2->getX1Maximum();
         x2minMeshB2 = meshBand2->getX2Minimum(); x2maxMeshB2 = meshBand2->getX2Maximum();
         x3minMeshB2 = meshBand2->getX3Minimum(); x3maxMeshB2 = meshBand2->getX3Maximum();
         double dBandX2=0.5/100.0;//1.29;
         double dBandY2=0.5/100.0;
         double dBandZ2=0.223/100.0;//0.344;//...
         double offsetXBridgeB2=x1minMesh+dBandX2;//originBridgeX1;
         double offsetYBridgeB2=originBridgeX2+dBandY2;//originBridgeX2;
         double offsetZBridgeB2=originBridgeX3+dBandZ2;//originBridgeX3;//-0.5*(x3minMesh-x3maxMesh);
         meshBand2->translate(-x1minMeshB2+offsetXBridgeB2, -x2minMeshB2+offsetYBridgeB2-coarseNodeDx, -x3minMeshB2+offsetZBridgeB2);//-(x3maxMeshB2-x3minMeshB2)*0.5); 

         x1minMeshB2 = meshBand2->getX1Minimum(); x1maxMeshB2 = meshBand2->getX1Maximum();
         x2minMeshB2 = meshBand2->getX2Minimum(); x2maxMeshB2 = meshBand2->getX2Maximum();
         x3minMeshB2 = meshBand2->getX3Minimum(); x3maxMeshB2 = meshBand2->getX3Maximum();

         if(myid == 0) GbSystem3D::writeGeoObject( meshBand2.get(), pathname+"/geo/Band2", WbWriterVtkXmlASCII::getInstance() );
         //////////////////////////////////////////////////////////////////////////

         //////////////////////////////////////////////////////////////////////////
         // refine
         //////////////////////////////////////////////////////////////////////////

         ///////////platte ausmessen:
         x1minMesh = mesh->getX1Minimum(); x1maxMesh = mesh->getX1Maximum();
         x2minMesh = mesh->getX2Minimum(); x2maxMesh = mesh->getX2Maximum();
         x3minMesh = mesh->getX3Minimum(); x3maxMesh = mesh->getX3Maximum();
         double deltaX3Platte=(x3maxMesh-x3minMesh);

         GbCuboid3DPtr refine2PlatteCube(new GbCuboid3D(  originX1-geoOverlap   , originX2-geoOverlap  , x3minMesh-H3*0.5
            , x1maxMesh+H3*5.0, originX2+geoOverlap+geoLength[1], x3maxMesh+H3));
         //if(myid == 0) GbSystem3D::writeGeoObject(refine2PlatteCube.get(), pathname+"/geo/refine2PlatteCube", WbWriterVtkXmlASCII::getInstance());
         //RefineCrossAndInsideGbObjectBlockVisitor refineAdapterP2(refine2PlatteCube, baseLevel, refineLevel-5);
         //grid->accept(refineAdapterP2);

         GbCuboid3DPtr refine3PlatteCube(new GbCuboid3D(   x1minMesh+H3*2.0  , originX2-geoOverlap  , x3minMesh+H3*0.8
            , x1maxMesh-H3*0.2, originX2+geoOverlap+geoLength[1], x3maxMesh+H3*0.1));
         //RefineCrossAndInsideGbObjectBlockVisitor refineAdapterP3(refine3PlatteCube, baseLevel, refineLevel-4);
         //grid->accept(refineAdapterP3);

         GbCuboid3DPtr refine4PlatteCube(new GbCuboid3D(   x1minMesh-H3*2.0  , originX2-geoOverlap  , x3minMesh+deltaX3Platte*0.04
            ,  x1maxMesh+H3*2.0, originX2+geoOverlap+geoLength[1], x3maxMesh+H3*0.25));
         //if(myid == 0) GbSystem3D::writeGeoObject(refine4PlatteCube.get(), pathname+"/geo/refine4PlatteCube", WbWriterVtkXmlASCII::getInstance());
         //RefineCrossAndInsideGbObjectBlockVisitor refineAdapterP4(refine4PlatteCube, baseLevel, refineLevel-3);
         //grid->accept(refineAdapterP4);

         GbCuboid3DPtr refine5PlatteCube(new GbCuboid3D(   originX1-geoOverlap , originX2-geoOverlap  ,x3minMesh-deltaX3Platte/*x3minMesh+deltaX3Platte*0.8*//* x3minMesh+deltaX3Platte*0.8*/
            ,  x1maxMesh+H3*5.0, originX2+geoOverlap+geoLength[1], x3maxMesh+H3));
         //if(myid == 0) GbSystem3D::writeGeoObject(refine5PlatteCube.get(), pathname+"/geo/refine5PlatteCube", WbWriterVtkXmlASCII::getInstance());
         //RefineCrossAndInsideGbObjectBlockVisitor refineAdapterP5(refine5PlatteCube, baseLevel, refineLevel-2);
         //grid->accept(refineAdapterP5);

         GbCuboid3DPtr refine6PlatteCube(new GbCuboid3D(   originX1-geoOverlap   , originX2-geoOverlap  , x3minMesh-deltaX3Platte*3.0/*x3minMesh+deltaX3Platte*0.9*/
            ,  x1maxMesh+H3*7.0, originX2+geoOverlap+geoLength[1], x3maxMesh+deltaX3Platte*3.0));
         //if(myid == 0) GbSystem3D::writeGeoObject(refine6PlatteCube.get(), pathname+"/geo/refine6PlatteCube", WbWriterVtkXmlASCII::getInstance());
         //RefineCrossAndInsideGbObjectBlockVisitor refineAdapterP6(refine6PlatteCube, baseLevel, refineLevel-1);
         //grid->accept(refineAdapterP6);

         //GbCuboid3DPtr wallsX1X2minRef4(new GbCuboid3D(  originX1-3.0*geoOverlap   , originX2-3.0*geoOverlap  , originX1-3.0*geoOverlap
         //  , originX1+geoLength[0]+geoOverlap, originX2+geoOverlap+geoLength[1], kanalhoeheSI*0.1));


         GbCuboid3DPtr refinePlatteBox(new GbCuboid3D(mesh->getX1Minimum(), mesh->getX2Minimum(), mesh->getX3Minimum()+(mesh->getX3Maximum()-mesh->getX3Minimum())/2.0, 
                                                      mesh->getX1Maximum(), mesh->getX2Maximum(), mesh->getX3Maximum()));
         if(myid == 0) GbSystem3D::writeGeoObject( refinePlatteBox.get(), pathname+"/geo/refinePlatteBox", WbWriterVtkXmlASCII::getInstance() );

         /////////////////////////////////////////////////
         ///interactoren
         int bbOption1 = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption1));
         D3Q27TriFaceMeshInteractorPtr triPlateInteractor( new D3Q27TriFaceMeshInteractor(mesh, grid, bcObst,Interactor3D::SOLID));
         D3Q27TriFaceMeshInteractorPtr triBandInteractor( new D3Q27TriFaceMeshInteractor( meshBand, grid, bcObst,Interactor3D::SOLID) );
         D3Q27TriFaceMeshInteractorPtr triBand2Interactor( new D3Q27TriFaceMeshInteractor( meshBand2, grid, bcObst,Interactor3D::SOLID) );

         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));	
        
         //////////////////////////////////////////////////////////////////////////
         //refinement
         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            //refineHelper.addGbObject( refine6PlatteCube, refineLevel-3);
            //refineHelper.addGbObject( refine5PlatteCube, refineLevel-2);
            //refineHelper.addGbObject( refine4PlatteCube, refineLevel-1);
            //refineHelper.addGbObject( refine3PlatteCube, refineLevel);
            refineHelper.addGbObject(refinePlatteBox, refineLevel);
            refineHelper.refine();

            //RefineAroundGbObjectHelper refineHelper(grid, refineLevel, boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triPlateInteractor), 0.0, hReal/4.0);
            //refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         }

         //BlocksPostprocessorPtr ppblocks1(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
         ////if(myid == 0) 
         //ppblocks1->update(0);

         //return;

         //GbCuboid3DPtr testBox(new GbCuboid3D(0.2, -1, 0.1, 1.6, 0.04, 0.5));
         //if(myid == 0) GbSystem3D::writeGeoObject(testBox.get(), pathname+"/geo/testBox", WbWriterVtkXmlASCII::getInstance());
         //D3Q27InteractorPtr testBoxInt(new D3Q27Interactor(testBox, grid, bcObst,Interactor3D::SOLID));

         ////////////////////////////////////////////
         /////delete solid blocks
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(triPlateInteractor);
         intHelper.addInteractor(triBandInteractor);
         intHelper.addInteractor(triBand2Interactor);
         //intHelper.addInteractor(testBoxInt);
         intHelper.addInteractor(densInteractor);
         intHelper.addInteractor(velBCInteractor);
         intHelper.selectBlocks();
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	 
         //////////////////////////////////////

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);


         if(myid == 0) UBLOG(logINFO,"Write blocks - start");
         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
         if(myid == 0) 
            ppblocks->update(0);
         if(myid == 0) UBLOG(logINFO,"Write blocks - end");

         

         unsigned long nob = grid->getNumberOfBlocks();
         unsigned long nod = nob * blocknx[0]*blocknx[1]*blocknx[2];
         unsigned long nod_real = nob * (blocknx[0]+3)*(blocknx[1]+3)*(blocknx[2]+3);

         double needMemAll  = double(nod_real*(27*sizeof(double) + sizeof(int)));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());

         if(myid == 0)
         {
            UBLOG(logINFO,"Number of blocks = " << nob);
            UBLOG(logINFO,"Number of nodes  = " << nod);
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
            UBLOG(logINFO,"Available memory per node/8.0 = " << (availMem/8.0) << " bytes");
         }
         ////////////////////////////
         LBMKernel3DPtr kernel;
         //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));

         //with sponge layer
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLBWithSpongeLayer(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));
         kernel->setWithSpongeLayer(true);
         kernel->setSpongeLayer(spongeLayer);

         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
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
         //////////////////////////////////////////
         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept( setConnsVisitor );

         intHelper.setBC();

         //initialization of decompositions
         D3Q27ETInitDistributionsBlockVisitor initVisitor( nuLB,rhoInit);
         initVisitor.setVx1(vx1Init);
         grid->accept(initVisitor);

         //Postprozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), 
            unitConverter, true));
         ppgeo->update(0);
         //grid->doPostProcess(0);
         ppgeo.reset();
         geoSch.reset();

         if(myid == 0) UBLOG(logINFO,"Preprozess - end");      
         
         //return;
      }
      else
      {
         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept( setConnsVisitor );
         SetSpongeLayerBlockVisitor ssp(spongeLayer);
         grid->accept(ssp);
         if(myid == 0) UBLOG(logINFO,"Restart - end"); 
      }
      UbSchedulerPtr visSch(new UbScheduler());
      //visSch->addSchedule(1,0,3);
      //visSch->addSchedule(100,100,1000);
      //visSch->addSchedule(1000,1000,5000);
      //visSch->addSchedule(5000,5000,100000);
      //visSch->addSchedule(100000,100000,10000000);

      visSch->addSchedule(10000,10000,10000000);
      //visSch->addSchedule(100,100,100000000);

      //UbSchedulerPtr resSchRMS(new UbScheduler());
      //resSchRMS->addSchedule(100000,0,10000000);
      //UbSchedulerPtr resSchMeans(new UbScheduler());
      //resSchMeans->addSchedule(100000,0,10000000);
      //UbSchedulerPtr stepAvSch(new UbScheduler());
      //int averageInterval=1000;
      //stepAvSch->addSchedule(averageInterval,0,10000000);

      //AverageValuesPostprocessor Avpp(grid, pathname + "/steps/stepAV", WbWriterVtkXmlBinary::getInstance(), visSch/*wann wird rausgeschrieben*/, stepAvSch/*wann wird gemittelt*/, resSchMeans,resSchRMS/*wann wird resettet*/);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), unitConverter);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 10, 30));
      nupsSch->addSchedule(1000, 1000, 1000000000);
      NUPSCounterPostprocessor npr(grid, nupsSch, comm);

      //mu::Parser decrViscFunc;
      //decrViscFunc.SetExpr("nue0+c0/(t+1)/(t+1)");
      //decrViscFunc.DefineConst("nue0", nueLB);
      //decrViscFunc.DefineConst("c0", 0.1);
      //UbSchedulerPtr DecrViscSch(new UbScheduler());
      //DecrViscSch->addSchedule(10,10,5000);
      //DecreaseViscosityPostprocessor decrViscPPPtr(grid, DecrViscSch,&decrViscFunc, comm);

      if(myid == 0)
      {
         UBLOG(logINFO,"PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      //double endTime = 80001;
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
int main(int argc, char* argv[])
{

   run(argv[1], UbSystem::stringTo<double>(argv[2]));

   return 0;
}


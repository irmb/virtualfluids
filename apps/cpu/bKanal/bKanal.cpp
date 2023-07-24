

#include <iostream>
#include <string>
#include <math.h> 

#include <sys/types.h> //mkdir rights
#include <sys/stat.h> //mkdir
#include <vfluids.h>

using namespace std;


void run(const char *cstr)
{
   try
   {
      string pathname; 
      string pathnameRestart;
      int numOfThreads =1;
      bool logfile = false;
      stringstream logFilename;
      double availMem = 0;

      UbLog::reportingLevel() = logINFO;

      CommunicatorPtr comm = vf::parallel::MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      string machine = string(cstr);

      if(machine == "my") 
      {
         pathname = "d:/temp/BKanal";
         numOfThreads = 1;
         logfile = false;
         availMem = 10.0e9;
      }
      else if(machine == "Ludwig")      
      {
         pathname =        "/work/koskuche/SFB880/BKanal";
         pathnameRestart = "/work/koskuche/SFB880/BKanal";
         numOfThreads = 8;
         availMem = 1.0e9;
         logfile = true;

         if (myid==0) 
         {
            const char* str = pathname.c_str();
#if defined(__unix__)
            int status=mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif 
         }

         if(myid ==0)
         {
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(myid)+".txt";
         }

      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      if(myid ==0 && logfile)
      {
         UbLog::output_policy::setStream(logFilename.str());
      }

      UBLOG(logINFO,"Testcase BreugemChannel");

      int baseLevel, refineLevel,nx[3],blocknx[3];
      double Re,velocity,rhoInit,vx1Init;//,vx2Init,vx3Init;

      //////////////////////////////////////////////////////////////////////////
      //physik
      //////////////////////////////////////////////////////////////////////////
      Re            = 5500;// 13286;//13286;//gemessen 18.98 m/s...*5.0 zum  testen ob was passiert
      velocity      = 0.01;  
      vx1Init       = 0.01;  
      rhoInit       = 0.0;
      SimulationParametersPtr param = SimulationParameters::getInstanz();

      int H=200;//200;//392;

      //  nx[0]      =8;//ok mit 8// (int)(3.0*(double)H/8.0/8.0);//2;// (int)(0.3*(double)H/6.0/4.0);//das "/4" hier ist wegen der verfeinerung da! //l�nge
      //  nx[1]      =8;//ok mit 8// (int)(2.0*(double)H/8.0/8.0);//2;// (int)(0.2*(double)H/6.0/4.0);//  //breite
      nx[2]      = (int)(2.0*(double)H/5.0/8.0);// //h�he gebiet

      //(3/2/2)-ratio:
      nx[1]=nx[2];
      nx[0]=15;

      blocknx[0] = 15;//10;//5;
      blocknx[1] = 15;//10;//5;
      blocknx[2] = 15;//10;//5;

      baseLevel   = 0;
      refineLevel = 2;//1;////3;//3 soll 1 test;

      ///////////////Weltabmessungen:
      //double kanallaengeSI = ( 2.0*(double)H);
      // double kanalbreiteSI = ( 1.0*(double)H);
      double kanalhoeheSI  = ( 2.0*(double)H);

      // double refinewidth1=kanalhoeheSI/10.0;

      double fineNodeDx   = (kanalhoeheSI) / (double)( blocknx[2]*nx[2]*(1<<refineLevel)+1 ); //+1--> gitter liegt jeweils 0.5dx innerhalb
      double coarseNodeDx = fineNodeDx * (double)(1<<refineLevel);//geowerte

      double blockLengthx1 = blocknx[0]*coarseNodeDx; //geowerte
      double blockLengthx2 = blockLengthx1;
      double blockLengthx3 = blockLengthx1;

      double originX1 = 0.0;//-50.0*propellerDurchmesser;  //geowerte
      double originX2 = 0.0;//-0.5*blockLengthx2*nx2;
      double originX3 = 0.0;// minX3 + 0.5*fineNodeDx;

      double geoLength[]   = {  nx[0]*blockLengthx1, nx[1]*blockLengthx2, nx[2]*blockLengthx3}; 

      bool periodicx1 = true;
      bool periodicx2 = true;
      bool periodicx3 = false;


      //##########################################################################
      //## physical parameters
      //##########################################################################
      double smagorinskiConstant = 0.18;


      double rhoLB         = 0.0;

      double rhoReal       = 1.0;
      double nueReal  = 1;//0.000016;//0.015;

      double hReal         = blocknx[2]*nx[2];//H*0.9;//0.0105;//<-m     1.05;//Plattendicke in cm(! cm nicht m !)
      double uReal         = Re*nueReal/hReal;

      //##Machzahl:
      //#Ma     = uReal/csReal
      double csReal=343.0;
      double Ma      = uReal/csReal;//Ma-Real!
      //double csReal  = uReal/Ma;
      double hLB     = hReal/coarseNodeDx;

      //LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter(hReal, csReal, rhoReal, hLB));

      LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter());
      double uLB          = velocity;
      double nuLB         = (uLB*hLB)/Re;

      //double uLB           = uReal   * unitConverter->getFactorVelocityWToLb();
      //double nueLB         = nueReal * unitConverter->getFactorViscosityWToLb();
      //double timestep      = unitConverter->getFactorTimeLbToW(coarseNodeDx);

      //velocity = uLB;
      // double viscosity =nueLB*1000.0;

      Grid3DPtr grid(new Grid3D(comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(10000,10000,10000000));
      RestartPostprocessor rp(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY);
      grid = rp.restart(-1);
      //////////////////////////////////////////////////////////////////////////

       if (grid->getTimeStep() == 0)
       {
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



         //bottom and top solid bc
         //iteractors
         int bbOption1 = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption1));
         double geoOverlap = coarseNodeDx;
         GbCuboid3DPtr bottomBCCuboid(new GbCuboid3D(originX1-geoOverlap, originX2-geoOverlap, originX3-geoOverlap, 
            originX1+geoLength[0]+coarseNodeDx, originX2+geoLength[1]+geoOverlap, originX3));
         if(myid == 0) GbSystem3D::writeGeoObject(bottomBCCuboid.get(), pathname+"/geo/bottomBCCuboid", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr bottomBCInteractor(new D3Q27Interactor(bottomBCCuboid,grid,bcObst,Interactor3D::SOLID)); 

         GbCuboid3DPtr topBCCuboid(new GbCuboid3D(originX1-geoLength[0]-coarseNodeDx, originX2-geoOverlap, originX3+geoLength[2],//-coarseNodeDx*0.5, 
            originX1+geoLength[0]+coarseNodeDx, originX2+geoLength[1]+geoOverlap, originX3+geoLength[2]+geoOverlap));
         if(myid == 0) GbSystem3D::writeGeoObject(topBCCuboid.get(), pathname+"/geo/topBCCuboid", WbWriterVtkXmlASCII::getInstance());
         D3Q27InteractorPtr topBCInteractor(new D3Q27Interactor(topBCCuboid,grid,bcObst,Interactor3D::SOLID)); 
         //grid->addAndInitInteractor( bottomBCInteractor ); 
         // grid->addAndInitInteractor( topBCInteractor ); 
         //////////////////////////////////////////////////////////////////////////
         if(myid == 0)
         {
            UBLOG(logINFO, "*****************************************");
            UBLOG(logINFO, "* Parameters                            *");
            UBLOG(logINFO, "* Re            ="<<Re);
            UBLOG(logINFO, "* Ma            ="<<Ma);
            UBLOG(logINFO, "* uReal         ="<<uReal);
            UBLOG(logINFO, "* nueReal       ="<<nueReal);
            UBLOG(logINFO, "* nue           ="<<nuLB);
            UBLOG(logINFO, "* velocity      ="<<uLB);
            // UBLOG(logINFO, "* LX1 (world/LB)="<<kanallaengeSI<<"/"<<kanallaengeSI/coarseNodeDx);
            //  UBLOG(logINFO, "* LX2 (world/LB)="<<kanalbreiteSI<<"/"<<kanalbreiteSI/coarseNodeDx);
            UBLOG(logINFO, "* LX3 (world/LB)="<<kanalhoeheSI<<"/"<<kanalhoeheSI/coarseNodeDx);
            UBLOG(logINFO, "* cdx           ="<<coarseNodeDx);
            UBLOG(logINFO, "* fdx           ="<<fineNodeDx);
            UBLOG(logINFO, "* dx_base       ="<<coarseNodeDx<<" == "<<coarseNodeDx);
            UBLOG(logINFO, "* dx_refine     ="<<fineNodeDx<<" == "<<fineNodeDx );
            UBLOG(logINFO, "* nx1/2/3       ="<<nx[0]<<"/"<<nx[1]<<"/"<<nx[2]);
            UBLOG(logINFO, "* blocknx1/2/3  ="<<blocknx[0]<<"/"<<blocknx[1]<<"/"<<blocknx[2]);
            UBLOG(logINFO, "* x2Periodic    ="<<periodicx2);
            UBLOG(logINFO, "* x3Periodic    ="<<periodicx3);
            UBLOG(logINFO, "*****************************************");
/*            UBLOGML(logINFO, "UnitConverter:"<<unitConverter->toString());
            UBLOG(logINFO, "*****************************************");  */   
         }

         if(myid == 0) UBLOG(logINFO,"Refinement - start");	

         //////////////////////////////////////////////////////////////////////////
         // refine
         //////////////////////////////////////////////////////////////////////////
         GbCuboid3DPtr wallsX1X2minRef3(new GbCuboid3D(  originX1-3.0*geoOverlap   , originX2-3.0*geoOverlap  , originX3-3.0*geoOverlap
            , originX1+geoLength[0]+geoOverlap, originX2+geoOverlap+geoLength[1], kanalhoeheSI*0.6/*0.55*/));

         GbCuboid3DPtr wallsX1X2minRef4(new GbCuboid3D(  originX1-3.0*geoOverlap   , originX2-3.0*geoOverlap  ,kanalhoeheSI*0.49
            , originX1+geoLength[0]+geoOverlap, originX2+geoOverlap+geoLength[1], kanalhoeheSI*0.53));

         GbCuboid3DPtr wallsX1X2maxRef2(new GbCuboid3D(  originX1-3.0*geoOverlap   , originX2-3.0*geoOverlap  ,kanalhoeheSI*0.9
            , originX1+geoLength[0]+geoOverlap, originX2+geoOverlap+geoLength[1], originX3+geoOverlap+geoLength[2]));

         GbCuboid3DPtr wallsX1X2maxRef1(new GbCuboid3D(  originX1-3.0*geoOverlap   , originX2-3.0*geoOverlap  ,kanalhoeheSI*0.95
            , originX1+geoLength[0]+geoOverlap, originX2+geoOverlap+geoLength[1], originX3+geoOverlap+geoLength[2]));

         if (refineLevel > 0)
         {

            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            refineHelper.addGbObject(wallsX1X2minRef3, refineLevel-1);
            refineHelper.addGbObject(wallsX1X2minRef4, refineLevel);
            refineHelper.addGbObject(wallsX1X2maxRef2, refineLevel-1);
            refineHelper.addGbObject(wallsX1X2maxRef1, refineLevel);

            refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         }

         ///interactoren
         //int bbOption1 = 0; //0=simple Bounce Back, 1=quadr. BB
         //D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption1));
         ///////w�rfel unten version ende
         ////////////////////////////////////////////////////////////////////////////////
         ////////PM grid
         //Tempor�r:
         //double  H=1.0;

         vector<D3Q27InteractorPtr> D3Q27InteractorPtrarray;
         ////////////////////////////////////////////////////////////////////////////////
         double inflowCubeDx = coarseNodeDx;///(double)(1<<inflowCubeLevel);
         double dpCubes=(double)H/20.0;//100.0; //30zum testen 100real
         double offsetZgrid=H+0.5*inflowCubeDx;
         double epschoch1drittel= 0.928318;
         double abstandIn=2.0*dpCubes;
         double c1oAbst=1.0/abstandIn;
         for (int Nflowdir=0;Nflowdir<((nx[0]*blocknx[0]*c1oAbst)*coarseNodeDx); Nflowdir++)
         {
            // for (int Nhorizon=((nx[2]*blocknx[2]*c1oAbst)*coarseNodeDx)*0.5-2; Nhorizon<((nx[2]*blocknx[2]*c1oAbst)*coarseNodeDx)*0.5-1-1; Nhorizon++)
            // {
            //  for (int Nhorizon=0;  Nhorizon<(((nx[2]*blocknx[2]+1)*c1oAbst)*coarseNodeDx)*0.1; Nhorizon++)//Nhorizon<((nx[2]*blocknx[2]*c1oAbst)*coarseNodeDx)*0.5; Nhorizon++)
            for (int Nhorizon=0;  Nhorizon<(((nx[2]*blocknx[2]+1)*c1oAbst)*coarseNodeDx)*0.5-1; Nhorizon++)//Nhorizon<((nx[2]*blocknx[2]*c1oAbst)*coarseNodeDx)*0.5; Nhorizon++)

            {
               for (int Nspanw=0; Nspanw<((nx[1]*blocknx[1]*c1oAbst)*coarseNodeDx); Nspanw++)
               {
                  // stringstream ss;
                  //     ss<<"cubeH"<<Nflowdir<<"x"<<Nhorizon<<"x"<<Nspanw;
                  ////     //   //geoOrigin ist Mitte, nicht vordere Ecke -> korrigieren
                  // int Nflowdir=1;
                  //int Nhorizon=0;
                  //int Nspanw=1;
                  double xminCubes1=originX1+(Nflowdir*abstandIn)-0.5*dpCubes+0.5*inflowCubeDx+3.0*coarseNodeDx/pow(2.0,refineLevel-1);
                  double xmaxCubes1=originX1+(Nflowdir*abstandIn)+0.5*dpCubes+0.5*inflowCubeDx+3.0*coarseNodeDx/pow(2.0,refineLevel-1);
                  double xminCubes=std::max(xminCubes1,2.0*coarseNodeDx/pow(2.0,refineLevel));
                  double xmaxCubes=std::min(xmaxCubes1,originX1+geoLength[0]-coarseNodeDx/pow(2.0,refineLevel));
                  double yminCubes=std::max(originX2+(Nspanw*abstandIn)-0.5*dpCubes+0.5*inflowCubeDx+3.0*coarseNodeDx/pow(2.0,refineLevel-1),2.0*coarseNodeDx/pow(2.0,refineLevel));
                  double ymaxCubes=std::min(originX2+(Nspanw*abstandIn)+0.5*dpCubes+0.5*inflowCubeDx+3.0*coarseNodeDx/pow(2.0,refineLevel-1),originX2+geoLength[1]-coarseNodeDx/pow(2.0,refineLevel));
                  double zminCubes=std::max(originX3+(Nhorizon*abstandIn)+4.0*coarseNodeDx/pow(2.0,refineLevel-1),2.0*coarseNodeDx/pow(2.0,refineLevel));
                  double zmaxCubes=std::min(originX3+(Nhorizon*abstandIn)+dpCubes+4.0*coarseNodeDx/pow(2.0,refineLevel-1),originX3+geoLength[2]-coarseNodeDx/pow(2.0,refineLevel));
                  ////     /*GbCuboid3D  *rectTemp = new GbCuboid3D(originX1+(Nflowdir*abstandIn)-0.5*dpCubes+0.5*inflowCubeDx, originX2+(Nspanw*abstandIn)-0.5*dpCubes+0.5*inflowCubeDx, originX3+(Nhorizon*abstandIn)-0.5*dpCubes+0.5*inflowCubeDx, 
                  ////										 originX1+(Nflowdir*abstandIn)+0.5*dpCubes+0.5*inflowCubeDx, originX2+(Nspanw*abstandIn)+0.5*dpCubes+0.5*inflowCubeDx, originX3+(Nhorizon*abstandIn)+0.5*dpCubes+0.5*inflowCubeDx );
                  ////*/
                  ////  
                  GbCuboid3DPtr rectTemp(new GbCuboid3D(xminCubes, yminCubes, zminCubes, 
                     xmaxCubes, ymaxCubes, zmaxCubes));
                  ////
                  //     ostringstream ostrcubes;
                  //	 ostrcubes<<pathname <<"/cubeH"<<Nflowdir<<"x"<<Nhorizon<<"x"<<Nspanw;
                  ////       
                  ////   
                  //// // GbSystem3D::writeGeoObject(rectTemp,outpath+cubeschar,WbWriterAvsASCII::getInstance());
                  ////  GbSystem3D::writeGeoObject(rectTemp,ostrcubes.str(),WbWriterAvsASCII::getInstance()); //??
                  //        ostrcubes.str("");
                  //         ostrcubes.clear();

                  ////  boost::shared_ptr<D3Q19AMRInteractor> interactorTemp( new D3Q19AMRInteractor( rectTemp,new D3Q19NoSlipBCAdapter(),AMR3DInteractor::SOLID,ss.str()) );
                  //  //  interactorService.addInteractor(interactorTemp);
                  D3Q27BoundaryConditionAdapterPtr cubeBCAdapter(new D3Q27NoSlipBCAdapter());                   //D3Q27DensityBCAdapter(rhoInit));
                  D3Q27InteractorPtr cubeInteractor( new D3Q27Interactor(rectTemp,grid,cubeBCAdapter,Interactor3D::SOLID));
                  D3Q27InteractorPtrarray.push_back(cubeInteractor);  


               }
            }}
         ////////////////
         //ende cubes
         //////////
         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));

         ////////////////////////////////////////////
         /////delete solid blocks
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(topBCInteractor);
         intHelper.addInteractor(bottomBCInteractor);
         for(size_t i=0; i<D3Q27InteractorPtrarray.size(); ++i)
         {
            intHelper.addInteractor(D3Q27InteractorPtrarray[i]);
         }
         intHelper.selectBlocks();
         if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	 
         //////////////////////////////////////

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

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
         }

         LBMKernel3DPtr kernel;
         kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx[0], blocknx[1], blocknx[2], LBMKernelETD3Q27CCLB::NORMAL));

         // LBMKernel3DPtr kernel(new LBMKernelETD3Q27CascadedTI(blocknx[0], blocknx[1], blocknx[2]));
         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27BGK(blocknx[0], blocknx[1], blocknx[2],1));
         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);
         //	  //scheint neuerdings fuer absturz zu sorgen:
         mu::Parser fctForcingX1;
         fctForcingX1.SetExpr("Fx1*dx");
         fctForcingX1.DefineConst("Fx1", 0.6*5.0e-6);//9.99685e-7);

         kernel->setForcingX1(fctForcingX1);
         kernel->setWithForcing(true); 

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
         intHelper.setBC();

         for(size_t i=0; i<D3Q27InteractorPtrarray.size(); ++i)
         {
            grid->addAndInitInteractor( D3Q27InteractorPtrarray[i] ); 
            char numstr[21];
            sprintf(numstr, "%f", (double)i);
            std::string pathObstCube = pathname+"/geo/obstBCCuboid"+ numstr;
            if(myid == 0) GbSystem3D::writeGeoObject(D3Q27InteractorPtrarray[i]->getGbObject3D().get(),
               /* rectTemp.get(),*/ pathObstCube, WbWriterVtkXmlASCII::getInstance());
         }

         UbTimer timer;
         timer.start();
         grid->accept( metisVisitor );

         if(myid == 0) UBLOG(logINFO,"Write blocks - start");
         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
         if(myid == 0) ppblocks->update(0);
         if(myid == 0) UBLOG(logINFO,"Write blocks - end");

         if(myid == 0) UBLOG(logINFO,"Write blocks - start");
         grid->accept( metisVisitor );
         if(myid == 0) ppblocks->update(1);
         ppblocks.reset();
         if(myid == 0) UBLOG(logINFO,"Write blocks - end");

         //inflow
         double uLB2=uLB;
         double raiseVelSteps = 0;
         vector<D3Q27BCFunction> velcX1BCs,dummy;

         mu::Parser inflowProfile;
         inflowProfile.SetExpr("uLB*0.9"); 

         inflowProfile.DefineConst("uLB",uLB2);
         velcX1BCs.push_back(D3Q27BCFunction(inflowProfile,raiseVelSteps,D3Q27BCFunction::INFCONST));


         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept( setConnsVisitor );

         //domain decomposition

         //initialization of decompositions
         D3Q27ETInitDistributionsBlockVisitor initVisitor( nuLB,rhoInit);
         initVisitor.setVx1(inflowProfile);
         grid->accept(initVisitor);

         //Postrozess
         //LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), 
            unitConverter,true));



         grid->doPostProcess(0);
         ppgeo.reset();
         geoSch.reset();

         if(myid == 0) UBLOG(logINFO,"Preprozess - end");      

      }
            else
            {
               //set forcing
               mu::Parser fctForcingX1, fctForcingX2, fctForcingX3;
               fctForcingX1.SetExpr("Fx1*dx");
               fctForcingX1.DefineConst("Fx1", 0.6*5.0e-6);
               SetForcingBlockVisitor forcingVisitor(fctForcingX1, fctForcingX2, fctForcingX3);
               grid->accept(forcingVisitor);

               //set connectors
               D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
               D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
               grid->accept( setConnsVisitor );
               if(myid == 0) UBLOG(logINFO,"Restart - end"); 
      }


      UbSchedulerPtr visSch(new UbScheduler());
      //visSch->addSchedule(100,1,1000);
      //visSch->addSchedule(1000,1000,10000);
      //visSch->addSchedule(10000,10000,100000);
      //visSch->addSchedule(20000,20000,800000);
      //visSch->addSchedule(50,350000,350500);
      //visSch->addSchedule(50,420000,420500);
      //visSch->addSchedule(50000,420500,10000000);
      visSch->addSchedule(2250,140000,450001);
      UbSchedulerPtr resSch(new UbScheduler());
      resSch->addSchedule(20000,20,10000000);
      UbSchedulerPtr resSchRMS(new UbScheduler());
      resSchRMS->addSchedule(40000,420000,10000000);
      UbSchedulerPtr resSchMeans(new UbScheduler());
      resSchMeans->addSchedule(40000,0,10000000);
      UbSchedulerPtr stepAvSch(new UbScheduler());
      stepAvSch->addSchedule(20,0,10000000);
      AverageValuesPostprocessor Avpp(grid, pathname + "/steps/stepAV", WbWriterVtkXmlBinary::getInstance(), 
                                      visSch/*wann wird rausgeschrieben*/, stepAvSch/*wann wird gemittelt*/, resSchMeans, resSchRMS/*wann wird resettet*/);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), unitConverter, comm);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 90050, 90080));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);


      UbSchedulerPtr AdjForcSch(new UbScheduler());
      AdjForcSch->addSchedule(100,20,20000000);
      D3Q27IntegrateValuesHelperPtr IntValHelp(new D3Q27IntegrateValuesHelper(grid, comm, 
         originX1, originX2, kanalhoeheSI*0.55/*0.501*/, 
         nx[0]*blockLengthx1, nx[1]*blockLengthx2, kanalhoeheSI*0.999));

      double vxZiel=uLB;
      //D3Q27AdjustForcingPostprocessor AdjForcPPPtr(grid, AdjForcSch,IntValHelp, vxZiel*0.6, comm);//da am 11.3.2013 velo um 1/0.6 zu hoch
      D3Q27AdjustForcingPostprocessor AdjForcPPPtr(grid, AdjForcSch,IntValHelp, vxZiel, comm);//dies sollte zu Re=5500 fuehren..

      UbSchedulerPtr visQSch(new UbScheduler());
      visQSch->addSchedule(10,90100,90130);
      QCriterionPostprocessor QKritPtr(grid,pathname+"/steps/Qq",WbWriterVtkXmlBinary::getInstance(),visQSch, comm);

      mu::Parser decrViscFunc;
      decrViscFunc.SetExpr("nue0+c0/(t+1)/(t+1)");
      decrViscFunc.DefineConst("nue0", nuLB);
      decrViscFunc.DefineConst("c0", 0.1);
      UbSchedulerPtr DecrViscSch(new UbScheduler());
      DecrViscSch->addSchedule(10,10,1000);
      DecreaseViscosityPostprocessor decrViscPPPtr(grid, DecrViscSch,&decrViscFunc, comm);

      cout << "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()<<endl;
      cout << "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed()<<endl;
      cout << "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()<<endl;

      double endTime = 2000000;//20000001;
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

   run(argv[1]);

   return 0;
}


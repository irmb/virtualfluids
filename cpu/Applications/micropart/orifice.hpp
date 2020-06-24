#include <vfluids.h>
using namespace std;

void orifice(const char *cstr)
{
   try
   {
      CommunicatorPtr comm(new MPICommunicator());
      int myid = comm->getProcessID();
      int numprocs = comm->getNumberOfProcesses();

      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
      double availMem = 0;
      string geoFile;
      int numOfThreads = 1;

      if(machine == "BOMBADIL") 
      {
    //     pathname = "/work/ehsan/orifice";
		 pathname = "d:/temp/orifice";
         availMem = 6.0e9;
		  int numOfThreads = 1;
         //geoFile = "c:/Data/micropart/DK19_7_02_Martin.stl";
         //geoFile = "c:/Data/micropart/ktoolcav.stl";
         //geoFile = "c:/Data/micropart/boxN.stl";
        geoFile = "d:/Data/Ehsan/orifice.stl";
      }
      else if(machine == "M01" || machine == "M02")      
      {
        // pathname = "/work/koskuche/scratch/mcpart/out";
		  pathname = "/work/ehsan/orifice";
         availMem = 12.0e9;
		  geoFile = "/work/ehsan/data/orifice.stl";
         //geoFile = "/home/koskuche/data/micropart/DK19_7_02_Martin.stl";

         numOfThreads = 1;
         if(myid ==0)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(myid)+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      UbLog::reportingLevel() = logINFO;
      //UbLog::reportingLevel() = logDEBUG1;

      int nodePerBlockX1 =8; //Anzahl an Knoten pro Block
      int nodePerBlockX2 =8;//(int)16;
      int nodePerBlockX3 =8;//8; //(int)16;

      double bH = nodePerBlockX1;    //gewuenschte Rand- und Blockbreite

      //Simulation Parameters
      const int baseLevel = 0;
      const int refineLevel =1;
      //length [m]
      double lSI = 1.55;//223.2;
      //length [LB]
      double lLB = 15;

	  
      double dx =lSI/lLB *2;

      double left_offset = 0;//*0.5;
      double right_offset  = 159;//0.5;//2*0.5
      double front_offset = 750;//0.15;
      double back_offset  = 750;//0.15;
      double top_offset = 250;//0.0;
      double bottom_offset  =750;// 70;//0.07;
	  
	   LBMReal vLB =0.00016103/5.0*sqrt(2.0);//0.00016103;
       LBMReal Re;
       LBMReal rhoLB = 0.0;
       LBMReal nueLB = 0.0000249;//(vLB*lLB)/Re;
       Re = (vLB*(500/dx))/nueLB;
       double dp_Ph=200.0*100000;//
	   double dp_lb=dp_Ph*0.001*(nueLB*dx)*(nueLB*dx);//nue_ph=10e-6 and dx is in micrometer
      // LBMReal nueLB = 0.000016103;
      // LBMReal Re=15000;
      // LBMReal rhoLB = 0.0;
      // LBMReal vLB =nueLB*Re/(500.0/dx);
     // // Re = (vLB*(0.303/dx))/nueLB;
	   // //Re = (vLB*lLB)/nueLB;
	  
      // LBMReal rhoWord = 1e-15;//kg/micrometre^3;//1000.0;
	  // LBMReal nueRE = 1e6;//micromter^2/s;//0.000001;
	  // LBMReal  vWorld=300*1e6;//micrometer/s;//nueRE*Re/ (lSI*4.0/9.0);
	  LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());
      //conv->init(lSI*1e-6,30000,rhoWord,vWorld,lLB,1.0/*rhoLB*/,vLB);
      
	 
 //////////////////////////////////////////////////////////////////////////
      GbObject3DPtr refineCube1(new  GbCuboid3D(78.0,-1.0,-1.0, 81/*370.0*/,20.0/*354.0*/,20.0));//-530.0,-280.0, -72.0, 530.0,354.0,70.0));
      if(myid == 0) GbSystem3D::writeGeoObject(refineCube1.get(), pathname+"/geo/refineCube1", WbWriterVtkXmlASCII::getInstance());

   //   GbObject3DPtr refineCube2(new  GbCuboid3D(-230.0,-90.0, -684.0/*-72.0*/, 600,100.0,70.0));
   //   if(myid == 0) GbSystem3D::writeGeoObject(refineCube2.get(), pathname+"/geo/refineCube2", WbWriterVtkXmlASCII::getInstance());
	  //
	  // GbObject3DPtr refineCube3(new  GbCuboid3D(-350.0,-957.0/*-120.0*/,-684.0/*-684.0*//* -72.0*/, 1700,957.0/*120.0*/,70.0));
   //   if(myid == 0) GbSystem3D::writeGeoObject(refineCube3.get(), pathname+"/geo/refineCube3", WbWriterVtkXmlASCII::getInstance());
	  //
	  // GbObject3DPtr refineCube4(new  GbCuboid3D(-170.0,-60.0, -684.0/*-72.0*/, 200,60.0,70.0));
   //   if(myid == 0) GbSystem3D::writeGeoObject(refineCube4.get(), pathname+"/geo/refineCube4", WbWriterVtkXmlASCII::getInstance());
	  //
	  // GbObject3DPtr refineCubeInlet(new  GbCuboid3D(-10600.0,-600.0, -600.0/*-72.0*/, -9000,600.0,60.0));
   //   if(myid == 0) GbSystem3D::writeGeoObject(refineCubeInlet.get(), pathname+"/geo/refineCubeInlet", WbWriterVtkXmlASCII::getInstance());
	  //
	  //GbObject3DPtr refineCubeOutlet(new  GbCuboid3D(9000,-600.0, -600.0/*-72.0*/,10550.0 ,600.0,60.0));
   //   if(myid == 0) GbSystem3D::writeGeoObject(refineCubeOutlet.get(), pathname+"/geo/refineCubeOutlet", WbWriterVtkXmlASCII::getInstance());
      //////////////////////////////////////////////////////////////////////////
      D3Q27TriFaceMeshInteractorPtr geoInt;
	  /////////////////
      //Grid3DPtr grid(new Grid3D());
        Grid3DPtr grid(new Grid3D(comm));

      UbSchedulerPtr rSch(new UbScheduler());
      rSch->addSchedule(100, 200, 20000);
      RestartPostprocessorPtr rp(new RestartPostprocessor(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY));

      std::string opt;
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());

      if(cstr!= NULL)
         opt = std::string(cstr);

      if/*(cstr== NULL)*/(cstr!= NULL)
      {
         if(myid==0) UBLOG(logINFO,"Restart step: " << opt);
         grid = rp->restart(UbSystem::stringTo<int>(opt));
         rp->reconnect(grid);

         // SetForcingBlockVisitor forcingVisitor(0.0, 0.0, 0.0);
         // grid->accept(forcingVisitor);

         //D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
		  if(myid==0) UBLOG(logINFO,"Restart finish: " << opt);
	 
      }
      else
      {
         if(myid ==0)
         {
            UBLOG(logINFO,"L = " <<lLB );
            UBLOG(logINFO,"v = " <<vLB );
            UBLOG(logINFO,"rho = " <<rhoLB );
            UBLOG(logINFO,"nue = " << nueLB );
			UBLOG(logINFO,"dx = " << dx );
            UBLOG(logINFO,"Re = " << Re );
			 UBLOG(logINFO,"dp_lb = " << dp_lb );
            UBLOG(logINFO,"Preprozess - start");
         }


         ////////////////////////////////////////////////////////////////////////
         //Grid
         //////////////////////////////////////////////////////////////////////////
         grid->setDeltaX(dx);
         grid->setBlockNX(nodePerBlockX1, nodePerBlockX2, nodePerBlockX3);

         ////////////////////////////////////////////////////////////////////////////
         //// Geometrie
         ////////////////////////////////////////////////////////////////////////////
         GbTriFaceMesh3DPtr geo (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(geoFile,"geo"));

         if(myid == 0) GbSystem3D::writeGeoObject(geo.get(), pathname+"/geo/geo", WbWriterVtkXmlASCII::getInstance());

         ////////////////////////////////////////////////////////////////////////////
         //// Randgeometrien erstellen
         ////////////////////////////////////////////////////////////////////////////
         double shiftForMG=grid->getDeltaX(refineLevel)*nodePerBlockX1 / 3.0*2.0;
          GbCuboid3DPtr plate1  = GbCuboid3DPtr( new GbCuboid3D( -7.5, -1.515e-1, -6.831e-2, 7.5, 1.515e-1, 0.0 ));

           GbCuboid3DPtr plate2  = GbCuboid3DPtr( new GbCuboid3D( -1.5e-1, -16.51e-1, -16.831e-2, 1.5e-1, -1.6e-2, 1.0 ));
           GbCuboid3DPtr plate3  = GbCuboid3DPtr( new GbCuboid3D( -1.5e-1, 1.6e-2, -16.831e-2, 1.5e-1, 16.515e-1, 1.0 ));

          // GbCuboid3DPtr plate1_1  = GbCuboid3DPtr( new GbCuboid3D( -7.5, -2.515e-1, -1.0e-1, 7.5, 2.515e-1, -6.831e-2 ));
          // GbCuboid3DPtr plate1_2  = GbCuboid3DPtr( new GbCuboid3D( -7.5, -2.515e-1, -0.0000001, 7.5, 2.515e-1, 1.0e-1 ));
          // GbCuboid3DPtr plate1_3  = GbCuboid3DPtr( new GbCuboid3D( -7.5, 1.515e-1, -6.831e-2, 7.5, 2.515e-1, 0.0  ));
          // GbCuboid3DPtr plate1_4  = GbCuboid3DPtr( new GbCuboid3D( -7.5, -2.515e-1, 0.0, 7.5, -1.515e-1, -1.0e-1 ));

          // GbCuboid3DPtr inflow  = GbCuboid3DPtr( new GbCuboid3D( -8.0, -1.0, -1.0, -7.5, 1.0, 1.0 ));
          // GbCuboid3DPtr outflow = GbCuboid3DPtr( new GbCuboid3D( 7.5, -1.0, -1.0, 8.0, 1.0, 1.0 ));
		  
		   // GbCuboid3DPtr plate2  = GbCuboid3DPtr( new GbCuboid3D( -1.5e-1, -16.51e-1, -16.831e-2, 1.5e-1, -1.6e-2, 1.0 ));
          // GbCuboid3DPtr plate3  = GbCuboid3DPtr( new GbCuboid3D( -1.5e-1, 1.6e-2, -16.831e-2, 1.5e-1, 16.515e-1, 1.0 ));

          // GbCuboid3DPtr plate1_1  = GbCuboid3DPtr( new GbCuboid3D( -left_offset-bH*dx, back_offset, -bottom_offset-bH*dx, right_offset+bH*dx, back_offset+bH*dx, top_offset+bH*dx ));
          // GbCuboid3DPtr plate1_2  = GbCuboid3DPtr( new GbCuboid3D( -left_offset-bH*dx, -front_offset-bH*dx, -bottom_offset-bH*dx, right_offset+bH*dx, -front_offset, top_offset+bH*dx ));
          // GbCuboid3DPtr plate1_3  = GbCuboid3DPtr( new GbCuboid3D( -left_offset-bH*dx, -front_offset-bH*dx, top_offset, right_offset+bH*dx, back_offset+bH*dx, top_offset+bH*dx+2.0*dx ));
          // GbCuboid3DPtr plate1_4  = GbCuboid3DPtr( new GbCuboid3D( -left_offset-bH*dx, -front_offset-bH*dx, -bottom_offset-bH*dx, right_offset+bH*dx, back_offset+bH*dx, -bottom_offset ));

          //GbCuboid3DPtr inflow  = GbCuboid3DPtr( new GbCuboid3D( -left_offset-5*bH*dx, -front_offset-5*bH*dx, -bottom_offset-5*bH*dx, -left_offset, back_offset+5*bH*dx, top_offset+5*bH*dx ));
          //GbCuboid3DPtr outflow = GbCuboid3DPtr( new GbCuboid3D( right_offset, -front_offset-5*bH*dx, -bottom_offset-5*bH*dx, right_offset+5.0*bH*dx, back_offset+5*bH*dx, top_offset+5*bH*dx ));
		  GbCuboid3DPtr inflow  = GbCuboid3DPtr( new GbCuboid3D( -5.0,-1.5, -1.5, 1.5, 20.0, 20.0 ));
		  GbCuboid3DPtr outflow = GbCuboid3DPtr( new GbCuboid3D( 157.50,-1.5, -1.5, 160.5, 20.0, 20.0));


		   GbObject3DPtr gridCube(new GbCuboid3D(inflow->getX1Maximum()-4.0*dx,inflow->getX2Minimum()-4.0*dx ,inflow->getX3Minimum()-4.0*dx,
			   outflow->getX1Minimum()-4.0*dx,outflow->getX2Maximum()-4.0*dx ,outflow->getX3Maximum()-4.0*dx
                                               ));

         GenBlocksGridVisitor genBlocks;
         genBlocks.addGeoObject(gridCube);
         grid->accept(genBlocks);
		  
         if(myid == 0)
         {
            GbSystem3D::writeGeoObject(gridCube.get(),pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());
            //GbSystem3D::writeGeoObject(plate2.get(),pathname+"/geo/plate2", WbWriterVtkXmlASCII::getInstance());
            //GbSystem3D::writeGeoObject(plate3.get(),pathname+"/geo/plate3", WbWriterVtkXmlASCII::getInstance());

            // GbSystem3D::writeGeoObject(plate1_1.get(),pathname+"/geo/plate1_1", WbWriterVtkXmlASCII::getInstance());
            // GbSystem3D::writeGeoObject(plate1_2.get(),pathname+"/geo/plate1_2", WbWriterVtkXmlASCII::getInstance());
            // GbSystem3D::writeGeoObject(plate1_3.get(),pathname+"/geo/plate1_3", WbWriterVtkXmlASCII::getInstance());
            // GbSystem3D::writeGeoObject(plate1_4.get(),pathname+"/geo/plate1_4", WbWriterVtkXmlASCII::getInstance());

            GbSystem3D::writeGeoObject(inflow.get(),pathname+"/geo/inflow", WbWriterVtkXmlASCII::getInstance());
            GbSystem3D::writeGeoObject(outflow.get(),pathname+"/geo/outflow", WbWriterVtkXmlASCII::getInstance());
         }
   

         if (refineLevel > 0)
         {
		  if(myid == 0) UBLOG(logINFO,"Refinement - start");   
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            refineHelper.addGbObject(refineCube1, refineLevel);
    //        refineHelper.addGbObject(refineCube3, 2);
			 //refineHelper.addGbObject(refineCube2, 3);
			 //refineHelper.addGbObject(refineCube4, 4);
			 //
			 //refineHelper.addGbObject(refineCubeInlet, 1);
			 //refineHelper.addGbObject(refineCubeOutlet, 1);
            refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");   
		 
		 
           // RefineCrossAndInsideGbObjectBlockVisitor refVisitor1(refineCube1, refineLevel-4);
            // grid->accept(refVisitor1);

			// RefineCrossAndInsideGbObjectBlockVisitor refVisitor3(refineCube3, refineLevel-3);
            // grid->accept(refVisitor3);

            // RefineCrossAndInsideGbObjectBlockVisitor refVisitor2(refineCube2, refineLevel-2);
            // grid->accept(refVisitor2);
			
			 // RefineCrossAndInsideGbObjectBlockVisitor refVisitor4(refineCube4, refineLevel-1);
            // grid->accept(refVisitor4);

            // RatioBlockVisitor ratioVisitor(refineLevel);
            // grid->accept(ratioVisitor);

            // RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
            // grid->accept(ratioSmoothVisitor);

            // OverlapBlockVisitor overlapVisitor(refineLevel);
            // grid->accept(overlapVisitor);

            // std::vector<int> dirs;
            // D3Q27System::getLBMDirections(dirs);
            // SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
            // grid->accept(interDirsVisitor);
            // if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         }

         //////////////////////////////////////////////////////////////////////////
         //INTERAKTOREN SETZEN (=Randbedingungen)
         //////////////////////////////////////////////////////////////////////////
         //oben/unten = Haftrand
         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
         geoInt = D3Q27TriFaceMeshInteractorPtr( new D3Q27TriFaceMeshInteractor(geo, grid, D3Q27BoundaryConditionAdapterPtr(new D3Q27NoSlipBCAdapter(bbOption)),Interactor3D::INVERSESOLID, Interactor3D::SIMPLE));
	     geoInt->setUseHalfSpaceCheck(true);
         geoInt->setRegardPointInObjectTest(true);
         if(myid == 0) UBLOG(logINFO,"stl - end"); 
         //D3Q27InteractorPtr plate1Int(new D3Q27Interactor(plate1, grid, bcObst,Interactor3D::INVERSESOLID));
         // D3Q27InteractorPtr plate2Int(new D3Q27Interactor(plate2, grid, bcObst,Interactor3D::SOLID));
         // D3Q27InteractorPtr plate3Int(new D3Q27Interactor(plate3, grid, bcObst,Interactor3D::SOLID));

         // D3Q27InteractorPtr plate1_1Int(new D3Q27Interactor(plate1_1, grid, bcObst,Interactor3D::SOLID));
         // D3Q27InteractorPtr plate1_2Int(new D3Q27Interactor(plate1_2, grid, bcObst,Interactor3D::SOLID));
         // D3Q27InteractorPtr plate1_3Int(new D3Q27Interactor(plate1_3, grid, bcObst,Interactor3D::SOLID));
         // D3Q27InteractorPtr plate1_4Int(new D3Q27Interactor(plate1_4, grid, bcObst,Interactor3D::SOLID));

         //links: geschwindigkeits-einfluss
         //Velocity-BC
         //////////////////////////////////////////////////////////////////////////
         mu::Parser fct;
         fct.DefineConst("vx1"  , vLB*9.0/4.0 );
         //fct = MathUtil::getDuctParaboloidX(0, 250*2.0, -51.08/2, 51.08, vLB*9.0/4.0);
         fct.SetExpr("vx1");
         //////////////////////////////////////////////////////////////////////////

         //////////////////////////////////////////////////////////////////////////
             D3Q27BoundaryConditionAdapterPtr velBCAdapter = D3Q27BoundaryConditionAdapterPtr(new D3Q27VelocityBCAdapter (false, false ,true ,fct, 0, D3Q27BCFunction::INFCONST));
    //     D3Q27BoundaryConditionAdapterPtr velBCAdapter = D3Q27BoundaryConditionAdapterPtr(new D3Q27VelocityBCAdapter (false, false ,true ,fct,fct,fct, 0, D3Q27BCFunction::INFCONST));
			 // velBCAdapter->setSecondaryBcOption(2);
            // D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(inflow, grid, velBCAdapter, Interactor3D::SOLID));

		 D3Q27BoundaryConditionAdapterPtr denBCAdapterInlet(new D3Q27DensityBCAdapter(3.0*(dp_lb-rhoLB)));
        denBCAdapterInlet->setSecondaryBcOption(1);
        D3Q27InteractorPtr inflowInt = D3Q27InteractorPtr( new D3Q27Interactor(inflow, grid, denBCAdapterInlet,Interactor3D::SOLID));
		 
         //rechts: druckrand
         //Density-BC
         //fuer Kompressibles Modell  rho = 1.0
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         denBCAdapter->setSecondaryBcOption(1);
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(outflow, grid, denBCAdapter,Interactor3D::SOLID));

         MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B);
         grid->accept( metisVisitor );
         
      //   BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
        // if(myid == 0) ppblocks->update(0);
         
         SolidBlocksHelper sd(grid, comm);
         sd.addInteractor(geoInt);
         sd.addInteractor(inflowInt);
         sd.addInteractor(outflowInt);
         // sd.addInteractor(plate1_1Int);
         // sd.addInteractor(plate1_2Int);
         // sd.addInteractor(plate1_3Int);
         // sd.addInteractor(plate1_4Int);
         // sd.addInteractor(plate2Int);
         // sd.addInteractor(plate3Int);
		   if(myid == 0) UBLOG(logINFO,"line"<<__LINE__); 
         sd.deleteSolidBlocks();     
         if(myid == 0) UBLOG(logINFO,"line"<<__LINE__); 
         grid->accept( metisVisitor );
         if(myid == 0) UBLOG(logINFO,"line"<<__LINE__);

         sd.setTransBlocks();
         if(myid == 0) UBLOG(logINFO,"line"<<__LINE__);
         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
         if(myid == 0) ppblocks->update(0);
         if(myid == 0) ppblocks.reset();

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nod_temp = nob * (nodePerBlockX1+gl) * (nodePerBlockX2+gl) * (nodePerBlockX3+gl);
         unsigned long nod = nob * (nodePerBlockX1) * (nodePerBlockX2) * (nodePerBlockX3);
         double needMemAll  = double(nod_temp*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());


         if(myid == 0)
         {
            UBLOG(logINFO,"Number of blocks = " << nob);
            UBLOG(logINFO,"Number of nodes  = " << nod);
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
         }  

         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(nodePerBlockX1, nodePerBlockX2, nodePerBlockX2));
         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(nodePerBlockX1, nodePerBlockX2, nodePerBlockX2));

		
		  int option = 0;
		// LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(nodePerBlockX1, nodePerBlockX2, nodePerBlockX3,option));
		  LBMKernel3DPtr kernel;
		kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(nodePerBlockX1, nodePerBlockX2, nodePerBlockX3, LBMKernelETD3Q27CCLB::MAGIC));
		//  
		// kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(nodePerBlockX1, nodePerBlockX2, nodePerBlockX3, 1));
		 
		 
         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }
		 
		  if(myid == 0) UBLOG(logINFO,"intractor - start"); 
          //inflow
         grid->addAndInitInteractor(inflowInt);

         //outflow
         grid->addAndInitInteractor(outflowInt);
         //canal
         grid->addAndInitInteractor(geoInt);
         // grid->addAndInitInteractor(plate1_1Int);
         // grid->addAndInitInteractor(plate1_2Int);
         // grid->addAndInitInteractor(plate1_3Int);
         // grid->addAndInitInteractor(plate1_4Int);
         // grid->addAndInitInteractor(plate2Int);
         // grid->addAndInitInteractor(plate3Int);

       

         //////////////////////////////////////////////////////////////////////////
         //connectoren setzen:

         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );
         //////////////////////////////////////////////////////////////////////////	 
         //domain decomposition
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);
			 
         //////////////////////////////////////////////////////////////////////////     
	   //Stroemungsfeld initialisieren
         //////////////////////////////////////////////////////////////////////////
         D3Q27ETInitDistributionsBlockVisitor initVisitor(rhoLB); //1.0
         initVisitor.setVx1(0); 
		   initVisitor.setVx1(0); 

         grid->accept(initVisitor);

         if(myid == 0)
         {
            //Abstände "q" als Linien rausschreiben
            std::vector< UbTupleFloat3 > nodes;
            std::vector< UbTupleInt2 >   lines;
            geoInt->addQsLineSet(nodes, lines);
            WbWriterVtkXmlBinary::getInstance()->writeLines(pathname+"/grid/qs",nodes,lines);
         }

          if(myid == 0) UBLOG(logINFO,"Preprozess - end");
		 
		 	  ////////////////////////
           //Set Postprozessors
           //////////////////////////////////////////////////////////////////////////
           {
            UbSchedulerPtr geoSch(new UbScheduler(1));
            D3Q27MacroscopicQuantitiesPostprocessor ppgeo(grid,geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv,  comm, true);
            grid->doPostProcess(0);
           }
	    

      }

      //////////////////////////////////////////////////////////////////////////
	   // UbSchedulerPtr visSchAv(new UbScheduler());
		UbSchedulerPtr visSchAv(new UbScheduler(100,100));
      // visSchAv->addSchedule(100,10,1000);
      // UbSchedulerPtr resSchAv(new UbScheduler());
	   UbSchedulerPtr resSchAv(new UbScheduler(100,100));
      // resSchAv->addSchedule(20,20,1000);
      AverageValuesPostprocessor       Avpp(grid,  pathname + "/Turbulence/stepAV", WbWriterVtkXmlBinary::getInstance(), visSchAv/*wann wird rausgeschrieben*/,resSchAv/*wann wird resettet*/,comm);
	  
	   D3Q27ShearStressPostprocessor  shear(grid,  pathname + "/shear/step", WbWriterVtkXmlBinary::getInstance(), visSchAv/*wann wird rausgeschrieben*/,resSchAv/*wann wird resettet*/,comm,iProcessor); 
	   //D3Q27ShearStressPostprocessor  shear(grid,  pathname + "/shear/step", WbWriterVtkXmlBinary::getInstance(), visSchAv/*wann wird rausgeschrieben*/,resSchAv/*wann wird resettet*/,comm);
	   shear.addInteractor(geoInt);
	   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  

      UbSchedulerPtr nupsSch(new UbScheduler(1, 5, 10));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      double outTime = 100.0;
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid,stepSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv,  comm);
      //////////////////////////////////////////////////////////////////////////
      //PathLine
       UbSchedulerPtr plSch(new UbScheduler(5000, 5000));
      const int numberofparticle=20;
	
	  std::vector<UbTupleDouble3 > potisions;
	  double randomx[numberofparticle];
	  double randomy[numberofparticle];
	  double randomz[numberofparticle];
	  double lowestx,highestx,lowesty,highesty,lowestz,highestz;
	  if(myid==0)
	  {
		  for(int i = 0; i < numberofparticle; i++)
		  {
			  double random; 
	        lowestx =-10300.0;  lowesty =-230;          lowestz =-250;
	        highestx=-9792.0;  highesty=-330;          highestz=-250; 
		  
	      double rangex=(highestx-lowestx),rangey=(highesty-lowesty),rangez=(highestz-lowestz);	
           randomx[i] = lowestx+(rangex*rand()/(RAND_MAX + 1.0));
		   randomy[i] = lowesty+(rangey*rand()/(RAND_MAX + 1.0));
	       randomz[i] = lowestz+(rangez*rand()/(RAND_MAX + 1.0));
		  //val<1>(potisions[i])= 0.506983973456;
		  //val<2>(potisions[i]) = lowesty+(rangey*rand()/(RAND_MAX + 1.0));
		   //val<3>(potisions[i]) = lowestz+(rangez*rand()/(RAND_MAX + 1.0));
		  }
		  for (int i=0;i<comm->getNumberOfProcesses();i++)
		  {
			  if (i!=0)
			  {
			      MPI_Send(randomx,numberofparticle, MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD);
				  MPI_Send(randomy,numberofparticle, MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD);
				  MPI_Send(randomz,numberofparticle, MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD);
			  }
		  }
	  }
	  if (myid!=0)
	  {
		  MPI_Status status; 
		  MPI_Recv(randomx,numberofparticle, MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		  MPI_Recv(randomy,numberofparticle, MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		  MPI_Recv(randomz,numberofparticle, MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	  }
	  for(int i = 0; i < numberofparticle; i++)
	  {	
		  potisions.push_back( makeUbTuple(randomx[i],randomy[i],randomz[i]) );
		  //val<1>(potisions[i])= 0.506983973456;
		  //val<2>(potisions[i]) = randomy[i];
		  //val<3>(potisions[i]) = randomz[i];
	  }
	 //  UBLOG(logINFO,"Rank="<<myid<<" positions  = " <<val<1>(potisions)<< " "<<val<2>(potisions)<<" "<< val<3>(potisions));
	 // D3Q27InterpolationProcessorPtr iProcessor2;
     // D3Q27PathLinePostprocessorMcpart pathLine(grid, pathname + "/pathLine/pathLine", WbWriterVtkXmlASCII::getInstance(), conv, plSch, comm,potisions, nueLB, iProcessor);
      //////////////////////////////////////////////////////////////////////////
      //Simulation
      //////////////////////////////////////////////////////////////////////////

	  double endTime = 1000.0;
      UbSchedulerPtr visSch(stepSch);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      if(myid == 0) UBLOG(logINFO,"Simulation-start");
      calculation->calculate();
      if(myid == 0) UBLOG(logINFO,"Simulation-end");

   }
   catch(std::exception& e)
   {
      cerr << e.what() << endl;
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
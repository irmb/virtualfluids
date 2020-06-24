#include <iostream>
#include <string>

#include <vfluids.h>
using namespace std;


void run(const char *cstr)
{
   CommunicatorPtr comm(new MPICommunicator());
   try
   {
      //Sleep(30000);
      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
	  string geosphere;	  
      int numOfThreads = 1;
      double availMem = 0;

      
      int myid = comm->getProcessID();

      if(machine == "BOMBADIL") 
      {
         pathname = "d:/temp/shear";
         numOfThreads = 1;
         availMem = 3.0e9;
         geosphere = "d:/Data/Ehsan/Agglomerat_n_00020_fd_1.882437_r_0033.930997.txt";
      }
       else if(machine == "M01" || machine == "M02")      
      {
         pathname = "/work/koskuche/scratch/smallAgg80";
		// geosphere = "/work/ehsan/data/Agglomerat4.txt";
	   // geosphere = "/work/ehsan/data/Agglomerat_n_00060_fd_1.858514_r_0061.500327.txt";
		// geosphere = "/work/ehsan/data/Agglomerat_n_00080_fd_1.855984_r_0071.870085.txt";
		//geosphere = "/work/ehsan/data/Agglomerat_n_00040_fd_1.864231_r_0049.358563.txt";
		//geosphere = "/work/ehsan/data/Agglomerat_n_00020_fd_1.882437_r_0033.930997.txt";
		geosphere = "/work/ehsan/data/Agglomerat_n_00500_fd_1.850643_r_0193.702967.txt";
      
		
		
		
         numOfThreads = 1;
         availMem =1.0e10;// 12.0e9;

         if(myid ==0)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
	   }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      double dx =0.1*4.0;

	  double eq_Diameter=2.0*38.0;//55.3586;//61.5003;//80;//71.8701;//61.5003;
      double L1 =35.0*eq_Diameter;
      double L2, L3, H;
      L2 = L3 = H =35.0*eq_Diameter;//1.0;//0.42*3.9;

      LBMReal radius = 6.0;
      LBMReal rhoReal = 1.0; //kg/m^3
      //LBMReal uReal = 0.45;//m/s
   //   LBMReal uLB = 0.05;
      LBMReal Re = 0.1;
      LBMReal rhoLB = 0.0;
      LBMReal l = L2 / dx;

      //LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter(1.0, 1/sqrt(3.0)*(uReal/uLB), 1.0, 1.0/dx, dx*dx*dx));
      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = 5;  
  
     

      //bounding box
      double d_minX1 = 0.0;
      double d_minX2 = 0.0;
      double d_minX3 = 0.0;

      double d_maxX1 = L1;
      double d_maxX2 = L2;
      double d_maxX3 = L3;

	  const int blocknx1 = 8;
      const int blocknx2 = 8;
      const int blocknx3 = 8;

      dx =14.4*2.0;// ( 0.9) *(double)(1<<refineLevel);

	  double area =/* radius/dx;*/radius*radius*PI/(dx/(double)(1<<refineLevel))/(dx/(double)(1<<refineLevel));//2.0*radius*H;
	  double nueReal=1e-6;//water
	  double uReal=Re*nueReal/(2.0*radius);//real velocity
	  double F_stokss=6*PI*.001/*water*/*radius*uReal;
      //LBMReal nueLB = (((4.0/9.0)*uLB)*2.0*(radius/dx))/Re;
      LBMReal nueLB =.75/((double)(1<<refineLevel));// ((uLB)*2.0*(radius/dx))/Re;
	//  LBMReal uLB  = ((nueLB)*Re)/ ( 2.0*(radius/dx));//base on the coarsest level
	  LBMReal uLB  = ((0.75)*Re)/ ( 2.0*(radius/(dx/((double)(1<<refineLevel)))));//base on the coarsest level if nueLB =.75/((double)(1<<refineLevel)) and  dx = ( 0.9) *(double)(1<<refineLevel)
  //     LBMReal uLB  = ((0.75)*Re)/ ((eq_Diameter/(dx/((double)(1<<refineLevel)))));//base on the coarsest level if nueLB =.75/((double)(1<<refineLevel)) and  dx = ( 0.9) *(double)(1<<refineLevel)
	 double blockLength = blocknx1*dx;
	

		double xsphere=1.0*L1/2.0;//0.5;
		double ysphere=L2/2.0;//0.75;
		double zsphere=L3/2.0;//0.75;
 //obstacle
    	 //////////////////////////////////////////////////////////////////////////
	  UbFileInputASCII file;
	  file.open(geosphere);
	  //file.skipLine();file.skipLine();//2line skiped
	  std::string NOP=file.readString();
	  std::string NOP2=file.readString();
	  const int numberOfParticles=file.readDouble();
	  if(myid == 0){UBLOG(logINFO,__FILE__<<" " <<__LINE__<<" number of particles="<<numberOfParticles);}	
	  //std::string Dia=file.readString();
	  double diameter=2.0*radius;//12;//file.readDouble();
	  file.skipLine();file.skipLine();file.skipLine();file.skipLine();file.skipLine();file.skipLine();file.skipLine();//7 line skiped
			GbSphere3DPtr *sphereP=new GbSphere3DPtr[numberOfParticles];
			  
			for (int i=0;i<numberOfParticles;i++)
			{
			double x=file.readDouble();
			double y=file.readDouble();
			double z=file.readDouble();
		///0degree in x direction		
			    double x_rotation= x;
			    double y_rotation= y;
			    double z_rotation= z;
///180degree in x direction		
			   // double x_rotation= x;
			   // double y_rotation= -y;
			   // double z_rotation= -z;			   
		///90degree in y direction	
			  // double x_rotation=-z;
			  // double y_rotation= y;
			  // double z_rotation=x;			
	   // ///90degree in z axis	
			   // double x_rotation=-y;
			   // double y_rotation=x;
			   // double z_rotation=z;
		//transfer	
			double x_final=x_rotation/*/1450*/  +xsphere;
			double y_final=y_rotation/*/1450*/  +ysphere;
			double z_final=z_rotation/*/1450*/  +zsphere;
				sphereP[i]=GbSphere3DPtr(new GbSphere3D(x_final, y_final, z_final, diameter/2.0/*/1450*/));
				if(myid == 0)GbSystem3D::writeGeoObject(sphereP[i].get(),pathname + "/sphere/sphere"+ "_" + UbSystem::toString(i), WbWriterVtkXmlASCII::getInstance());
			}
			file.close();
///////////////////////////////	
       D3Q27InteractorPtr *spherePInt=new D3Q27InteractorPtr[numberOfParticles];	
      double offs = dx;

      //double g_minX1 = d_minX1-offs-0.499999*dx;
      double g_minX1 = d_minX1-offs;
      double g_minX2 = d_minX2-offs;
      double g_minX3 = d_minX3-offs;

      double g_maxX1 = d_maxX1+offs;
      double g_maxX2 = d_maxX2+offs;
      double g_maxX3 = d_maxX3+offs;
			if(myid == 0){UBLOG(logINFO,__FILE__<<" " <<__LINE__);}	  
      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
			if(myid == 0){UBLOG(logINFO,__FILE__<<" " <<__LINE__);}	  
      

      //refinement area
      double rf = 0.50*blockLength;
      // GbObject3DPtr refineCube(new  GbCuboid3D(sphereP[0]->getX1Minimum()-rf*3.0/4.0, sphereP[3]->getX2Minimum()-rf*3.0/4.0, sphereP[5]->getX3Minimum()-rf*1.0/2.0, 
         // sphereP[2]->getX1Maximum()+rf*3.0/4.0, sphereP[4]->getX2Maximum()+rf*3.0/4.0, sphereP[6]->getX3Maximum()+rf*1.0/2.0));
		 
		 //////////
	   double level5=xsphere-(xsphere-eq_Diameter/2-0.50*(blocknx1*dx/pow(2.0,5)));//0.065;//.085;
	   double level4=level5+.1*(blocknx1*dx/pow(2.0,4));//0.015;//0.1;
	   double level3=level4+0.50*(blocknx1*dx/pow(2.0,3));//0.015;//0.115;
	   double level2=level3+1.0*(blocknx1*dx/pow(2.0,2));//.035;//0.15;
	   double level1=level2+1.0*(blocknx1*dx/pow(2.0,1));//.05;//0.2;
	   
	    GbCuboid3DPtr refineCube1(new GbCuboid3D(  xsphere-level1,ysphere-level1, zsphere-level1,xsphere+level1,ysphere+level1, zsphere+level1));
	    GbCuboid3DPtr refineCube2(new GbCuboid3D(  xsphere-level2,ysphere-level2, zsphere-level2,xsphere+level2,ysphere+level2, zsphere+level2));
		GbCuboid3DPtr refineCube3(new GbCuboid3D(  xsphere-level3,ysphere-level3, zsphere-level3,xsphere+level3,ysphere+level3, zsphere+level3));
		GbCuboid3DPtr refineCube4(new GbCuboid3D(  xsphere-level4,ysphere-level4, zsphere-level4,xsphere+level4,ysphere+level4, zsphere+level4));
		GbCuboid3DPtr refineCube5(new GbCuboid3D(  xsphere-level5,ysphere-level5, zsphere-level5,xsphere+level5,ysphere+level5, zsphere+level5));
		 ///////////

      Grid3DPtr grid(new Grid3D(comm));

      UbSchedulerPtr rSch(new UbScheduler(100000, 100000));
      //RestartPostprocessorPtr rp(new RestartPostprocessor(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY));

      //UbSchedulerPtr emSch(new UbScheduler(1000, 1000));
      //EmergencyExitPostprocessor em(grid, emSch, pathname+"/checkpoints/emex.txt", rp, comm);

      std::string opt;

      if(cstr!= NULL)
         opt = std::string(cstr);

      if/*(cstr== NULL)*/(cstr!= NULL)
      {
         opt = std::string(cstr);

         if(myid==0) UBLOG(logINFO,"Restart step: " << opt);

         //grid = rp->restart(UbSystem::stringTo<int>(opt));
         //rp->reconnect();

         //cylinderInt = 

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );

         //domain decomposition
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);
      }
      else
      {
         if(myid ==0)
         {
            UBLOG(logINFO,"L = " << L2/dx );
            UBLOG(logINFO,"v = " << uLB );
            UBLOG(logINFO,"rho = " << rhoLB );
            UBLOG(logINFO,"nue = " << nueLB );
            UBLOG(logINFO,"Re = " << Re );
			UBLOG(logINFO,"F_stokss = " << F_stokss );
			UBLOG(logINFO,"dx = " << dx );
            UBLOG(logINFO,conv->toString() );
            UBLOG(logINFO,"Preprozess - start");
         }

         grid->setDeltaX(dx);
         grid->setBlockNX(blocknx1, blocknx2, blocknx3);
		 grid->setPeriodicX1(false);
         grid->setPeriodicX2(true);
         grid->setPeriodicX3(true);

         // UbTupleDouble6 bouningBox(gridCube->getX1Minimum(),gridCube->getX2Minimum(),gridCube->getX3Minimum(),
         // gridCube->getX1Maximum(),gridCube->getX2Maximum(),gridCube->getX3Maximum());
         // UbTupleInt3 blockNx(blocknx1, blocknx2, blocknx3);
         // UbTupleInt3 gridNx(8, 16, 16);
         // grid = Grid3DPtr(new Grid3D(bouningBox, blockNx, gridNx));

         if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
 //        if(myid ==0) GbSystem3D::writeGeoObject(refineCube.get(),pathname + "/geo/refineCube", WbWriterVtkXmlBinary::getInstance());
 
         ////
         if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
         if(myid ==0) GbSystem3D::writeGeoObject(refineCube1.get(), pathname + "/geo/refineCube1", WbWriterVtkXmlBinary::getInstance());
		 if(myid ==0) GbSystem3D::writeGeoObject(refineCube2.get(),pathname + "/geo/refineCube2", WbWriterVtkXmlBinary::getInstance());
         if(myid ==0) GbSystem3D::writeGeoObject(refineCube3.get(),pathname + "/geo/refineCube3", WbWriterVtkXmlBinary::getInstance());
		 if(myid ==0) GbSystem3D::writeGeoObject(refineCube4.get(),pathname + "/geo/refineCube4", WbWriterVtkXmlBinary::getInstance());
		 if(myid ==0) GbSystem3D::writeGeoObject(refineCube5.get(),pathname + "/geo/refineCube5", WbWriterVtkXmlBinary::getInstance());
		 ////
		 
         GenBlocksGridVisitor genBlocks;
         genBlocks.addGeoObject(gridCube);
         grid->accept(genBlocks);

         //walls
         GbCuboid3DPtr addWallYmin (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_minX2, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_minX3));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_maxX2, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_maxX3, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_minX1, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow (new GbCuboid3D(d_maxX1, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));

        if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
			RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
			refineHelper.addGbObject(refineCube5, refineLevel);
			
		//	refineHelper.addGbObject(refineCube1, refineLevel);
			// refineHelper.addGbObject(refineCube2, refineLevel-1);
			// refineHelper.addGbObject(refineCube3, refineLevel-2);
			// refineHelper.addGbObject(refineCube4, refineLevel-3);
			//refineHelper.addGbObject(refineCube5, refineLevel-4);
			

            refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");   
		 
         }

         MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B);
         grid->accept( metisVisitor );

         SolidBlocksHelper sd(grid, comm);

         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
		 D3Q27BoundaryConditionAdapterPtr bcObst2(new D3Q27SlipBCAdapter(bbOption));
        // cylinderInt = D3Q27InteractorPtr ( new D3Q27Interactor(cylinder, grid, bcObst,Interactor3D::SOLID));
		
			for (int i=0;i<numberOfParticles;i++)
			{      
				spherePInt[i]= D3Q27InteractorPtr( new D3Q27Interactor(sphereP[i], grid, bcObst,Interactor3D::SOLID));
			}
         //walls
         D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcObst2,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcObst2,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcObst2,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcObst2,Interactor3D::SOLID));

		  //for shear strees =0
		 // D3Q27BoundaryConditionAdapterPtr velBCAdapter2(new D3Q27VelocityBCAdapter ());
         // velBCAdapter2->setSecondaryBcOption(1);
         // D3Q27InteractorPtr addWallYminInt  = D3Q27InteractorPtr( new D3Q27Interactor(addWallYmin, grid, velBCAdapter2, Interactor3D::SOLID));
		 // D3Q27InteractorPtr addWallZminInt  = D3Q27InteractorPtr( new D3Q27Interactor(addWallZmin, grid, velBCAdapter2, Interactor3D::SOLID));
		 // D3Q27InteractorPtr addWallYmaxInt  = D3Q27InteractorPtr( new D3Q27Interactor(addWallYmax, grid, velBCAdapter2, Interactor3D::SOLID));
		 // D3Q27InteractorPtr addWallZmaxInt  = D3Q27InteractorPtr( new D3Q27Interactor(addWallZmax, grid, velBCAdapter2, Interactor3D::SOLID));

		
		 
         mu::Parser fct;
         //fct.SetExpr("16*U*x2*x3*(H-x2)*(H-x3)/H^4");
         //fct.DefineConst("U", uLB);
         //fct.DefineConst("H", H);

		 fct.SetExpr("U");
         fct.DefineConst("U", uLB);
         
         //inflow
         D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
         velBCAdapter->setSecondaryBcOption(0);
         D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));
 
         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));

	for (int i=0;i<numberOfParticles;i++)
			{      
				sd.addInteractor(spherePInt[i]  );
			}
         //sd.addInteractor(cylinderInt);
         // sd.addInteractor(addWallYminInt);
         // sd.addInteractor(addWallZminInt);
         // sd.addInteractor(addWallYmaxInt);
         // sd.addInteractor(addWallZmaxInt);
         sd.addInteractor(inflowInt);
         sd.addInteractor(outflowInt);
if(myid == 0) UBLOG(logINFO,"delete - start"); 
         sd.deleteSolidBlocks();
if(myid == 0) UBLOG(logINFO,"delete - end"); 

         grid->accept( metisVisitor );

         sd.setTransBlocks();

         ppblocks->update(0);
         ppblocks.reset();

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nod = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);

         double needMemAll  = double(nod*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());

         if(myid == 0)
         {
            UBLOG(logINFO,"Number of blocks = " << nob);
            UBLOG(logINFO,"Number of nodes  = " << nod);
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
         }            

         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3));
         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));
         //option = 0 - ohne param., option = 1 - mit param.
         int option = 0;
       LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, option));
      //   LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB_Geier(blocknx1, blocknx2, blocknx3, option));
	 
//	     LBMKernel3DPtr kernel(new  LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3, option));
         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //walls
        // grid->addAndInitInteractor(addWallYminInt);
        // grid->addAndInitInteractor(addWallZminInt);
        // grid->addAndInitInteractor(addWallYmaxInt);
        // grid->addAndInitInteractor(addWallZmaxInt);

         //obstacle
         //grid->addAndInitInteractor(cylinderInt);
			for (int i=0;i<numberOfParticles;i++)
			{      
				grid->addAndInitInteractor(spherePInt[i]  );
			}

         //inflow
         grid->addAndInitInteractor(inflowInt);

         //outflow
         grid->addAndInitInteractor(outflowInt);

         //domain decomposition
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //initialization of distributions
         D3Q27ETInitDistributionsBlockVisitor initVisitor(rhoLB);
         initVisitor.setVx1(fct);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, comm, true));
 if(myid == 0) UBLOG(logINFO,"/grid/nodes");         
		ppgeo->update(0);
		 if(myid == 0) UBLOG(logINFO,"line"<<__LINE__); 
         ppgeo.reset();

         if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 
      }  
      double outTime = 5000.0;
      UbSchedulerPtr visSch(new UbScheduler(outTime));
      visSch->addSchedule(1000, 1000, 10000);
      visSch->addSchedule(10000, 10000, 50000);
      visSch->addSchedule(1000, 1000, 100000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv, comm);

      double fdx = grid->getDeltaX(grid->getFinestInitializedLevel());
      double point1[3] = {0.45, 0.20, 0.205};
      double point2[3] = {0.55, 0.20, 0.205};

      D3Q27IntegrateValuesHelperPtr h1(new D3Q27IntegrateValuesHelper(grid, comm, 
         point1[0]-1.0*fdx, point1[1]-1.0*fdx, point1[2]-1.0*fdx, 
         point1[0], point1[1], point1[2]));
      if(myid ==0) GbSystem3D::writeGeoObject(h1->getBoundingBox().get(),pathname + "/geo/iv1", WbWriterVtkXmlBinary::getInstance());
      D3Q27IntegrateValuesHelperPtr h2(new D3Q27IntegrateValuesHelper(grid, comm, 
         point2[0], point2[1]-1.0*fdx, point2[2]-1.0*fdx, 
         point2[0]+1.0*fdx, point2[1], point2[2]));
      if(myid ==0) GbSystem3D::writeGeoObject(h2->getBoundingBox().get(),pathname + "/geo/iv2", WbWriterVtkXmlBinary::getInstance());
      //D3Q27PressureDifferencePostprocessor rhopp(grid, visSch, pathname + "/results/rho_diff.txt", h1, h2, conv, comm);
      D3Q27PressureDifferencePostprocessor rhopp(grid, visSch, pathname + "/results/rho_diff.txt", h1, h2, rhoReal, uReal, uLB, comm);
    
      
      double v    = uLB;//4.0*uLB/9.0;
    //  D3Q27ForcesPostprocessor fp(grid, visSch, pathname + "/results/forces.txt", comm, rhoLB, v, area, D3Q27ForcesPostprocessor::X, D3Q27ForcesPostprocessor::Y, D3Q27ForcesPostprocessor::Z);
    //      for (int i=0;i<numberOfParticles;i++)
			 //{      
				// fp.addInteractor(spherePInt[i]  );
			 //}
	  
      UbSchedulerPtr nupsSch(new UbScheduler(10, 10, 40));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      double endTime = 65001.0;
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


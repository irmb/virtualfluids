#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;


void run(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathname = config.getString("pathname");
	  string		  pathGeo = config.getString("pathGeo");
	  string		  geoFile = config.getString("geoFile");
      int             numOfThreads = config.getInt("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
	  vector<double>  boundingBox = config.getVector<double>("boundingBox");
      //vector<double>  length = config.getVector<double>("length");
	  double          uLB = config.getDouble("uLB");
	  double          uF2 = config.getDouble("uF2");
	  double		  nuL = config.getDouble("nuL");
	  double		  nuG = config.getDouble("nuG");
	  double		  densityRatio = config.getDouble("densityRatio");
	  double		  sigma = config.getDouble("sigma");
	  int		      interfaceThickness = config.getInt("interfaceThickness");
	  double		  radius = config.getDouble("radius");
	  double		  theta = config.getDouble("contactAngle");
	  double		  gr = config.getDouble("gravity");
	  double		  phiL = config.getDouble("phi_L");
	  double		  phiH = config.getDouble("phi_H");
	  double		  tauH = config.getDouble("Phase-field Relaxation");
	  double		  mob = config.getDouble("Mobility");


      double          endTime = config.getDouble("endTime");
      double          outTime = config.getDouble("outTime");
      double          availMem = config.getDouble("availMem");
      int             refineLevel = config.getInt("refineLevel");
      double          Re = config.getDouble("Re");
      double          dx = config.getDouble("dx");
      bool            logToFile = config.getBool("logToFile");
      double          restartStep = config.getDouble("restartStep");
      double          cpStart = config.getValue<double>("cpStart");
      double          cpStep = config.getValue<double>("cpStep");
      bool            newStart = config.getValue<bool>("newStart");


	  bool            eastBoundary   = config.getValue<bool>("eastBoundary");
	  bool            westBoundary   = config.getValue<bool>("westBoundary");
	  bool            northBoundary  = config.getValue<bool>("northBoundary");
	  bool            southBoundary  = config.getValue<bool>("southBoundary");
	  bool            topBoundary    = config.getValue<bool>("topBoundary");
	  bool            bottomBoundary = config.getValue<bool>("bottomBoundary");
	  
	  int             eastBoundaryType   = config.getInt("eastBoundaryType");
	  int             westBoundaryType   = config.getInt("westBoundaryType");
	  int             northBoundaryType  = config.getInt("northBoundaryType");
	  int             southBoundaryType  = config.getInt("southBoundaryType");
	  int             topBoundaryType    = config.getInt("topBoundaryType");
	  int             bottomBoundaryType = config.getInt("bottomBoundaryType");


      double beta  = 12*sigma/interfaceThickness;
	  double kappa = 1.5*interfaceThickness*sigma;
	  
	  CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (logToFile)
      {
#if defined(__unix__)
         if (myid == 0)
         {
            const char* str = pathname.c_str();
            mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         }
#endif 

         if (myid == 0)
         {
            stringstream logFilename;
            logFilename << pathname + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }

      //Sleep(30000);

      LBMReal dLB; // = length[1] / dx;
      LBMReal rhoLB = 0.0;
      LBMReal nuLB = nuL; //(uLB*dLB) / Re;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;

 
      
      Grid3DPtr grid(new Grid3D(comm));
      //grid->setPeriodicX1(true);
	  //grid->setPeriodicX2(true);
	  //grid->setPeriodicX3(true);
      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(cpStep, cpStart));
      //RestartCoProcessor rp(grid, rSch, comm, pathname, RestartCoProcessor::TXT);
      MPIIORestart1CoProcessor rcp(grid, rSch, pathname, comm);
      //////////////////////////////////////////////////////////////////////////
	  

	  
      
	  
	  mu::Parser fctF1;
	  //fctF1.SetExpr("vy1*(1-((x1-x0)^2+(x3-z0)^2)/(R^2))");
	  fctF1.SetExpr("vy1");
	  fctF1.DefineConst("vy1", -uLB);
	  fctF1.DefineConst("R", 7.5);
	  fctF1.DefineConst("x0", 60.0);
	  fctF1.DefineConst("z0", 60.0);
	  

	  if (newStart)
      {

         //bounding box
         /*double g_minX1 = 0.0;
         double g_minX2 = -length[1] / 2.0;
         double g_minX3 = -length[2] / 2.0;

         double g_maxX1 = length[0];
         double g_maxX2 = length[1] / 2.0;
         double g_maxX3 = length[2] / 2.0;*/

		 double g_minX1 = boundingBox[0];
		 double g_minX2 = boundingBox[2];
		 double g_minX3 = boundingBox[4];

		 double g_maxX1 = boundingBox[1];
		 double g_maxX2 = boundingBox[3];
		 double g_maxX3 = boundingBox[5];

         //geometry

         //GbObject3DPtr innerCube(new GbCuboid3D(g_minX1+2, g_minX2+2, g_minX3+2, g_maxX1-2, g_maxX2-2, g_maxX3-2));

		 //GbObject3DPtr cylinder1(new GbCylinder3D(g_minX1 - 2.0*dx, g_maxX2/2, g_maxX3/2, g_minX1 + 12.0*dx, g_maxX2/2, g_maxX3/2, radius));
		 //GbObject3DPtr cylinder2(new GbCylinder3D(g_minX1 + 12.0*dx, g_maxX2/2, g_maxX3/2, g_maxX1 + 2.0*dx, g_maxX2/2, g_maxX3/2, dLB / 2.0));
		 
		 //GbObject3DPtr cylinder(new GbCylinder3D(g_minX1 - 2.0*dx, g_maxX2/2, g_maxX3/2, g_maxX1 + 2.0*dx, g_maxX2/2, g_maxX3/2, dLB / 2.0));
		 //GbObject3DPtr cylinders(new GbObject3DManager());
		 //GbObject3DPtr cylinders1(new GbObjectGroup3D());
		 

		 
		 
		 GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

		 //GbTriFaceMesh3DPtr cylinder;
		 //if (myid==0) UBLOG(logINFO, "Read geoFile:start");
		 //cylinder = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+geoFile, "geoCylinders", GbTriFaceMesh3D::KDTREE_SAHPLIT));
		 //GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/Stlgeo", WbWriterVtkXmlBinary::getInstance());
		 
		 /*GbObject3DPtr cylinder(new GbCuboid3D(g_minX1 + 2.0, g_minX2 + 2.0, g_minX3 + 2.0, g_maxX1 - 2.0, g_maxX2 - 2.0, g_maxX3 -2.0));
		 if (myid == 0) GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/solidWall", WbWriterVtkXmlBinary::getInstance());*/

		 
		 
		 bool isDefined_SlipBC   = false;
		 bool isDefined_noSlipBC = false;
		 bool isDefined_Outflow  = false;

		 if (eastBoundary)
		 {
			 switch (eastBoundaryType)
			 {
			 case 0:

			 	break;
			 case 1:
				 break;
			 case 2:
				 break;
			 case 3:
				 break;
			 }
		 }

		 /*GbObject3DPtr geoWallEast(new GbCuboid3D(g_maxX1 + 2.0, g_maxX2 + 2.0, g_maxX3 + 2.0, g_maxX1 - 2.0, g_minX2 - 2.0, g_minX3 -2.0));
		 if (myid == 0) GbSystem3D::writeGeoObject(geoWallEast.get(), pathname + "/geo/geoWallEast", WbWriterVtkXmlBinary::getInstance());

		 GbObject3DPtr geoWallWest(new GbCuboid3D(g_minX1 + 2.0, g_maxX2 + 2.0, g_maxX3 + 2.0, g_minX1 - 2.0, g_minX2 - 2.0, g_minX3 -2.0));
		 if (myid == 0) GbSystem3D::writeGeoObject(geoWallWest.get(), pathname + "/geo/geoWallWest", WbWriterVtkXmlBinary::getInstance());

		 GbObject3DPtr geoWallTop(new GbCuboid3D(g_maxX1 + 2.0, g_maxX2 + 2.0, g_maxX3 + 2.0, g_minX1 - 2.0, g_maxX2 - 2.0, g_minX3 -2.0));
		 if (myid == 0) GbSystem3D::writeGeoObject(geoWallTop.get(), pathname + "/geo/geoWallTop", WbWriterVtkXmlBinary::getInstance());

		 GbObject3DPtr geoWallBottom(new GbCuboid3D(g_maxX1 + 2.0, g_minX2 + 2.0, g_maxX3 + 2.0, g_minX1 - 2.0, g_minX2 - 2.0, g_minX3 -2.0));
		 if (myid == 0) GbSystem3D::writeGeoObject(geoWallBottom.get(), pathname + "/geo/geoWallBottom", WbWriterVtkXmlBinary::getInstance());

		 GbObject3DPtr geoWallNorth(new GbCuboid3D(g_maxX1 + 2.0, g_maxX2 + 2.0, g_minX3 + 2.0, g_minX1 - 2.0, g_minX2 - 2.0, g_minX3 -2.0));
		 if (myid == 0) GbSystem3D::writeGeoObject(geoWallNorth.get(), pathname + "/geo/geoWallNorth", WbWriterVtkXmlBinary::getInstance());

		 GbObject3DPtr geoWallSouth(new GbCuboid3D(g_maxX1 + 2.0, g_maxX2 + 2.0, g_maxX3 + 2.0, g_minX1 - 2.0, g_minX2 - 2.0, g_maxX3 -2.0));
		 if (myid == 0) GbSystem3D::writeGeoObject(geoWallSouth.get(), pathname + "/geo/geoWallSouth", WbWriterVtkXmlBinary::getInstance());*/
		 
		 
		 //inflow
		 //GbCuboid3DPtr geoInflowF1(new GbCuboid3D(40.0, 628.0, 40.0, 80, 631.0, 80.0));  // Original
		 GbCuboid3DPtr geoInflowF1(new GbCuboid3D(g_minX1, g_minX2-0.5*dx, g_minX3, g_maxX1, g_minX2 - 1.0*dx, g_maxX3));
		 if (myid==0) GbSystem3D::writeGeoObject(geoInflowF1.get(), pathname+"/geo/geoInflowF1", WbWriterVtkXmlASCII::getInstance());



		 //outflow
		 //GbCuboid3DPtr geoOutflow(new GbCuboid3D(-1.0, 499+80, -1.0, 121.0, 501.0+80, 121.0)); // Original
		 GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_minX1, g_maxX2-1*dx, g_minX3, g_maxX1, g_maxX2, g_maxX3));
		 if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());
		 
		 
		 double blockLength = blocknx[0] * dx;




         if (myid == 0)
         {
            UBLOG(logINFO, "uLb = " << uLB);
            UBLOG(logINFO, "rho = " << rhoLB);
            UBLOG(logINFO, "nuLb = " << nuLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << dx);
            UBLOG(logINFO, "Preprocess - start");
         }

         grid->setDeltaX(dx);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

		 grid->setPeriodicX1(true);
		 grid->setPeriodicX2(true);
		 grid->setPeriodicX3(true);

         

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);




		 //BC Adapter
		 //////////////////////////////////////////////////////////////////////////////
		 BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
		 noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithmMultiphase()));


		BCAdapterPtr denBCAdapter(new DensityBCAdapter(rhoLB));
		denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithmMultiphase()));
		
		double r = 5.0; //boost::dynamic_pointer_cast<GbCylinder3D>(cylinder)->getRadius();
		double cx1 = g_minX1;
		double cx2 = 0.0; //cylinder->getX2Centroid();
		double cx3 = 0.0; //cylinder->getX3Centroid();


		
		mu::Parser fctPhi_F1;
		fctPhi_F1.SetExpr("phiH");
		fctPhi_F1.DefineConst("phiH", phiH);

		mu::Parser fctPhi_F2;
		fctPhi_F2.SetExpr("phiL");
		fctPhi_F2.DefineConst("phiL", phiL);
		
		mu::Parser fctvel_F2_init;
		fctvel_F2_init.SetExpr("U");
		fctvel_F2_init.DefineConst("U", 0);

		//fct.SetExpr("U");
		//fct.DefineConst("U", uLB);
		//BCAdapterPtr velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
		
		BCAdapterPtr velBCAdapterF1       (new VelocityBCAdapterMultiphase(false, true, false, fctF1  , phiH, 0.0, endTime));
		
		//BCAdapterPtr velBCAdapterF2_1_init(new VelocityBCAdapterMultiphase(false, false, true, fctF2_1, phiH, 0.0, endTime));
		//BCAdapterPtr velBCAdapterF2_2_init(new VelocityBCAdapterMultiphase(false, false, true, fctF2_2, phiH, 0.0, endTime));

		//BCAdapterPtr velBCAdapterF2_1_init(new VelocityBCAdapterMultiphase(false, false, true, fctvel_F2_init, phiL, 0.0, endTime));
		//BCAdapterPtr velBCAdapterF2_2_init(new VelocityBCAdapterMultiphase(false, false, true, fctvel_F2_init, phiL, 0.0, endTime));
		
		velBCAdapterF1->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithmMultiphase()));
		//velBCAdapterF2_1_init->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithmMultiphase()));
		//velBCAdapterF2_2_init->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithmMultiphase()));
		 
		 
		 //velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityWithDensityBCAlgorithm()));
		 //mu::Parser fct;
		 //fct.SetExpr("U");
		 //fct.DefineConst("U", uLB);
		 //BCAdapterPtr velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
		 //velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingVelocityBCAlgorithm()));


		 //////////////////////////////////////////////////////////////////////////////////
		 //BC visitor
		 BoundaryConditionsBlockVisitorMultiphase bcVisitor;
		 bcVisitor.addBC(noSlipBCAdapter);
		 bcVisitor.addBC(denBCAdapter);
		 bcVisitor.addBC(velBCAdapterF1);
		 //bcVisitor.addBC(velBCAdapterF2_1_init);
		 //bcVisitor.addBC(velBCAdapterF2_2_init);



         WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         ppblocks->process(0);

         //Interactor3DPtr tubes(new D3Q27TriFaceMeshInteractor(cylinder, grid, noSlipBCAdapter, Interactor3D::SOLID));
		 //D3Q27InteractorPtr tubes = D3Q27InteractorPtr(new D3Q27Interactor(cylinder, grid, noSlipBCAdapter, Interactor3D::INVERSESOLID));

	     
		 
		 /*D3Q27InteractorPtr wallEast   = D3Q27InteractorPtr(new D3Q27Interactor(geoWallEast, grid, noSlipBCAdapter, Interactor3D::SOLID));
		 D3Q27InteractorPtr wallWest   = D3Q27InteractorPtr(new D3Q27Interactor(geoWallWest, grid, noSlipBCAdapter, Interactor3D::SOLID));
		 D3Q27InteractorPtr wallTop    = D3Q27InteractorPtr(new D3Q27Interactor(geoWallTop, grid, noSlipBCAdapter, Interactor3D::SOLID));
		 D3Q27InteractorPtr wallBottom = D3Q27InteractorPtr(new D3Q27Interactor(geoWallBottom, grid, noSlipBCAdapter, Interactor3D::SOLID));
		 D3Q27InteractorPtr wallNorth  = D3Q27InteractorPtr(new D3Q27Interactor(geoWallNorth, grid, noSlipBCAdapter, Interactor3D::SOLID));
		 D3Q27InteractorPtr wallSouth  = D3Q27InteractorPtr(new D3Q27Interactor(geoWallSouth, grid, noSlipBCAdapter, Interactor3D::SOLID));*/
		 
		 D3Q27InteractorPtr inflowF1Int = D3Q27InteractorPtr(new D3Q27Interactor(geoInflowF1, grid, velBCAdapterF1, Interactor3D::SOLID));

         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));

		 //SetSolidBlockVisitor visitor1(inflowF2_1Int, SetSolidBlockVisitor::BC);
		 //grid->accept(visitor1);
		 //SetSolidBlockVisitor visitor2(inflowF2_2Int, SetSolidBlockVisitor::BC);
		 //grid->accept(visitor2);


         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW));
         InteractorsHelper intHelper(grid, metisVisitor);
		 //intHelper.addInteractor(tubes);
		 //intHelper.addInteractor(inflowF1Int);
		 //intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();


         ppblocks->process(0);
         ppblocks.reset();

         unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
         int ghostLayer = 3;
         unsigned long long numberOfNodesPerBlock = (unsigned long long)(blocknx[0])* (unsigned long long)(blocknx[1])* (unsigned long long)(blocknx[2]);
         unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
         unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
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

         LBMKernelPtr kernel;

         kernel = LBMKernelPtr(new MultiphaseCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2], MultiphaseCumulantLBMKernel::NORMAL));

         kernel->setWithForcing(true);
		 kernel->setForcingX1(0.0);
		 kernel->setForcingX2(gr);
		 kernel->setForcingX3(0.0);

		 kernel->setPhiL(phiL);
		 kernel->setPhiH(phiH);
		 kernel->setPhaseFieldRelaxation(tauH);
		 kernel->setMobility(mob);

         BCProcessorPtr bcProc(new BCProcessor());
         //BCProcessorPtr bcProc(new ThinWallBCProcessor());

         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitorMultiphase kernelVisitor(kernel, nuL, nuG, densityRatio, beta, kappa, theta, availMem, needMem);
         
		 grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

		 //inflowF2_1Int->initInteractor();
		 //inflowF2_2Int->initInteractor();

         intHelper.setBC();
		 
        
         grid->accept(bcVisitor);

         //initialization of distributions
		 LBMReal x1c =  (g_maxX1+g_minX1)/2; //radius; //g_minX1; //radius; //19; //(g_maxX1+g_minX1)/2;
		 LBMReal x2c = (g_maxX2+g_minX2)/2; //g_minX2 + 2;
		 LBMReal x3c = (g_maxX3+g_minX3)/2;
		 mu::Parser fct1;
		 
		 //fct1.SetExpr("0.5-0.5*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
		 //fct1.SetExpr("phiM-phiM*tanh((sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/(interfaceThickness*phiM))");
		 
		 fct1.SetExpr("0.5*(phiH + phiL) - 0.5*(phiH - phiL)*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
		 
		 
		 //fct1.SetExpr("0.5*(phiH + phiL) + 0.5*(phiH - phiL)*tanh(2*((x2-radius))/interfaceThickness)");
		 //fct1.SetExpr("phiL");
		 fct1.DefineConst("x1c", x1c);
		 fct1.DefineConst("x2c", x2c);
		 fct1.DefineConst("x3c", x3c);
		 fct1.DefineConst("phiL", phiL);
		 fct1.DefineConst("phiH", phiH);
		 fct1.DefineConst("radius", radius);
		 fct1.DefineConst("interfaceThickness", interfaceThickness);
		 
		 mu::Parser fct2;
		 //fct2.SetExpr("vx1*(1-((x2-y0)^2+(x3-z0)^2)/(R^2))");
		 /*fct2.SetExpr("vx1");
		 fct2.DefineConst("R", 10.0);
		 fct2.DefineConst("vx1", uLB);
		 fct2.DefineConst("y0", 1.0);
		 fct2.DefineConst("z0", 31.0);*/
		 fct2.SetExpr("0.5*uLB*(phiH + phiL) - 0.5*uLB*(phiH - phiL)*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
		 fct2.DefineConst("uLB", uLB);
		 fct2.DefineConst("x1c", x1c);
		 fct2.DefineConst("x2c", x2c);
		 fct2.DefineConst("x3c", x3c);
		 fct2.DefineConst("phiL", phiL);
		 fct2.DefineConst("phiH", phiH);
		 fct2.DefineConst("radius", radius);
		 fct2.DefineConst("interfaceThickness", interfaceThickness);


		 mu::Parser fct3;
		 fct3.SetExpr("0.5*sigma*(phiH + phiL)/radius - 0.5*sigma/radius*(phiH - phiL)*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
		 fct3.DefineConst("sigma", sigma);
		 fct3.DefineConst("x1c", x1c);
		 fct3.DefineConst("x2c", x2c);
		 fct3.DefineConst("x3c", x3c);
		 fct3.DefineConst("phiL", phiL);
		 fct3.DefineConst("phiH", phiH);
		 fct3.DefineConst("radius", radius);
		 fct3.DefineConst("interfaceThickness", interfaceThickness);

		 InitDistributionsBlockVisitorMultiphase initVisitor(densityRatio, interfaceThickness, radius);
         initVisitor.setPhi(fct1);
         //initVisitor.setVx1(fct2);
		 //initVisitor.setRho(fct3);
         grid->accept(initVisitor);

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         //InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         //ConnectorFactoryPtr factory(new Block3DConnectorFactory());
         //ConnectorBlockVisitor setConnsVisitor(comm, nuLB, iProcessor, factory);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);




         //boundary conditions grid
         {
            UbSchedulerPtr geoSch(new UbScheduler(1));
            WriteBoundaryConditionsCoProcessorPtr ppgeo(
               new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
            ppgeo->process(0);
            ppgeo.reset();
         }

         if (myid == 0) UBLOG(logINFO, "Preprocess - end");
      }
      else
      {
         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "uLb = " << uLB);
            UBLOG(logINFO, "rho = " << rhoLB);
            UBLOG(logINFO, "nuLb = " << nuLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << dx);
            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
         }

         rcp.restart((int)restartStep);
         grid->setTimeStep(restartStep);

         //BCAdapterPtr velBCAdapter(new VelocityBCAdapter());
         //velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithm()));
         //velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityWithDensityBCAlgorithm()));
         //bcVisitor.addBC(velBCAdapter);
         //grid->accept(bcVisitor);

         //set connectors
         //InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      UbSchedulerPtr visSch(new UbScheduler(outTime));
      WriteMacroscopicQuantitiesCoProcessor pp(grid, visSch, pathname, WbWriterVtkXmlASCII::getInstance(), conv, comm);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

	  
	  
	  

	  
	  //UbSchedulerPtr bcSch(new UbScheduler(1, 12000, 12000));
	  //TimeDependentBCCoProcessorPtr inflowF2 (new TimeDependentBCCoProcessor(grid,bcSch));
	  //inflowF2->addInteractor(inflowF2_1Int);
	  //inflowF2->addInteractor(inflowF2_2Int);

      //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch,CalculationManager::MPI));
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
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
   //Sleep(30000);
	if (argv != NULL)
   {
      if (argv[1] != NULL)
      {
         run(string(argv[1]));
      }
      else
      {
         cout << "Configuration file is missing!" << endl;
      }
   }

}


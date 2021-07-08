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


	  bool            isEastBoundary   = config.getValue<bool>("isEastBoundary");
	  bool            isWestBoundary   = config.getValue<bool>("isWestBoundary");
	  bool            isNorthBoundary  = config.getValue<bool>("isNorthBoundary");
	  bool            isSouthBoundary  = config.getValue<bool>("isSouthBoundary");
	  bool            isTopBoundary    = config.getValue<bool>("isTopBoundary");
	  bool            isBottomBoundary = config.getValue<bool>("isBottomBoundary");
	  bool            isPeriodicX1     = config.getValue<bool>("isPeriodicX1");
	  bool            isPeriodicX2     = config.getValue<bool>("isPeriodicX2");
	  bool            isPeriodicX3     = config.getValue<bool>("isPeriodicX3");
	  
	  int             eastBoundaryType   = config.getInt("eastBoundaryType");
	  int             westBoundaryType   = config.getInt("westBoundaryType");
	  int             northBoundaryType  = config.getInt("northBoundaryType");
	  int             southBoundaryType  = config.getInt("southBoundaryType");
	  int             topBoundaryType    = config.getInt("topBoundaryType");
	  int             bottomBoundaryType = config.getInt("bottomBoundaryType");

	  bool            isInitialVelocity     = config.getValue<bool>("isInitialVelocity");
	  string          phaseFieldProfile = config.getString("phaseFieldProfile");
	  string          velocityProfile = config.getString("velocityProfile");

	  double          p_in = config.getDouble("p_in");
	  double          p_out = config.getDouble("p_out");



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

		 grid->setPeriodicX1(isPeriodicX1);
		 grid->setPeriodicX2(isPeriodicX2);
		 grid->setPeriodicX3(isPeriodicX3);



		 GenBlocksGridVisitor genBlocks(gridCube);
		 grid->accept(genBlocks);

		 WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

		 ppblocks->process(0);

		 BoundaryConditionsBlockVisitorMultiphase bcVisitor;
		 
		 Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW));
		 InteractorsHelper intHelper(grid, metisVisitor);


		 
		 		 
		 
		 if (isEastBoundary)
		 {
			 GbCuboid3DPtr geoEastBoundary(new GbCuboid3D(g_maxX1 + 1.0*dx, g_minX2 - 1.0*dx, g_minX3 - 1.0*dx, g_maxX1 - 1.0*dx, g_maxX2 + 1.0*dx, g_maxX3 + 1.0*dx));
			 if (myid==0) GbSystem3D::writeGeoObject(geoEastBoundary.get(), pathname+"/geo/geoEastBoundary", WbWriterVtkXmlASCII::getInstance());
			 
			 if (eastBoundaryType == 1)
			 {
				 BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
				 noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(noSlipBCAdapter);
				 D3Q27InteractorPtr eastBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoEastBoundary, grid, noSlipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(eastBoundaryInt);

			 }

			 if (eastBoundaryType == 4)
			 {
				 BCAdapterPtr slipBCAdapter(new SlipBCAdapter());
				 slipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new SlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(slipBCAdapter);
				 D3Q27InteractorPtr eastBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoEastBoundary, grid, slipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(eastBoundaryInt);
			 }

			 if (eastBoundaryType == 2)
			 {
				 BCAdapterPtr velBCAdapter (new VelocityBCAdapterMultiphase(true, false, false, fctF1  , phiH, 0.0, endTime));
				 velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithmMultiphase()));
				 bcVisitor.addBC(velBCAdapter);
				 D3Q27InteractorPtr eastBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoEastBoundary, grid, velBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(eastBoundaryInt);
			 }

			 if (eastBoundaryType == 3)
			 {
				 BCAdapterPtr denBCAdapter(new DensityBCAdapter(rhoLB));
				 denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithmMultiphase()));
				 bcVisitor.addBC(denBCAdapter);
				 D3Q27InteractorPtr eastBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoEastBoundary, grid, denBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(eastBoundaryInt);
			 }
			 if (eastBoundaryType == 5)
			 {
				 BCAdapterPtr denBCAdapter(new DensityBCAdapter(p_out));
				 denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new EqDensityBCAlgorithm()));
				 bcVisitor.addBC(denBCAdapter);
				 D3Q27InteractorPtr eastBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoEastBoundary, grid, denBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(eastBoundaryInt);
			 }
			
		 }

		 
		 if (isWestBoundary)
		 {
			 GbCuboid3DPtr geoWestBoundary(new GbCuboid3D(g_minX1 + 1.0*dx, g_minX2 - 1.0*dx, g_minX3 - 1.0*dx, g_minX1 - 1.0*dx, g_maxX2 + 1.0*dx, g_maxX3 + 1.0*dx));
			 if (myid==0) GbSystem3D::writeGeoObject(geoWestBoundary.get(), pathname+"/geo/geoWestBoundary", WbWriterVtkXmlASCII::getInstance());

			 if (westBoundaryType == 1)
			 {
				 BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
				 noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(noSlipBCAdapter);
				 D3Q27InteractorPtr westBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoWestBoundary, grid, noSlipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(westBoundaryInt);

			 }

			 if (westBoundaryType == 4)
			 {
				 BCAdapterPtr slipBCAdapter(new SlipBCAdapter());
				 slipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new SlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(slipBCAdapter);
				 D3Q27InteractorPtr westBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoWestBoundary, grid, slipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(westBoundaryInt);
			 }

			 if (westBoundaryType == 2)
			 {
				 BCAdapterPtr velBCAdapter (new VelocityBCAdapterMultiphase(true, false, false, fctF1  , phiH, 0.0, endTime));
				 velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithmMultiphase()));
				 bcVisitor.addBC(velBCAdapter);
				 D3Q27InteractorPtr westBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoWestBoundary, grid, velBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(westBoundaryInt);
			 }

			 if (westBoundaryType == 3)
			 {
				 BCAdapterPtr denBCAdapter(new DensityBCAdapter(rhoLB));
				 denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithmMultiphase()));
				 bcVisitor.addBC(denBCAdapter);
				 D3Q27InteractorPtr westBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoWestBoundary, grid, denBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(westBoundaryInt);
			 }
			 if (westBoundaryType == 5)
			 {
				 BCAdapterPtr denBCAdapter(new DensityBCAdapter(p_in));
				 denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new EqDensityBCAlgorithm()));
				 bcVisitor.addBC(denBCAdapter);
				 D3Q27InteractorPtr westBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoWestBoundary, grid, denBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(westBoundaryInt);
			 }

		 }		 

		 
		 if (isNorthBoundary)
		 {
			 GbCuboid3DPtr geoNorthBoundary(new GbCuboid3D(g_minX1 - 1.0*dx, g_maxX2 + 1.0*dx, g_minX3 + 1.0*dx, g_maxX1 + 1.0*dx, g_minX2 - 1.0*dx, g_minX3 - 1.0*dx));
			 if (myid==0) GbSystem3D::writeGeoObject(geoNorthBoundary.get(), pathname+"/geo/geoNorthBoundary", WbWriterVtkXmlASCII::getInstance());

			 if (northBoundaryType == 1)
			 {
				 BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
				 noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(noSlipBCAdapter);
				 D3Q27InteractorPtr northBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoNorthBoundary, grid, noSlipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(northBoundaryInt);

			 }

			 if (northBoundaryType == 4)
			 {
				 BCAdapterPtr slipBCAdapter(new SlipBCAdapter());
				 slipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new SlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(slipBCAdapter);
				 D3Q27InteractorPtr northBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoNorthBoundary, grid, slipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(northBoundaryInt);
			 }

			 if (northBoundaryType == 2)
			 {
				 BCAdapterPtr velBCAdapter (new VelocityBCAdapterMultiphase(true, false, false, fctF1  , phiH, 0.0, endTime));
				 velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithmMultiphase()));
				 bcVisitor.addBC(velBCAdapter);
				 D3Q27InteractorPtr northBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoNorthBoundary, grid, velBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(northBoundaryInt);
			 }

			 if (northBoundaryType == 3)
			 {
				 BCAdapterPtr denBCAdapter(new DensityBCAdapter(rhoLB));
				 denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithmMultiphase()));
				 bcVisitor.addBC(denBCAdapter);
				 D3Q27InteractorPtr northBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoNorthBoundary, grid, denBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(northBoundaryInt);
			 }

		 }			 
		 
		 
		 if (isSouthBoundary)
		 {
			 GbCuboid3DPtr geoSouthBoundary(new GbCuboid3D(g_minX1 - 1.0*dx, g_maxX2 + 1.0*dx, g_maxX3 + 1.0*dx, g_maxX1 + 1.0*dx, g_minX2 - 1.0*dx, g_maxX3 - 1.0*dx));
			 if (myid==0) GbSystem3D::writeGeoObject(geoSouthBoundary.get(), pathname+"/geo/geoSouthBoundary", WbWriterVtkXmlASCII::getInstance());

			 if (southBoundaryType == 1)
			 {
				 BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
				 noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(noSlipBCAdapter);
				 D3Q27InteractorPtr southBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoSouthBoundary, grid, noSlipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(southBoundaryInt);

			 }

			 if (southBoundaryType == 4)
			 {
				 BCAdapterPtr slipBCAdapter(new SlipBCAdapter());
				 slipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new SlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(slipBCAdapter);
				 D3Q27InteractorPtr southBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoSouthBoundary, grid, slipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(southBoundaryInt);
			 }

			 if (southBoundaryType == 2)
			 {
				 BCAdapterPtr velBCAdapter (new VelocityBCAdapterMultiphase(true, false, false, fctF1  , phiH, 0.0, endTime));
				 velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithmMultiphase()));
				 bcVisitor.addBC(velBCAdapter);
				 D3Q27InteractorPtr southBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoSouthBoundary, grid, velBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(southBoundaryInt);
			 }

			 if (southBoundaryType == 3)
			 {
				 BCAdapterPtr denBCAdapter(new DensityBCAdapter(rhoLB));
				 denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithmMultiphase()));
				 bcVisitor.addBC(denBCAdapter);
				 D3Q27InteractorPtr southBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoSouthBoundary, grid, denBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(southBoundaryInt);
			 }

		 }




		 if (isTopBoundary)
		 {
			 GbCuboid3DPtr geoTopBoundary(new GbCuboid3D(g_minX1 - 1.0*dx, g_maxX2 + 1.0*dx, g_minX3 - 1.0*dx, g_maxX1 + 1.0*dx, g_maxX2 - 1.0*dx, g_maxX3 + 1.0*dx));
			 if (myid==0) GbSystem3D::writeGeoObject(geoTopBoundary.get(), pathname+"/geo/geoTopBoundary", WbWriterVtkXmlASCII::getInstance());

			 if (topBoundaryType == 1)
			 {
				 BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
				 noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(noSlipBCAdapter);
				 D3Q27InteractorPtr topBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoTopBoundary, grid, noSlipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(topBoundaryInt);

			 }

			 if (topBoundaryType == 4)
			 {
				 BCAdapterPtr slipBCAdapter(new SlipBCAdapter());
				 slipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new SlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(slipBCAdapter);
				 D3Q27InteractorPtr topBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoTopBoundary, grid, slipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(topBoundaryInt);
			 }

			 if (topBoundaryType == 2)
			 {
				 BCAdapterPtr velBCAdapter (new VelocityBCAdapterMultiphase(true, false, false, fctF1  , phiH, 0.0, endTime));
				 velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithmMultiphase()));
				 bcVisitor.addBC(velBCAdapter);
				 D3Q27InteractorPtr topBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoTopBoundary, grid, velBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(topBoundaryInt);
			 }

			 if (topBoundaryType == 3)
			 {
				 BCAdapterPtr denBCAdapter(new DensityBCAdapter(rhoLB));
				 denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithmMultiphase()));
				 bcVisitor.addBC(denBCAdapter);
				 D3Q27InteractorPtr topBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoTopBoundary, grid, denBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(topBoundaryInt);
			 }

		 }			 


		 if (isBottomBoundary)
		 {
			 GbCuboid3DPtr geoBottomBoundary(new GbCuboid3D(g_minX1 - 1.0*dx, g_minX2 + 1.0*dx, g_minX3 - 1.0*dx, g_maxX1 + 1.0*dx, g_minX2 - 1.0*dx, g_maxX3 + 1.0*dx));
			 if (myid==0) GbSystem3D::writeGeoObject(geoBottomBoundary.get(), pathname+"/geo/geoBottomBoundary", WbWriterVtkXmlASCII::getInstance());

			 if (bottomBoundaryType == 1)
			 {
				 BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
				 noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(noSlipBCAdapter);
				 D3Q27InteractorPtr bottomBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoBottomBoundary, grid, noSlipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(bottomBoundaryInt);

			 }

			 if (bottomBoundaryType == 4)
			 {
				 BCAdapterPtr slipBCAdapter(new SlipBCAdapter());
				 slipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new SlipBCAlgorithmMultiphase()));
				 bcVisitor.addBC(slipBCAdapter);
				 D3Q27InteractorPtr bottomBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoBottomBoundary, grid, slipBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(bottomBoundaryInt);
			 }

			 if (bottomBoundaryType == 2)
			 {
				 BCAdapterPtr velBCAdapter (new VelocityBCAdapterMultiphase(true, false, false, fctF1  , phiH, 0.0, endTime));
				 velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithmMultiphase()));
				 bcVisitor.addBC(velBCAdapter);
				 D3Q27InteractorPtr bottomBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoBottomBoundary, grid, velBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(bottomBoundaryInt);
			 }

			 if (bottomBoundaryType == 3)
			 {
				 BCAdapterPtr denBCAdapter(new DensityBCAdapter(rhoLB));
				 denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithmMultiphase()));
				 bcVisitor.addBC(denBCAdapter);
				 D3Q27InteractorPtr bottomBoundaryInt = D3Q27InteractorPtr(new D3Q27Interactor(geoBottomBoundary, grid, denBCAdapter, Interactor3D::SOLID));
				 intHelper.addInteractor(bottomBoundaryInt);
			 }

		 }





         



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
		 
		 //fct1.SetExpr("0.5*(phiH + phiL) - 0.5*(phiH - phiL)*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
		 fct1.SetExpr(phaseFieldProfile);
		 
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
		 //fct2.SetExpr("0.5*uLB*(phiH + phiL) - 0.5*uLB*(phiH - phiL)*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
		 fct2.SetExpr(velocityProfile);
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
		 if (isInitialVelocity) initVisitor.setVx1(fct2);
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


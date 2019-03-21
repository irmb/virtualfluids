#include <iostream>
#include <string>
#include "VirtualFluids.h"


using namespace std;
void run(string configname)
{
   try
   {
      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (myid == 0) UBLOG(logINFO, "Testcase organ pipe");

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      ConfigurationFile   config;
      config.load(configname);

      bool            newStart = config.getValue<bool>("newStart");
      double          restartStep = config.getValue<double>("restartStep");
      double          cpStart = config.getValue<double>("cpStart");
      double          cpStep = config.getValue<double>("cpStep");
      int             endTime = config.getValue<int>("endTime");
      double          outTimeStep = config.getValue<double>("outTimeStep");
      double          outTimeStart = config.getValue<double>("outTimeStart");
      double          availMem = config.getValue<double>("availMem");
      bool            logToFile = config.getValue<bool>("logToFile");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");
      string          pathOut = config.getValue<string>("pathOut");
      string          pathGeo = config.getValue<string>("pathGeo");

      string          opipeGeoFile = "/OrganPipeTransformed.stl";
      string          inletTubeGeoFile = "/tubeTransformed.stl";
      //string          extendedGeoFile = "/FluteExtension.stl";

      double  deltaXfine = 0.0000625;
      const int baseLevel = 0;
      int refineLevel = 9;
      double deltaXcoarse = deltaXfine*(double)(1 << refineLevel);

      LBMReal rho_LB = 0.0;
      double rhoReal = 1.2041; //(kg/m3)
      double uReal = 33.9; //m/s
      double L = 0.195; //m
      double hLB = L / deltaXcoarse;
      double csReal = 343.3;
      double nuReal = 1.51e-5; //m^2/s
      LBMUnitConverter unitConverter(L, csReal, rhoReal, hLB);
      if (myid == 0) UBLOG(logINFO, unitConverter.toString());

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

      double nu_LB = nuReal * unitConverter.getFactorViscosityWToLb();
      double u_LB = uReal * unitConverter.getFactorVelocityWToLb();


      vector<int> blocknx ={ 16, 16, 16 };

      if (myid == 0) UBLOG(logINFO, "Read organ pipe geometry:start");
      SPtr<GbTriFaceMesh3D> organPipeGeo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo + opipeGeoFile, "opipeGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
      organPipeGeo->translate(1.37, 0.0, 0.0);
      if (myid == 0) UBLOG(logINFO, "Read organ pipe geometry:end");
      if (myid == 0) GbSystem3D::writeGeoObject(organPipeGeo.get(), pathOut + "/geo/organPipeGeo", WbWriterVtkXmlBinary::getInstance());

      SPtr<Grid3D> grid(new Grid3D(comm));

      //bounding box
      vector<double> dim = {4.096, 4.096, 4.096};

      double g_minX1 = 0;
      double g_minX2 =  -dim[1]*0.5;
      double g_minX3 =  -dim[2]*0.5;

      double g_maxX1 =  dim[0];
      double g_maxX2 =  dim[1]*0.5;
      double g_maxX3 =  dim[2]*0.5;


      double d_minX2 = -dim[1]*0.5;
      double d_minX3 = -dim[2]*0.5;

      double d_maxX1 = dim[0];
      double d_maxX2 = dim[1]*0.5;
      double d_maxX3 = dim[2]*0.5;

      if (myid == 0)
      {
         UBLOG(logINFO, "Parameters:");
         UBLOG(logINFO, "u_Real              = " << uReal   << " [m/s]");
         UBLOG(logINFO, "rho_Real            = " << rhoReal << " [kg/m^3]");
         UBLOG(logINFO, "nu_Real             = " << nuReal  << " [m^2/s]");
         UBLOG(logINFO, "u_LB                = " << u_LB    << " [dx/dt]");
         UBLOG(logINFO, "rho_LB              = " << rho_LB  << " [mass/dx^3]");
         UBLOG(logINFO, "nu_LB               = " << nu_LB   << " [dx^2/dt]");
         UBLOG(logINFO, "dx coarse           = " << deltaXcoarse);
         UBLOG(logINFO, "dx fine             = " << deltaXfine);
         UBLOG(logINFO, "number of refinement levels = " << refineLevel);
         UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
         UBLOG(logINFO, "path = " << pathOut);
         UBLOG(logINFO, "Preprocess - start");
      }

      //BC adapters
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      SPtr<BCAdapter> slipBCAdapter(new SlipBCAdapter());
      slipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SlipBCAlgorithm()));

      double diameter_inlet = 0.005;
      double cx1 = g_minX1;
      double cx2 = 0.0;
      double cx3 = 0.0;

      mu::Parser fct;
      //fct.SetExpr("U");
      fct.SetExpr("U*(1-(((((x2-y0)^2+(x3-z0)^2)^0.5)/R)^19.538))");
      fct.DefineConst("x0", cx1);
      fct.DefineConst("y0", cx2);
      fct.DefineConst("z0", cx3);
      fct.DefineConst("R", diameter_inlet/2.0);
      fct.DefineConst("U", u_LB*1.16182766);

      SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

      SPtr<BCAdapter> outflowBCAdapter(new DensityBCAdapter(rho_LB));
      outflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(slipBCAdapter);
      bcVisitor.addBC(velBCAdapter);
      bcVisitor.addBC(outflowBCAdapter);

      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new BCProcessor());

      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      double bulckViscosity = 10.0*nu_LB;
      dynamicPointerCast<CompressibleCumulant4thOrderViscosityLBMKernel>(kernel)->setBulkViscosity(bulckViscosity);
      kernel->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIOMigrationCoProcessor> migCoProcessor(new MPIIOMigrationCoProcessor(grid, mSch, pathOut + "/mig", comm));
      migCoProcessor->setLBMKernel(kernel);
      migCoProcessor->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////

      if (newStart)
      {
         grid->setPeriodicX1(false);
         grid->setPeriodicX2(false);
         grid->setPeriodicX3(false);
         grid->setDeltaX(deltaXcoarse);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

         //generate block grid
         SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
         GenBlocksGridVisitor genBlocks(gridCube);

         grid->accept(genBlocks);

         //geometry
         if (myid == 0) UBLOG(logINFO, "Read inlet pipe geometry:start");
         SPtr<GbTriFaceMesh3D> inletTubeGeo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo + inletTubeGeoFile, "inPipeGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         inletTubeGeo->translate(1.37, 0.0, 0.0);
         if (myid == 0) UBLOG(logINFO, "Read inlet pipe geometry:end");
         if (myid == 0) GbSystem3D::writeGeoObject(inletTubeGeo.get(), pathOut + "/geo/inletTubeGeo", WbWriterVtkXmlBinary::getInstance());
         SPtr<Interactor3D> organPipeInter = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(organPipeGeo, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
         SPtr<Interactor3D> inletTubeInter = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(inletTubeGeo, grid, noSlipBCAdapter, Interactor3D::SOLID));

         double op_offset = organPipeGeo->getX1Minimum() - inletTubeGeo->getX1Minimum();
         double startX1it = inletTubeGeo->getX1Minimum();
         //////////////////////////////////////////////////////////////////////////
         //refinement

         SPtr<GbObject3D> refineBoxL6(new GbCuboid3D(startX1it, -0.0165, -0.0165, startX1it+0.2634+op_offset, 0.0165, 0.0165));
         if (myid==0) GbSystem3D::writeGeoObject(refineBoxL6.get(), pathOut+"/geo/refineBoxL6", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbObject3D> refineBoxL7(new GbCuboid3D(startX1it, -0.024, -0.024, startX1it+0.09+op_offset, 0.024, 0.024));
         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL7.get(), pathOut + "/geo/refineBoxL7", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbObject3D> refineBoxL81(new GbCuboid3D(startX1it, -0.005, -0.005, startX1it+0.02+op_offset, 0.005, 0.005));
         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL81.get(), pathOut + "/geo/refineBoxL81", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbObject3D> refineBoxL82(new GbCuboid3D(startX1it+0.02, -0.0165, -0.0165, startX1it+0.06+op_offset, 0.0165, 0.0165));
         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL82.get(), pathOut + "/geo/refineBoxL82", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbObject3D> refineBoxL83(new GbCuboid3D(startX1it+0.06, -0.0165, -0.0165, startX1it+0.09+op_offset, 0.0165, 0.024));
         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL83.get(), pathOut + "/geo/refineBoxL83", WbWriterVtkXmlBinary::getInstance());

         SPtr<GbObject3D> refineBoxL9(new GbCuboid3D(startX1it + 0.06, -0.0165, 0.01, startX1it + 0.09 + op_offset, 0.0165, 0.013));
         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL9.get(), pathOut + "/geo/refineBoxL9", WbWriterVtkXmlBinary::getInstance());

         if (refineLevel>0)
         {
            if (myid==0) UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
            refineHelper.addGbObject(refineBoxL6, refineLevel-3);
            refineHelper.addGbObject(refineBoxL7, refineLevel-2);
            refineHelper.addGbObject(refineBoxL81, refineLevel-1);
            refineHelper.addGbObject(refineBoxL82, refineLevel-1);
            refineHelper.addGbObject(refineBoxL83, refineLevel-1);
            refineHelper.addGbObject(refineBoxL9, refineLevel);
            refineHelper.refine();
            if (myid==0) UBLOG(logINFO, "Refinement - end");
         }
         //////////////////////////////////////////////////////////////////////////

         //walls
         double blockLength = blocknx[0]*deltaXcoarse;
         GbCuboid3DPtr addWallXmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallXmin.get(), pathOut + "/geo/addWallXmin", WbWriterVtkXmlASCII::getInstance());
         GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_minX2, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathOut + "/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());
         GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1-blockLength, g_maxX2, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathOut + "/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());
         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathOut + "/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());
         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathOut + "/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //wall interactors
         SPtr<D3Q27Interactor> addWallXminInt(new D3Q27Interactor(addWallXmin, grid, slipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, slipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, slipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, slipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, slipBCAdapter, Interactor3D::SOLID));

         //inflow
         GbCylinder3DPtr geoInflow(new GbCylinder3D(startX1it-deltaXfine*3.0, 0.0, 0.0, startX1it+deltaXfine*3.0, 0.0, 0.0, diameter_inlet));
         if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathOut + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());
         SPtr<D3Q27Interactor> inflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         GbCylinder3DPtr geoInflowCover(new GbCylinder3D(startX1it-deltaXfine*5.0, 0.0, 0.0, startX1it, 0.0, 0.0, diameter_inlet+deltaXfine*2.0));
         if (myid == 0) GbSystem3D::writeGeoObject(geoInflowCover.get(), pathOut + "/geo/geoInflowCover", WbWriterVtkXmlASCII::getInstance());
         SPtr<D3Q27Interactor> inflowCoverIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflowCover, grid, noSlipBCAdapter, Interactor3D::SOLID));

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - deltaXcoarse, g_minX3 - deltaXcoarse, g_maxX1 + deltaXcoarse, g_maxX2 + deltaXcoarse, g_maxX3 + deltaXcoarse));
         if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());
         SPtr<D3Q27Interactor> outflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowBCAdapter, Interactor3D::SOLID));


         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, MetisPartitioner::RECURSIVE));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(outflowIntr);
         intHelper.addInteractor(addWallXminInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(inflowIntr);
         intHelper.addInteractor(organPipeInter);
         intHelper.addInteractor(inletTubeInter);
         intHelper.addInteractor(inflowCoverIntr);
         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");

         {
            WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
            ppblocks.process(0);
         }

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

         SetKernelBlockVisitor kernelVisitor(kernel, nu_LB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (myid == 0) UBLOG(logINFO, "SetKernelBlockVisitor:end");

         if (myid == 0)
         {
            UBLOG(logINFO, "PID = " << myid << " Point 5");
            UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
         }

         if (refineLevel > 0)
         {
            bool twoTypeOfConnectorsCheck = false;
            SetUndefinedNodesBlockVisitor undefNodesVisitor(twoTypeOfConnectorsCheck);
            grid->accept(undefNodesVisitor);
         }

         if (myid == 0) UBLOG(logINFO, "SetUndefinedNodesBlockVisitor:end");

         //BC
         intHelper.setBC();
         if (myid == 0) UBLOG(logINFO, "intHelper.setBC():end");

         //initialization of distributions
         InitDistributionsBlockVisitor initVisitor;
         grid->accept(initVisitor);

         //Geometry and boundary conditions
         {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
            ppgeo.process(0);
         }

      }
      else
      {
         //restartCoProcessor->restart((int)restartStep);
         migCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);
      }
      ////set connectors
      //InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
      SPtr<InterpolationProcessor> iProcessor(new CompressibleOffsetMomentsInterpolationProcessor());
      dynamicPointerCast<CompressibleOffsetMomentsInterpolationProcessor>(iProcessor)->setBulkViscosity(nu_LB, bulckViscosity);
      SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
      grid->accept(setConnsVisitor);

      //bcVisitor should be accept after initialization!!!!
      grid->accept(bcVisitor);
      if (myid == 0) UBLOG(logINFO, "grid->accept(bcVisitor):end");



      ////sponge layer
      ////////////////////////////////////////////////////////////////////////////
      /*GbCuboid3DPtr spongeLayerX1max(new GbCuboid3D(g_maxX1 - 0.4, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
      if (myid == 0) GbSystem3D::writeGeoObject(spongeLayerX1max.get(), pathOut + "/geo/spongeLayerX1max", WbWriterVtkXmlASCII::getInstance());
      SpongeLayerBlockVisitor slVisitorX1max(spongeLayerX1max, kernel, nu_LB, D3Q27System::E);
      grid->accept(slVisitorX1max);

      GbCuboid3DPtr spongeLayerX1min(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_minX1 + 0.2, g_maxX2 + blockLength, g_maxX3 + blockLength));
      if (myid == 0) GbSystem3D::writeGeoObject(spongeLayerX1min.get(), pathOut + "/geo/spongeLayerX1min", WbWriterVtkXmlASCII::getInstance());
      SpongeLayerBlockVisitor slVisitorX1min(spongeLayerX1min, spKernel, nuLB, D3Q27System::W);
      grid->accept(slVisitorX1min);

      GbCuboid3DPtr spongeLayerX2max(new GbCuboid3D(g_maxX1 - 0.4, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
      if (myid == 0) GbSystem3D::writeGeoObject(spongeLayerX1max.get(), pathOut + "/geo/spongeLayerX1max", WbWriterVtkXmlASCII::getInstance());
      SpongeLayerBlockVisitor slVisitorX1max(spongeLayerX1max, kernel, nu_LB, D3Q27System::E);
      grid->accept(slVisitorX1max);

      GbCuboid3DPtr spongeLayerX2min(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_minX1 + 0.2, g_maxX2 + blockLength, g_maxX3 + blockLength));
      if (myid == 0) GbSystem3D::writeGeoObject(spongeLayerX1min.get(), pathOut + "/geo/spongeLayerX1min", WbWriterVtkXmlASCII::getInstance());
      SpongeLayerBlockVisitor slVisitorX1min(spongeLayerX1min, spKernel, nuLB, D3Q27System::W);
      grid->accept(slVisitorX1min);

      GbCuboid3DPtr spongeLayerX3min(new GbCuboid3D(g_minX1 + 0.2, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 - 0.4, g_maxX2 + blockLength, g_minX3 + 0.2));
      if (myid == 0) GbSystem3D::writeGeoObject(spongeLayerX3min.get(), pathOut + "/geo/spongeLayerX3min", WbWriterVtkXmlASCII::getInstance());
      SpongeLayerBlockVisitor slVisitorX3min(spongeLayerX3min, kernel, nuLB, D3Q27System::B);
      grid->accept(slVisitorX3min);

      GbCuboid3DPtr spongeLayerX3max(new GbCuboid3D(g_minX1 + 0.2, g_minX2 - blockLength, g_maxX3 - 0.2, g_maxX1 - 0.4, g_maxX2 + blockLength, g_maxX3 + blockLength));
      if (myid == 0) GbSystem3D::writeGeoObject(spongeLayerX3max.get(), pathOut + "/geo/spongeLayerX3max", WbWriterVtkXmlASCII::getInstance());
      SpongeLayerBlockVisitor slVisitorX3max(spongeLayerX3max, kernel, nuLB, D3Q27System::T);
      grid->accept(slVisitorX3max);*/
      /////////////////////////////////////////////////////////////////////////////

      if (myid==0) UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      if (myid==0) UBLOG(logINFO, "Preprozess - end");

      int numOfThreads = 1;
      SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      std::shared_ptr<NUPSCounterCoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> stepSch(new UbScheduler(outTimeStep, outTimeStart));

      SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      SPtr<UbScheduler> stepMV(new UbScheduler(1, 0, 1000000));
      SPtr<MicrophoneArrayCoProcessor> micCoProcessor(new MicrophoneArrayCoProcessor(grid, stepSch, pathOut, comm));
      std::vector<UbTupleFloat3> nodes;
      micCoProcessor->addMicrophone(Vector3D(1.43865003014, 0.0, organPipeGeo->getX3Maximum()+0.05));
      nodes.push_back(UbTupleFloat3(float(1.43865003014), float(0.0), float(organPipeGeo->getX3Maximum()+0.05)));
      micCoProcessor->addMicrophone(Vector3D(organPipeGeo->getX1Maximum()+0.05, 0.0, organPipeGeo->getX3Centroid()));
      nodes.push_back(UbTupleFloat3(float(organPipeGeo->getX1Maximum()+0.05), float(0.0), float(organPipeGeo->getX3Centroid())));
      if (myid==0) WbWriterVtkXmlBinary::getInstance()->writeNodes(pathOut+"/geo/mic", nodes);

      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(micCoProcessor);
      calculator->addCoProcessor(migCoProcessor);
      calculator->addCoProcessor(writeMQCoProcessor);
      /////////////////////////////////////////////////////////////////////////////////////

      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      if (myid == 0) UBLOG(logINFO, "Simulation-end");

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

//////////////////////////////////////////////////////////////////////////
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
}
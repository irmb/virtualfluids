#include <iostream>
#include <string>
#include "VirtualFluids.h"


using namespace std;
void run(string configname)
{
   SPtr<Communicator> comm = MPICommunicator::getInstance();
   int myid = comm->getProcessID();

   if (myid == 0) UBLOG(logINFO, "Testcase organ pipe");

   SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

   ConfigurationFile   config;
   config.load(configname);

   string          pathOut = config.getValue<string>("pathOut");
   string          pathGeo = config.getValue<string>("pathGeo");

   double  deltaXfine = 0.000125; //0.000036*2.0;
   const int baseLevel = 0;
   int refineLevel = 9;
   double deltaXcoarse = deltaXfine*(double)(1 << refineLevel);

   double availMem = 5e9;

   //string pathOut = "e:/temp/OrganPipe";
   //string pathGeo = "C:/Users/maini/Desktop/organflute";
   
   string opipeGeoFile = "/02_organf_scaled.stl";
   string inPpipeGeoFile = "/pipeScaled.stl";

   LBMReal rho_LB = 0.0;
   double rhoReal = 1.2041; //(kg/m3)
   double uReal = 3.0; //m/s
   double L = 0.195; //m
   double hLB = L / deltaXcoarse;
   double csReal = 343.3;
   double nuReal = 1.51e-5; //m^2/s
   LBMUnitConverter unitConverter(L, csReal, rhoReal, hLB);
   if (myid == 0) UBLOG(logINFO, unitConverter.toString());

   bool logToFile = false;

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

   vector<int> blocknx = { 16, 16, 16 };

   SPtr<Grid3D> grid(new Grid3D(comm));

   if (myid == 0) UBLOG(logINFO, "Read organ pipe geometry:start");
   SPtr<GbTriFaceMesh3D> opipeGeo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo + opipeGeoFile, "opipeGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
   if (myid == 0) UBLOG(logINFO, "Read organ pipe geometry:end");
   opipeGeo->rotate(0, 270, 0);
   opipeGeo->rotate(0, 0, 90);
   opipeGeo->translate(0.0, 0.0, -(0.0952+0.0165));
   if (myid == 0) GbSystem3D::writeGeoObject(opipeGeo.get(), pathOut + "/geo/opipeGeo", WbWriterVtkXmlBinary::getInstance());

   if (myid == 0) UBLOG(logINFO, "Read inlet pipe geometry:start");
   SPtr<GbTriFaceMesh3D> inPpipeGeo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo + inPpipeGeoFile, "inPipeGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
   if (myid == 0) UBLOG(logINFO, "Read inlet pipe geometry:end");
   if (myid == 0) GbSystem3D::writeGeoObject(inPpipeGeo.get(), pathOut + "/geo/inPpipeGeo", WbWriterVtkXmlBinary::getInstance());

   double offs1 = opipeGeo->getX1Minimum();

   //bounding box
   double g_minX1 = offs1;
   double g_minX2 = -2.048;
   double g_minX3 = -2.048;

   double g_maxX1 = 4.096 + offs1;
   double g_maxX2 = 2.048;
   double g_maxX3 = 2.048;
   double radius_inlet = 0.0055;

   if (myid == 0)
   {
      UBLOG(logINFO, "Parameters:");
      UBLOG(logINFO, "u_LB                = " << u_LB);
      UBLOG(logINFO, "rho_LB              = " << rho_LB);
      UBLOG(logINFO, "nu_LB               = " << nu_LB);
      UBLOG(logINFO, "dx coarse           = " << deltaXcoarse);
      UBLOG(logINFO, "dx fine             = " << deltaXfine);
      UBLOG(logINFO, "number of refinement levels    = " << refineLevel);
      UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
      UBLOG(logINFO, "path = " << pathOut);
      UBLOG(logINFO, "Preprocess - start");
   }

   //BC adapters
   SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
   noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

   SPtr<BCAdapter> slipBCAdapter(new SlipBCAdapter());
   slipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SlipBCAlgorithm()));

   mu::Parser fct;
   fct.SetExpr("U");
   fct.DefineConst("U", u_LB);

   SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
   velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

   SPtr<BCAdapter> outflowBCAdapter(new DensityBCAdapter(rho_LB));
   outflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

   BoundaryConditionsBlockVisitor bcVisitor;
   bcVisitor.addBC(noSlipBCAdapter);
   bcVisitor.addBC(velBCAdapter);
   bcVisitor.addBC(outflowBCAdapter);

   SPtr<BCProcessor> bcProc;
   bcProc = SPtr<BCProcessor>(new BCProcessor());

   SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());

   kernel->setBCProcessor(bcProc);
   //////////////////////////////////////////////////////////////////////////
   //restart
   double          cpStep = 100;
   double          cpStart = 100;
   SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
   SPtr<MPIIOMigrationCoProcessor> migCoProcessor(new MPIIOMigrationCoProcessor(grid, mSch, pathOut + "/mig", comm));
   migCoProcessor->setLBMKernel(kernel);
   migCoProcessor->setBCProcessor(bcProc);
   //////////////////////////////////////////////////////////////////////////

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

   //////////////////////////////////////////////////////////////////////////
   //refinement

   ////Sphere of Radius 3L  
   //SPtr<GbSphere3D> refineSphereL5(new GbSphere3D(g_minX1, 0.0, 0.0, 3.0*L));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineSphereL5.get(), pathOut + "/geo/refineSphereL5", WbWriterVtkXmlASCII::getInstance());

   //SPtr<GbObject3D> refineBoxL6(new GbCuboid3D(-0.12, -0.023, -0.023, 0.14, 0.023, 0.023));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL6.get(), pathOut + "/geo/refineBoxL6", WbWriterVtkXmlASCII::getInstance());

   ////full pipe refine 
   //SPtr<GbObject3D> refineBoxL7(new GbCuboid3D(-0.1107, -0.019, -0.019, 0.1342, 0.019, 0.019));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL7.get(), pathOut + "/geo/refineBoxL7", WbWriterVtkXmlASCII::getInstance());

   ////refine languid and upperlip
   //SPtr<GbObject3D> refineBoxL9(new GbCuboid3D(-0.1007, -0.0115, -0.0115, -0.0634, 0.0115, 0.01255));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL9.get(), pathOut + "/geo/refineBoxL9", WbWriterVtkXmlASCII::getInstance());
  
   //SPtr<GbSphere3D> refineSphereL9(new GbSphere3D(g_minX1+75e-3, 0.0, 0.015, 13e-3));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineSphereL9.get(), pathOut + "/geo/refineSphereL9", WbWriterVtkXmlASCII::getInstance());

   //SPtr<GbObject3D> refineBoxL1(new GbCuboid3D(-0.14, -2.048, -2.048, 3.9643, 2.048, 2.048));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL1.get(), pathOut + "/geo/refineBoxL1", WbWriterVtkXmlASCII::getInstance());

   //SPtr<GbObject3D> refineBoxL2(new GbCuboid3D(-0.14, -0.683, -0.683, 3.32, 0.683, 0.683));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL2.get(), pathOut + "/geo/refineBoxL2", WbWriterVtkXmlASCII::getInstance());

   //SPtr<GbObject3D> refineBoxL3(new GbCuboid3D(-0.14, -0.6, -0.6, 2.5, 0.6, 0.6));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL3.get(), pathOut + "/geo/refineBoxL3", WbWriterVtkXmlASCII::getInstance());

   //SPtr<GbObject3D> refineBoxL4(new GbCuboid3D(-0.14, -0.4315, -0.4315, 2.05, 0.4315, 0.4315));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL4.get(), pathOut + "/geo/refineBoxL4", WbWriterVtkXmlASCII::getInstance());

   //SPtr<GbObject3D> refineBoxL5(new GbCuboid3D(-0.1317, -0.0635, -0.0635, 0.1815, 0.0635, 0.0635));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL5.get(), pathOut + "/geo/refineBoxL5", WbWriterVtkXmlASCII::getInstance());

   //SPtr<GbObject3D> refineBoxL6(new GbCuboid3D(-0.1317, -0.0235, -0.0235, 0.1442, 0.0235, 0.028));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL6.get(), pathOut + "/geo/refineBoxL6", WbWriterVtkXmlASCII::getInstance());

   SPtr<GbObject3D> refineBoxL7(new GbCuboid3D(-0.1317, -0.027, -0.027, 0.1322, 0.027, 0.027));
   if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL7.get(), pathOut + "/geo/refineBoxL7", WbWriterVtkXmlASCII::getInstance());

   SPtr<GbObject3D> refineBoxL8(new GbCuboid3D(-0.1317, -0.022, -0.022, -0.04355, 0.022, 0.022));
   if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL8.get(), pathOut + "/geo/refineBoxL8", WbWriterVtkXmlASCII::getInstance());

   SPtr<GbObject3D> refineBoxL93(new GbCuboid3D(-0.0623, -0.0165, -0.0165, -0.0473, 0.0165, 0.02125));
   if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL93.get(), pathOut + "/geo/refineBoxL93", WbWriterVtkXmlASCII::getInstance());

   SPtr<GbObject3D> refineBoxL92(new GbCuboid3D(-0.1107, -0.0165, -0.0165, -0.0623, 0.0165, 0.0165 ));
   if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL92.get(), pathOut + "/geo/refineBoxL92", WbWriterVtkXmlASCII::getInstance());

   SPtr<GbObject3D> refineBoxL91(new GbCuboid3D(-0.1317, -0.0062, -0.0062, -0.1107, 0.0062, 0.0062));
   if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL91.get(), pathOut + "/geo/refineBoxL91", WbWriterVtkXmlASCII::getInstance());

   //SPtr<GbObject3D> refineBoxL10(new GbCuboid3D(-0.14, -0.0026, -0.0026, -0.12, 0.0026, 0.0026));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL10.get(), pathOut + "/geo/refineBoxL10", WbWriterVtkXmlASCII::getInstance());

   if (refineLevel > 0)
   {
      if (myid == 0) UBLOG(logINFO, "Refinement - start");
      RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
      //refineHelper.addGbObject(refineBoxL3, refineLevel - 6);
      //refineHelper.addGbObject(refineBoxL4, refineLevel - 5);
      //refineHelper.addGbObject(refineBoxL5, refineLevel - 4);
      //refineHelper.addGbObject(refineBoxL6, refineLevel - 3);
      refineHelper.addGbObject(refineBoxL7, refineLevel - 2);
      refineHelper.addGbObject(refineBoxL8, refineLevel - 1);
      refineHelper.addGbObject(refineBoxL93, refineLevel);
      refineHelper.addGbObject(refineBoxL92, refineLevel);
      refineHelper.addGbObject(refineBoxL91, refineLevel);
      //refineHelper.addGbObject(refineBoxL10, refineLevel);
      refineHelper.refine();
      if (myid == 0) UBLOG(logINFO, "Refinement - end");
   }
   //////////////////////////////////////////////////////////////////////////

   {
      WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
      ppblocks.process(0);
   }

   //interactors
   SPtr<Interactor3D> opipeInter = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(opipeGeo, grid, noSlipBCAdapter, Interactor3D::SOLID));
   SPtr<Interactor3D> InPipeInter = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(inPpipeGeo, grid, noSlipBCAdapter, Interactor3D::SOLID));
  
   //walls
   GbCuboid3DPtr addWallXmin(new GbCuboid3D(g_minX1 - deltaXcoarse, g_minX2 - deltaXcoarse, g_minX3 - deltaXcoarse, g_minX1, g_maxX2 + deltaXcoarse, g_maxX3 + deltaXcoarse));
   if (myid == 0) GbSystem3D::writeGeoObject(addWallXmin.get(), pathOut + "/geo/addWallXmin", WbWriterVtkXmlASCII::getInstance());
   GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1-0.001 , g_minX2-0.001, g_minX3-0.001, g_maxX1+0.001, g_minX2, g_maxX3+0.001));
   if (myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathOut + "/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());
   GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1-0.001, g_maxX2, g_minX3-0.001, g_maxX1+0.001, g_maxX2+0.001, g_maxX3+0.001));
   if (myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathOut + "/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());
   GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-0.001, g_minX2-0.001, g_minX3-0.001, g_maxX1+0.001, g_maxX2+0.001, g_minX3));
   if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathOut + "/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());
   GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-0.001, g_minX2-0.001, g_maxX3, g_maxX1+0.001, g_maxX2+0.001, g_maxX3+0.001));
   if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathOut + "/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

   //wall interactors
   SPtr<D3Q27Interactor> addWallXminInt(new D3Q27Interactor(addWallXmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
   SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
   SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
   SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
   SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));

   //inflow
   GbCylinder3DPtr geoInflow(new GbCylinder3D(g_minX1- deltaXcoarse, 0.0, 0.0, g_minX1, 0.0, 0.0, radius_inlet));
   if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathOut + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());
   //outflow
   GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - deltaXcoarse, g_minX3 - deltaXcoarse, g_maxX1 + deltaXcoarse, g_maxX2 + deltaXcoarse, g_maxX3 + deltaXcoarse));
   if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());
   //inflow
   SPtr<D3Q27Interactor> inflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));
   //outflow
   SPtr<D3Q27Interactor> outflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowBCAdapter, Interactor3D::SOLID));

   ////////////////////////////////////////////
   //METIS
   SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
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
   intHelper.addInteractor(opipeInter);
   intHelper.addInteractor(InPipeInter);
   intHelper.addInteractor(inflowIntr);

   intHelper.selectBlocks();
   if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");

   {
      WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
      ppblocks.process(1);
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
      SetUndefinedNodesBlockVisitor undefNodesVisitor;
      grid->accept(undefNodesVisitor);
   }

   if (myid == 0) UBLOG(logINFO, "SetUndefinedNodesBlockVisitor:end");

   //BC
   intHelper.setBC();
   if (myid == 0) UBLOG(logINFO, "intHelper.setBC():end");

   //initialization of distributions
   InitDistributionsBlockVisitor initVisitor1;
   //initVisitor1.setVx1(fct);
   grid->accept(initVisitor1);

   ////set connectors
   InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
   SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
   grid->accept(setConnsVisitor);

   //bcVisitor should be accept after initialization!!!!
   grid->accept(bcVisitor);
   if (myid == 0) UBLOG(logINFO, "grid->accept(bcVisitor):end");

   //Post process
   {
      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
      ppgeo.process(0);
   }

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


   vector<double>  nupsStep = { 10,10,100000 };
   int numOfThreads = 1;
   double outTimeStep = 100;
   double outTimeStart = 100;
   double endTime = 100;

   SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
   std::shared_ptr<NUPSCounterCoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

   SPtr<UbScheduler> stepSch(new UbScheduler(outTimeStep, outTimeStart));

   SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));

   SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
   SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
   calculator->addCoProcessor(nupsCoProcessor);
   calculator->addCoProcessor(migCoProcessor);
   calculator->addCoProcessor(writeMQCoProcessor);
   /////////////////////////////////////////////////////////////////////////////////////

   //SPtr<IntegrateValuesHelper> mic1(new IntegrateValuesHelper(grid, comm, - 0.0598, -0.004, -0.004, 0.1857, 0.004, 0.004));
   //SPtr<IntegrateValuesHelper> mic2(new IntegrateValuesHelper(grid, comm, 0.1777, -0.004, -0.004, 0.1857, 0.004, 0.004));
   ////for check in paraview
   //if (myid == 0) GbSystem3D::writeGeoObject(mic1->getBoundingBox().get(),pathOut + "/geo/mic1", WbWriterVtkXmlBinary::getInstance());
   //if (myid == 0) GbSystem3D::writeGeoObject(mic2->getBoundingBox().get(), pathOut + "/geo/mic2", WbWriterVtkXmlBinary::getInstance());
   //SPtr<UbScheduler> stepMV(new UbScheduler(1, 0, 1000000));
   //SPtr<TimeseriesCoProcessor> tsp1(new TimeseriesCoProcessor(grid, stepMV, mic1, pathOut + "/mic/mic1", comm));
   //SPtr<TimeseriesCoProcessor> tsp2(new TimeseriesCoProcessor(grid, stepMV, mic1, pathOut + "/mic/mic2", comm));
   //calculator->addCoProcessor(tsp1);
   //calculator->addCoProcessor(tsp2);
   if (myid == 0) UBLOG(logINFO, "Simulation-start");
   calculator->calculate();
   if (myid == 0) UBLOG(logINFO, "Simulation-end");

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
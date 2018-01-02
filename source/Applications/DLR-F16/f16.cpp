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

void setBC(Grid3DPtr grid, string pathGeo, string fngFileWhole, string zigZagTape, vector<double>  boundingBox, double uLB, double rhoLB, double blockLength, BCProcessorPtr bcProcessor)
{
   CommunicatorPtr comm = MPICommunicator::getInstance();
   int myid = comm->getProcessID();
   
   std::vector<std::vector<Block3DPtr> > blockVector;
   int minInitLevel;
   int maxInitLevel;
   int gridRank;

   gridRank = comm->getProcessID();
   minInitLevel = grid->getCoarsestInitializedLevel();
   maxInitLevel = grid->getFinestInitializedLevel();
   blockVector.resize(maxInitLevel+1);
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      grid->getBlocks(level, gridRank, true, blockVector[level]);
   }

   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      BOOST_FOREACH(Block3DPtr block, blockVector[level])
      {
         if (block)
         {
            LBMKernelPtr kernel = block->getKernel();
            kernel->setBCProcessor(bcProcessor->clone(kernel));
         }
      }
   }
   
   SetUndefinedNodesBlockVisitor undefNodesVisitor;
   grid->accept(undefNodesVisitor);
   
   BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
   noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithm()));
   if (myid==0) UBLOG(logINFO, "Read fngFileWhole:start");
   GbTriFaceMesh3DPtr fngMeshWhole = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileWhole, "fngMeshWhole", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
   if (myid==0) UBLOG(logINFO, "Read fngFileWhole:end");
   fngMeshWhole->rotate(0.0, 0.5, 0.0);
   D3Q27TriFaceMeshInteractorPtr fngIntrWhole = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshWhole, grid, noSlipBCAdapter, Interactor3D::SOLID));

   if (myid==0) UBLOG(logINFO, "Read zigZagTape:start");
   string ZckbndFilename = pathGeo+"/"+zigZagTape;
   GbTriFaceMesh3DPtr meshBand1(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape1"));
   meshBand1->rotate(0.0, 5, 0.0);
   meshBand1->translate(15, 0, -12.850);
   // Zackenband2
   GbTriFaceMesh3DPtr meshBand2(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape2"));
   meshBand2->rotate(0.0, 5, 0.0);
   meshBand2->translate(15, 5, -12.850);

   GbTriFaceMesh3DPtr meshBand5(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape5"));
   meshBand5->rotate(0.0, -1, 0.0);
   meshBand5->rotate(0.0, 0.0, 180.0);
   //meshBand5->translate(30, 0, -37.3);
   meshBand5->translate(30, 0, -37.2);
   
   // Zackenband6
   GbTriFaceMesh3DPtr meshBand6(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape6"));
   meshBand6->rotate(0.0, -1, 0.0);
   meshBand6->rotate(0.0, 0.0, 180.0);
   //meshBand6->translate(30, 5, -37.3);
   meshBand6->translate(30, 5, -37.2);

   D3Q27TriFaceMeshInteractorPtr triBand1Interactor(new D3Q27TriFaceMeshInteractor(meshBand1, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
   D3Q27TriFaceMeshInteractorPtr triBand2Interactor(new D3Q27TriFaceMeshInteractor(meshBand2, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
   D3Q27TriFaceMeshInteractorPtr triBand3Interactor(new D3Q27TriFaceMeshInteractor(meshBand5, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
   D3Q27TriFaceMeshInteractorPtr triBand4Interactor(new D3Q27TriFaceMeshInteractor(meshBand6, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));

   mu::Parser fct;
   fct.SetExpr("U");
   fct.DefineConst("U", uLB);
   BCAdapterPtr velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
   velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityWithDensityBCAlgorithm()));

   BCAdapterPtr outflowBCAdapter(new DensityBCAdapter(rhoLB));
   outflowBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithm()));

   double g_minX1 = boundingBox[0]*1000.0;
   double g_minX2 = boundingBox[2]*1000.0;
   double g_minX3 = boundingBox[4]*1000.0;

   double g_maxX1 = boundingBox[1]*1000.0;
   double g_maxX2 = boundingBox[3]*1000.0;
   double g_maxX3 = boundingBox[5]*1000.0;
   
   GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
   
   GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
   
   D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, velBCAdapter, Interactor3D::SOLID));
   D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, velBCAdapter, Interactor3D::SOLID));
   //inflow
   GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
   //outflow
   GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
   //inflow
   D3Q27InteractorPtr inflowIntr = D3Q27InteractorPtr(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));
   //outflow
   D3Q27InteractorPtr outflowIntr = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, outflowBCAdapter, Interactor3D::SOLID));

   SetSolidBlockVisitor v1(fngIntrWhole, SetSolidBlockVisitor::BC);
   grid->accept(v1);

   SetSolidBlockVisitor v2(triBand1Interactor, SetSolidBlockVisitor::BC);
   grid->accept(v2);

   SetSolidBlockVisitor v3(triBand1Interactor, SetSolidBlockVisitor::BC);
   grid->accept(v3);

   SetSolidBlockVisitor v4(triBand2Interactor, SetSolidBlockVisitor::BC);
   grid->accept(v4);
   
   SetSolidBlockVisitor v5(triBand3Interactor, SetSolidBlockVisitor::BC);
   grid->accept(v5);

   SetSolidBlockVisitor v6(triBand4Interactor, SetSolidBlockVisitor::BC);
   grid->accept(v6);

   SetSolidBlockVisitor v7(addWallZminInt, SetSolidBlockVisitor::BC);
   grid->accept(v7);

   SetSolidBlockVisitor v8(addWallZmaxInt, SetSolidBlockVisitor::BC);
   grid->accept(v8);

   SetSolidBlockVisitor v9(inflowIntr, SetSolidBlockVisitor::BC);
   grid->accept(v9);

   SetSolidBlockVisitor v10(outflowIntr, SetSolidBlockVisitor::BC);
   grid->accept(v10);
   
   inflowIntr->initInteractor();
   outflowIntr->initInteractor();
   addWallZminInt->initInteractor();
   addWallZmaxInt->initInteractor();
   fngIntrWhole->initInteractor();
   triBand1Interactor->initInteractor();
   triBand2Interactor->initInteractor();
   triBand3Interactor->initInteractor();
   triBand3Interactor->initInteractor();
   triBand4Interactor->initInteractor();
}

void run(string configname)
{
   try
   {

      ConfigurationFile   config;
      config.load(configname);

      string          pathOut = config.getValue<string>("pathOut");
      string          pathGeo = config.getValue<string>("pathGeo");
      string          fngFileWhole1 = config.getValue<string>("fngFileWhole1");
      string          fngFileWhole2 = config.getValue<string>("fngFileWhole2");
      string          fngFileTrailingEdge = config.getValue<string>("fngFileTrailingEdge");
      string          fngFileBodyPart = config.getValue<string>("fngFileBodyPart");
      string          zigZagTape = config.getValue<string>("zigZagTape");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blockNx = config.getVector<int>("blockNx");
      vector<double>  boundingBox = config.getVector<double>("boundingBox");
      double          uLB = config.getValue<double>("uLB");
      double          restartStep = config.getValue<double>("restartStep");
      double          cpStart = config.getValue<double>("cpStart");
      double          cpStep = config.getValue<double>("cpStep");
      double          endTime = config.getValue<double>("endTime");
      double          outTimeStep = config.getValue<double>("outTimeStep");
      double          outTimeStart = config.getValue<double>("outTimeStart");
      double          availMem = config.getValue<double>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      bool            logToFile = config.getValue<bool>("logToFile");
      bool            porousTralingEdge = config.getValue<bool>("porousTralingEdge");
      double          deltaXfine = config.getValue<double>("deltaXfine")*1000.0;
      bool            thinWall = config.getValue<bool>("thinWall");
      double          refineDistance = config.getValue<double>("refineDistance");
      double          startDistance = config.getValue<double>("startDistance");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");
      bool            newStart = config.getValue<bool>("newStart");
      bool            writeBlocks = config.getValue<bool>("writeBlocks");
      string          sampleFilename = config.getValue<string>("sampleFilename");
      string          pathReInit = config.getValue<string>("pathReInit");
      int             stepReInit = config.getValue<int>("stepReInit");
      
      double          pcpStart = config.getValue<double>("pcpStart");
      double          pcpStop  = config.getValue<double>("pcpStop");
      //double          p_inf    = config.getValue<double>("p_inf");
      
      double          timeAvStart       = config.getValue<double>("timeAvStart");
      double          timeAvStop        = config.getValue<double>("timeAvStop");
      
      int             chunk = config.getValue<int>("chunk");


      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (logToFile)
      {
#if defined(__unix__)
         if (myid==0)
         {
            const char* str = pathOut.c_str();
            mkdir(str, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
         }
#endif 

         if (myid==0)
         {
            stringstream logFilename;
            logFilename<<pathOut+"/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 1");
         UBLOG(logINFO, "PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }


      //the geometry is in mm

      double g_minX1 = boundingBox[0]*1000.0;
      double g_minX2 = boundingBox[2]*1000.0;
      double g_minX3 = boundingBox[4]*1000.0;

      double g_maxX1 = boundingBox[1]*1000.0;
      double g_maxX2 = boundingBox[3]*1000.0;
      double g_maxX3 = boundingBox[5]*1000.0;

      //////////////////////////////////////////////////////////////////////////
      double deltaXcoarse = deltaXfine*(double)(1<<refineLevel);
      //double nx2_temp = floor((g_maxX2 - g_minX2) / (deltaXcoarse*(double)blockNx[0]));

      //deltaXcoarse = (g_maxX2 - g_minX2) / (nx2_temp*(double)blockNx[0]);
      //UBLOG(logINFO, "nx2_temp:"<<nx2_temp);
      //g_maxX2 -= 0.5* deltaXcoarse;
      //////////////////////////////////////////////////////////////////////////
      double blockLength = (double)blockNx[0]*deltaXcoarse;

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
      double lLB = lReal*1000.0/deltaXcoarse;
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
      velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityWithDensityBCAlgorithm()));
      //velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithm()));

      BCAdapterPtr outflowBCAdapter(new DensityBCAdapter(rhoLB));
      outflowBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingOutflowBCAlgorithm()));
      //denBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonEqDensityBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      //bcVisitor.addBC(slipBCAdapter);
      bcVisitor.addBC(velBCAdapter);
      bcVisitor.addBC(outflowBCAdapter);

      LBMKernelPtr kernel = LBMKernelPtr(new CompressibleCumulant2LBMKernel(blockNx[0], blockNx[1], blockNx[2], CompressibleCumulant2LBMKernel::NORMAL));
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
      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(cpStep, cpStart));
      //MPIIORestartCoProcessor rcp(grid, rSch, pathOut, comm);
      //rcp.setChunk(chunk);
      
      UbSchedulerPtr rSch2(new UbScheduler(restartStep));
      //RestartCoProcessor rp(grid, rSch2, comm, pathOut, RestartCoProcessor::BINARY);
     
      //MPIIORestart2CoProcessor rcp2(grid, rSch2, pathOut+"/mpiio2", comm);      
      //rcp2.setLBMKernel(kernel);
      //rcp2.setBCProcessor(bcProc);
      
      
      MPIIORestart1CoProcessor rcp3(grid, rSch, pathOut+"/mpiio3", comm);
      rcp3.setLBMKernel(kernel);
      rcp3.setBCProcessor(bcProc);
      
      //MPIIORestart3CoProcessor rcp4(grid, rSch, pathOut+"/mpiio4", comm);
      //rcp4.setLBMKernel(kernel);
      //rcp4.setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////


      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 2");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }

      //if (grid->getTimeStep() == 0)
      if (newStart)
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

         if (myid==0)
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
            UBLOG(logINFO, "* blocknx             = "<<blockNx[0]<<"x"<<blockNx[1]<<"x"<<blockNx[2]);
            UBLOG(logINFO, "* refineDistance      = "<<refineDistance);
            UBLOG(logINFO, "* number of levels    = "<<refineLevel+1);
            UBLOG(logINFO, "* number of threads   = "<<numOfThreads);
            UBLOG(logINFO, "* number of processes = "<<comm->getNumberOfProcesses());
            UBLOG(logINFO, "Preprozess - start");
         }

         GbTriFaceMesh3DPtr fngMeshWhole1;
         GbTriFaceMesh3DPtr fngMeshWhole2;
         GbTriFaceMesh3DPtr fngMeshBodyPart;
         GbTriFaceMesh3DPtr fngMeshTrailingEdge;
         GbTriFaceMesh3DPtr porousTrailingEdge;
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
            fngMeshTrailingEdge->translate(-0.05, 0, 1.3);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshTrailingEdge.get(), pathOut+"/geo/fngMeshTrailingEdge", WbWriterVtkXmlBinary::getInstance());

            if (myid==0) UBLOG(logINFO, "Read porousTrailingEdge:start");
            //porousTrailingEdge = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+sampleFilename, "porousTrailingEdge", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            porousTrailingEdge = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+sampleFilename, "porousTrailingEdge", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read porousTrailingEdge:end");
            porousTrailingEdge->rotate(90, -90, 0.0);
            porousTrailingEdge->rotate(0, -4.3, 0.0);
            //porousTrailingEdge->translate(280.5, 40.0, 3.5);
            porousTrailingEdge->translate(276, 15.95, 3.5);
            if (myid==0) GbSystem3D::writeGeoObject(porousTrailingEdge.get(), pathOut+"/geo/porousTrailingEdge", WbWriterVtkXmlASCII::getInstance());

            //string samplePathname = pathGeo+"/"+sampleFilename;

            //int pmNX1 = 669;
            //int pmNX2 = 2945;
            //int pmNX3 = 100;

            //double deltaVoxelX1 = 13393e-6;
            //double deltaVoxelX2 = 13393e-6;
            //double deltaVoxelX3 = 13393e-6;

            //double lthreshold = 11538;
            //double uthreshold = 65535;

            //if (myid==0) UBLOG(logINFO, "read voxel matrix: start");
            //GbVoxelMatrix3DPtr porousTrailingEdge(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0.0 , lthreshold, uthreshold));
            //porousTrailingEdge->readMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
            //porousTrailingEdge->setVoxelMatrixDelta((float)deltaVoxelX1, (float)deltaVoxelX2, (float)deltaVoxelX3);
            //porousTrailingEdge->setVoxelMatrixMininum(0.0, 0.0, 0.0);
            //if (myid==0) UBLOG(logINFO, "read voxel matrix: end");

            //if (myid==0) UBLOG(logINFO, "rotate voxel matrix: start");
            //porousTrailingEdge->rotate90aroundZ();
            //porousTrailingEdge->rotate90aroundZ();
            //porousTrailingEdge->rotate90aroundZ();
            //porousTrailingEdge->rotate90aroundX();
            //porousTrailingEdge->rotateAroundY(0.07);
            //porousTrailingEdge->translate(276, 15.95, 3.26);
            //
            //if (myid==0) UBLOG(logINFO, "rotate voxel matrix: end");

            ////if (myid==0) porousTrailingEdge->writeToVTKImageDataASCII(pathOut+"/geo/PorousTrailingEdge");
            //if (myid==0) porousTrailingEdge->writeToVTKImageDataAppended(pathOut+"/geo/PorousTrailingEdge");

         }
         else
         {
            //if (myid==0) UBLOG(logINFO, "Read fngFileWhole1:start");
            //fngMeshWhole1 = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileWhole1, "fngMeshWhole1", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            //if (myid==0) UBLOG(logINFO, "Read fngFileWhole1:end");
            //fngMeshWhole1->rotate(0.0, 0.5, 0.0);
            //if (myid==0) GbSystem3D::writeGeoObject(fngMeshWhole1.get(), pathOut+"/geo/fngMeshWhole1", WbWriterVtkXmlBinary::getInstance());
            
            if (myid==0) UBLOG(logINFO, "Read fngFileWhole2:start");
            fngMeshWhole2 = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileWhole2, "fngMeshWhole2", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read fngFileWhole2:end");
            fngMeshWhole2->rotate(0.0, 0.5, 0.0);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshWhole2.get(), pathOut+"/geo/fngMeshWhole2", WbWriterVtkXmlBinary::getInstance());
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
         meshBand5->rotate(0.0, 0.0, 180.0);
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


         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 3");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }
         //////////////////////////////////////////////////////////////////////////
         //return;

         //Interactor3DPtr fngIntrWhole1;
         Interactor3DPtr fngIntrWhole2;
         Interactor3DPtr fngIntrBodyPart;
         Interactor3DPtr fngIntrTrailingEdge;
         Interactor3DPtr porousIntrTrailingEdge;

         if (porousTralingEdge)
         {
            fngIntrBodyPart = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshBodyPart, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
            fngIntrTrailingEdge = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshTrailingEdge, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
            porousIntrTrailingEdge = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(porousTrailingEdge, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
         }
         else
         {
            //fngIntrWhole1 = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshWhole1, grid, noSlipBCAdapter, Interactor3D::SOLID));
            fngIntrWhole2 = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshWhole2, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::POINTS));
         }

         D3Q27TriFaceMeshInteractorPtr triBand1Interactor(new D3Q27TriFaceMeshInteractor(meshBand1, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
         D3Q27TriFaceMeshInteractorPtr triBand2Interactor(new D3Q27TriFaceMeshInteractor(meshBand2, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
         D3Q27TriFaceMeshInteractorPtr triBand3Interactor(new D3Q27TriFaceMeshInteractor(meshBand5, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
         D3Q27TriFaceMeshInteractorPtr triBand4Interactor(new D3Q27TriFaceMeshInteractor(meshBand6, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));

         if (refineLevel>0&&myid==0&&writeBlocks)
         {
            if (myid==0) UBLOG(logINFO, "Refinement - start");
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
               boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrBodyPart)->refineBlockGridToLevel(refineLevel-1, startDistance, refineDistance);
            }
            else
            {
               boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrWhole2)->refineBlockGridToLevel(refineLevel, startDistance, refineDistance);
            }

            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand1Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand2Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand3Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);
            //boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(triBand4Interactor)->refineBlockGridToLevel(refineLevel, 0.0, refineDistance);


            //GbObject3DPtr fngBox(new GbCuboid3D(fngMeshWhole->getX1Minimum(), fngMeshWhole->getX2Minimum(), fngMeshWhole->getX3Minimum(),
            //                                    fngMeshWhole->getX1Maximum(), fngMeshWhole->getX2Maximum(), fngMeshWhole->getX3Maximum()));
            //if (myid==0) GbSystem3D::writeGeoObject(fngBox.get(), pathOut+"/geo/fngBox", WbWriterVtkXmlASCII::getInstance());

            //RefineCrossAndInsideGbObjectBlockVisitor refVisitor0(fngBox, refineLevel);
            //grid->accept(refVisitor0);


            //GbObject3DPtr bandTopBox(new GbCuboid3D(meshBand1->getX1Minimum(), meshBand1->getX2Minimum(), meshBand1->getX3Minimum(),
            //   meshBand1->getX1Maximum(), meshBand1->getX2Maximum(), meshBand1->getX3Maximum()));
            //if (myid==0) GbSystem3D::writeGeoObject(bandTopBox.get(), pathOut+"/geo/bandTopBox", WbWriterVtkXmlASCII::getInstance());

            //RefineCrossAndInsideGbObjectBlockVisitor refVisitor1(bandTopBox, refineLevel-1);
            //grid->accept(refVisitor1);

            //GbObject3DPtr bandBottomBox(new GbCuboid3D(meshBand5->getX1Minimum(), meshBand5->getX2Minimum(), meshBand5->getX3Minimum(),
            //   meshBand5->getX1Maximum(), meshBand5->getX2Maximum(), meshBand5->getX3Maximum()));
            //if (myid==0) GbSystem3D::writeGeoObject(bandBottomBox.get(), pathOut+"/geo/bandBottomBox", WbWriterVtkXmlASCII::getInstance());

            //RefineCrossAndInsideGbObjectBlockVisitor refVisitor2(bandBottomBox, refineLevel-1);
            //grid->accept(refVisitor2);

            //GbObject3DPtr teBox1(new GbCuboid3D(269.0, 0.0, 1.0, 270.0, 100.0, 8.5));
            // for porous teY
            //GbObject3DPtr teBox1(new GbCuboid3D(269.0, 0.0, -10.0, 310.0, 100.0, 20.5));
            //GbObject3DPtr teBox1(new GbCuboid3D(200.0, 0.0, -20.0, 400.0, 100.0, 20.0));
            //if (myid==0) GbSystem3D::writeGeoObject(teBox1.get(), pathOut+"/geo/teBox1", WbWriterVtkXmlASCII::getInstance());

            //RefineCrossAndInsideGbObjectBlockVisitor refVisitor3(teBox1, 5);
            //grid->accept(refVisitor3);

            //GbObject3DPtr teBox2(new GbCuboid3D(271.0, 0.0, 3.0, 279.0, 100.0, 5.7));
            //if (myid==0) GbSystem3D::writeGeoObject(teBox2.get(), pathOut+"/geo/teBox2", WbWriterVtkXmlASCII::getInstance());

            //RefineCrossAndInsideGbObjectBlockVisitor refVisitor4(teBox2, refineLevel);
            //grid->accept(refVisitor4);

            //level 1
            GbObject3DPtr wakeBoxL1(new GbCuboid3D(200.0, 0.0, -20.0, 2000.0, 100.0, 20.0));
            if (myid==0) GbSystem3D::writeGeoObject(wakeBoxL1.get(), pathOut+"/geo/wakeBoxL1", WbWriterVtkXmlASCII::getInstance());
            RefineCrossAndInsideGbObjectBlockVisitor refVisitorWakeBoxL1(wakeBoxL1, 1);
            grid->accept(refVisitorWakeBoxL1);
            
            //level 4
            //GbObject3DPtr teBoxL5(new GbCuboid3D(200.0, 0.0, -20.0, 400.0, 100.0, 20.0));
            //if (myid==0) GbSystem3D::writeGeoObject(teBoxL5.get(), pathOut+"/geo/teBoxL5", WbWriterVtkXmlASCII::getInstance());
            //RefineCrossAndInsideGbObjectBlockVisitor refVisitorTeBoxL5(teBoxL5, 4);
            //grid->accept(refVisitorTeBoxL5);
            
            //level 5
            //GbObject3DPtr teBoxL6(new GbCuboid3D(270.0, 0.0, -3.0, 320.0, 100.0, 10.0));
            //if (myid==0) GbSystem3D::writeGeoObject(teBoxL6.get(), pathOut+"/geo/teBoxL6", WbWriterVtkXmlASCII::getInstance());
            //RefineCrossAndInsideGbObjectBlockVisitor refVisitorTeBoxL6(teBoxL6, 5);
            //grid->accept(refVisitorTeBoxL6);

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
            if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");
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
               SetSolidBlockVisitor v(fngIntrBodyPart, SetSolidBlockVisitor::SOLID);
               grid->accept(v);
               std::vector<Block3DPtr>& sb = fngIntrBodyPart->getSolidBlockSet();
               BOOST_FOREACH(Block3DPtr block, sb)
               {
                  grid->deleteBlock(block);
               }
               fngIntrBodyPart->removeSolidBlocks();
               fngIntrBodyPart->removeBcBlocks();
            }
            else
            {
               SetSolidBlockVisitor v(fngIntrWhole2, SetSolidBlockVisitor::SOLID);
               grid->accept(v);
               std::vector<Block3DPtr>& sb = fngIntrWhole2->getSolidBlockSet();
               BOOST_FOREACH(Block3DPtr block, sb)
               {
                  grid->deleteBlock(block);
               }
               fngIntrWhole2->removeSolidBlocks();
               fngIntrWhole2->removeBcBlocks();
            }

            if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - end");
            //////////////////////////////////////

            //if (porousTralingEdge)
            //{
            //   grid->setRank(0);
            //   boost::dynamic_pointer_cast<D3Q27TriFaceMeshInteractor>(fngIntrTrailingEdge)->refineBlockGridToLevel(refineLevel, startDistance, refineDistance);
            //   grid->setRank(rank);

            //   //GbObject3DPtr trailingEdgeCube(new GbCuboid3D(fngMeshTrailingEdge->getX1Minimum()-blockLength, fngMeshTrailingEdge->getX2Minimum(), fngMeshTrailingEdge->getX3Minimum()-blockLength/2.0,
            //   //   fngMeshTrailingEdge->getX1Maximum()+blockLength, fngMeshTrailingEdge->getX2Maximum(), fngMeshTrailingEdge->getX3Maximum()+blockLength/2.0));
            //   //if (myid == 0) GbSystem3D::writeGeoObject(trailingEdgeCube.get(), pathOut + "/geo/trailingEdgeCube", WbWriterVtkXmlASCII::getInstance());

            //   //RefineCrossAndInsideGbObjectBlockVisitor refVisitor(trailingEdgeCube, refineLevel);
            //   //grid->accept(refVisitor);
            //}

            RatioBlockVisitor ratioVisitor(refineLevel);
            CheckRatioBlockVisitor checkRatio(refineLevel);
            int count = 0;

            do {
               grid->accept(ratioVisitor);
               checkRatio.resetState();
               grid->accept(checkRatio);
               if (myid==0) UBLOG(logINFO, "count = "<<count++<<" state = "<<checkRatio.getState());
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

            if (myid==0) UBLOG(logINFO, "Refinement - end");
         }
         grid->updateDistributedBlocks(comm);

         //if (writeBlocks)
         //{
         //   grid->updateDistributedBlocks(comm);
         //   rcp.writeBlocks(0);
         //}
         //else
         //{
           //rcp.readBlocks(restartStep);
           //grid->setTimeStep(restartStep);
         //}

         //return;

         //Sleep(1000*myid);

           


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
         //D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, slipBCAdapter, Interactor3D::SOLID));
         //D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, slipBCAdapter, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, velBCAdapter, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, velBCAdapter, Interactor3D::SOLID));

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathOut+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         //inflow
         D3Q27InteractorPtr inflowIntr = D3Q27InteractorPtr(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //outflow
         D3Q27InteractorPtr outflowIntr = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, outflowBCAdapter, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         //grid->deleteBlockIDs();
         //RenumberBlockVisitor renumber;
         //grid->accept(renumber);
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(inflowIntr);
         intHelper.addInteractor(outflowIntr);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(triBand1Interactor);
         intHelper.addInteractor(triBand2Interactor);
         intHelper.addInteractor(triBand3Interactor);
         intHelper.addInteractor(triBand4Interactor);
         

         if (porousTralingEdge)
         {
            intHelper.addInteractor(fngIntrBodyPart);
            intHelper.addInteractor(porousIntrTrailingEdge);

            //string samplePathname = pathGeo+"/"+sampleFilename;

            //double pmNX1 = 669;
            //double pmNX2 = 2945;
            //double pmNX3 = 1119;

            //double deltaVoxelX1 = 13393e-6;
            //double deltaVoxelX2 = 13393e-6;
            //double deltaVoxelX3 = 13393e-6;
            //
            //double lthreshold = 11538;
            //double uthreshold = 65535;

            //if (myid==0) UBLOG(logINFO, "read voxel matrix: start");
            //GbVoxelMatrix3DPtr porousTrailingEdge(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0, lthreshold, uthreshold));
            //porousTrailingEdge->readMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
            //porousTrailingEdge->setVoxelMatrixDelta((float)deltaVoxelX1, (float)deltaVoxelX2, (float)deltaVoxelX3);
            //porousTrailingEdge->setVoxelMatrixMininum(0.0, 0.0, 0.0);
            //if (myid==0) UBLOG(logINFO, "read voxel matrix: end");

            //if (myid==0) UBLOG(logINFO, "rotate voxel matrix: start");
            //porousTrailingEdge->rotate90aroundZ();
            //porousTrailingEdge->rotate90aroundX();
            //if (myid==0) UBLOG(logINFO, "rotate voxel matrix: end");

            //if (myid==0) porousTrailingEdge->writeToVTKImageDataASCII(pathOut+"/geo/PorousTrailingEdge");
         }
         else
         {
            intHelper.addInteractor(fngIntrWhole2);
         }

         //////////////////////////////////////////////////////////////////////////
         intHelper.selectBlocks();

         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 4");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }
         //////////////////////////////////////

         {
            WriteBlocksCoProcessor ppblocks(grid, UbSchedulerPtr(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
            ppblocks.process(2);
         }

         unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
         int ghostLayer = 3;
         unsigned long long numberOfNodesPerBlock = (unsigned long long)(blockNx[0])* (unsigned long long)(blockNx[1])* (unsigned long long)(blockNx[2]);
         unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
         unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blockNx[0]+ghostLayer) * (blockNx[1]+ghostLayer) * (blockNx[2]+ghostLayer);
         double needMemAll = double(numberOfNodesPerBlockWithGhostLayer*(27*sizeof(double)+sizeof(int)+sizeof(float)*4));
         double needMem = needMemAll/double(comm->getNumberOfProcesses());

         if (myid==0)
         {
            UBLOG(logINFO, "Number of blocks = "<<numberOfBlocks);
            UBLOG(logINFO, "Number of nodes  = "<<numberOfNodes);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for (int level = minInitLevel; level<=maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO, "Number of blocks for level "<<level<<" = "<<nobl);
               UBLOG(logINFO, "Number of nodes for level "<<level<<" = "<<nobl*numberOfNodesPerBlock);
            }
            UBLOG(logINFO, "Necessary memory  = "<<needMemAll<<" bytes");
            UBLOG(logINFO, "Necessary memory per process = "<<needMem<<" bytes");
            UBLOG(logINFO, "Available memory per process = "<<availMem<<" bytes");
         }

         //LBMKernelPtr kernel = LBMKernelPtr(new CompressibleCumulantLBMKernel(blockNx[0], blockNx[1], blockNx[2], CompressibleCumulantLBMKernel::NORMAL));
         ////LBMKernelPtr kernel = LBMKernelPtr(new IncompressibleCumulantLBMKernel(blockNx[0], blockNx[1], blockNx[2], IncompressibleCumulantLBMKernel::NORMAL));

         //BCProcessorPtr bcProc;

         //if (thinWall)
         //{
            //bcProc = BCProcessorPtr(new ThinWallBCProcessor());
         //}
         //else
         //{
            //bcProc = BCProcessorPtr(new BCProcessor());
         //}

         //kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (myid==0) UBLOG(logINFO, "SetKernelBlockVisitor:end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 5");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         if (refineLevel>0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         if (myid==0) UBLOG(logINFO, "SetUndefinedNodesBlockVisitor:end");

         //BC
         intHelper.setBC();
         if (myid==0) UBLOG(logINFO, "intHelper.setBC():end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 6");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }



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

         InitDistributionsBlockVisitor initVisitor1(nuLB, rhoLB);
         initVisitor1.setVx1(fct);
         ////initVisitor.setVx1(inflowProfileVx1);
         ////initVisitor.setVx2(inflowProfileVx2);
         ////initVisitor.setVx3(inflowProfileVx3);
         ////initVisitor.setNu(nuLB);
         grid->accept(initVisitor1);



         ////set connectors
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         //InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         
         //Grid3DPtr oldGrid(new Grid3D(comm));
         //
         ////with MPIIORestartCoProcessor
         //UbSchedulerPtr iSch(new UbScheduler());
         //MPIIORestart1CoProcessor rcpInit(oldGrid, iSch, pathReInit, comm);
         //rcpInit.setChunk(chunk);
         //rcpInit.restart(stepReInit);

         //////with MPIIORestart2CoProcessor
         ////UbSchedulerPtr iSch(new UbScheduler());
         ////MPIIORestart2CoProcessor rcp(oldGrid, iSch, pathReInit, comm);
         ////rcp.readBlocks(stepReInit);
         ////Grid3DVisitorPtr newMetisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////oldGrid->accept(newMetisVisitor);
         ////rcp.readDataSet(stepReInit);
         ////rcp.readBoundaryConds(stepReInit);

         //InitDistributionsWithInterpolationGridVisitor initVisitor(oldGrid, iProcessor, nuLB);
         //grid->accept(initVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //bcVisitor should be accept after initialization!!!!
         grid->accept(bcVisitor);
         if (myid==0) UBLOG(logINFO, "grid->accept(bcVisitor):end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 7");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }


         //Postrozess
         {
            UbSchedulerPtr geoSch(new UbScheduler(1));
            WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm);
            ppgeo.process(0);
         }

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 8");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         //fngIntrWhole1.reset();
         fngIntrWhole2.reset();

         ////UbSchedulerPtr rSch(new UbScheduler(cpStep, cpStart));
         ////MPIIORestartCoProcessor rcp(grid, rSch, pathOut, comm);

         GbCuboid3DPtr sponfeLayerBB1(new GbCuboid3D(g_maxX1-750, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitor(sponfeLayerBB1);

         if (myid==0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         
         //rcp2.process(restartStep);
         //return;
         //////////////////////////////////////////////////////////////////////////
         //////MPIIORestart2CoProcessor
         //UbSchedulerPtr iSch(new UbScheduler());
         //rcp2.readBlocks(restartStep);
         //grid->updateDistributedBlocks(comm);
         
         //Grid3DVisitorPtr newMetisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         //grid->accept(newMetisVisitor);
         
         //rcp2.restart((int)restartStep);
         //grid->setTimeStep(restartStep);
         
         //rcp.readBlocks(restartStep);
         //rcp.readDataSet(restartStep);
         //rcp.readBoundaryConds(restartStep);
         //grid->setTimeStep(restartStep);
         
         //setBC(grid, pathGeo, fngFileWhole2, zigZagTape, boundingBox, uLB, rhoLB, blockLength, bcProc);
         
         //rp.process(restartStep);

         rcp3.restart((int)restartStep);
         grid->setTimeStep(restartStep);
         
         //{
            //WriteBlocksCoProcessor ppblocks(grid, UbSchedulerPtr(new UbScheduler(1)), pathOut+"/mpiio3", WbWriterVtkXmlASCII::getInstance(), comm);
            //ppblocks.process(0);
         //}
         
         //{
            //UbSchedulerPtr stepSch(new UbScheduler(1));
            //WriteMacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathOut+"/mpiio3", WbWriterVtkXmlBinary::getInstance(), conv, comm);
            //pp.process(restartStep);
         //} 
 
         //{
           
            //UbSchedulerPtr geoSch(new UbScheduler(1));
            //WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut+"/mpiio3", WbWriterVtkXmlBinary::getInstance(), conv, comm);
            //ppgeo.process(0);
         //}

         //rcp3.process(restartStep);
         
         //return;
         
         
     
            
         
                 
         ////////////////////////////////////////////////////////////////////////////
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         grid->accept(bcVisitor);

         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);
         
      }

      UbSchedulerPtr nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      UbSchedulerPtr stepSch(new UbScheduler(outTimeStep,outTimeStart));



      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 9");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }

      ////////////////////////////////////////////////////////////////////////
      ////MPIIORestart2CoProcessor 
      //grid->deleteBlockIDs();
      //RenumberBlockVisitor renumber;
      //grid->accept(renumber);
      //UbSchedulerPtr iSch(new UbScheduler(1));
      //MPIIORestart2CoProcessor rcpInit(grid, iSch, pathOut+"/mpiio2", comm);
      //rcpInit.process(0);
      ////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      ////MPIIORestartCoProcessor 
      //UbSchedulerPtr iSch(new UbScheduler(1));
      //MPIIORestartCoProcessor rcpInit(grid, iSch, pathOut, comm);
      //rcpInit.process(0);
      //////////////////////////////////////////////////////////////////////////

      WriteMacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm);
      //pp.process(0);

      //rcp.process(0);

      //return;

      //////////////////////////////////////////////////////////////////////////
      ////Forces calculation
      //////////////////////////////////////////////////////////////////////////
      //if (myid==0) UBLOG(logINFO, "Read fngFileWhole2:start");
      //GbTriFaceMesh3DPtr fngMeshWhole2 = GbTriFaceMesh3DPtr(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo+"/"+fngFileWhole2, "fngMeshWhole2", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
      //if (myid==0) UBLOG(logINFO, "Read fngFileWhole2:end");
      //fngMeshWhole2->rotate(0.0, 0.5, 0.0);
      //D3Q27TriFaceMeshInteractorPtr fngIntrWhole2 = D3Q27TriFaceMeshInteractorPtr(new D3Q27TriFaceMeshInteractor(fngMeshWhole2, grid, noSlipBCAdapter, Interactor3D::SOLID));

      //SetSolidBlockVisitor fngVisitor(fngIntrWhole2, SetSolidBlockVisitor::BC);
      //grid->accept(fngVisitor);
      //fngIntrWhole2->initInteractor();
      
      //grid->accept(bcVisitor);

      //string ZckbndFilename = pathGeo+"/"+zigZagTape;
      //GbTriFaceMesh3DPtr meshBand1(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape1"));
      //meshBand1->rotate(0.0, 5, 0.0);
      //meshBand1->translate(15, 0, -12.850);
      //// Zackenband2
      //GbTriFaceMesh3DPtr meshBand2(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape2"));
      //meshBand2->rotate(0.0, 5, 0.0);
      //meshBand2->translate(15, 5, -12.850);

      //GbTriFaceMesh3DPtr meshBand5(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape5"));
      //meshBand5->rotate(0.0, -1, 0.0);
      //meshBand5->rotate(0.0, 0.0, 180.0);
      //meshBand5->translate(30, 0, -37.2);
      //// Zackenband6
      //GbTriFaceMesh3DPtr meshBand6(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(ZckbndFilename, "zigZagTape6"));
      //meshBand6->rotate(0.0, -1, 0.0);
      //meshBand6->rotate(0.0, 0.0, 180.0);
      //meshBand6->translate(30, 5, -37.2);

      //D3Q27TriFaceMeshInteractorPtr triBand1Interactor(new D3Q27TriFaceMeshInteractor(meshBand1, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
      //D3Q27TriFaceMeshInteractorPtr triBand2Interactor(new D3Q27TriFaceMeshInteractor(meshBand2, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
      //D3Q27TriFaceMeshInteractorPtr triBand3Interactor(new D3Q27TriFaceMeshInteractor(meshBand5, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));
      //D3Q27TriFaceMeshInteractorPtr triBand4Interactor(new D3Q27TriFaceMeshInteractor(meshBand6, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES));

      //SetSolidOrTransBlockVisitor band1Visitor(triBand1Interactor, SetSolidOrTransBlockVisitor::TRANS);
      //grid->accept(band1Visitor);
      //triBand1Interactor->initInteractor();

      //SetSolidOrTransBlockVisitor band2Visitor(triBand2Interactor, SetSolidOrTransBlockVisitor::TRANS);
      //grid->accept(band2Visitor);
      //triBand2Interactor->initInteractor();

      //SetSolidOrTransBlockVisitor band3Visitor(triBand3Interactor, SetSolidOrTransBlockVisitor::TRANS);
      //grid->accept(band3Visitor);
      //triBand3Interactor->initInteractor();

      //SetSolidOrTransBlockVisitor band4Visitor(triBand4Interactor, SetSolidOrTransBlockVisitor::TRANS);
      //grid->accept(band4Visitor);
      //triBand4Interactor->initInteractor();

      //double b    = 30; //wingspan
      //double t    = 300; //chord length
      //double area = (b*t)/(deltaXcoarse*deltaXcoarse);
      //double v    = uLB;

      //CalculateForcesCoProcessor fp(grid, stepSch, pathOut + "/forces/forces.txt", comm, v, area);
      //fp.addInteractor(fngIntrWhole2);
      //fp.addInteractor(triBand1Interactor);
      //fp.addInteractor(triBand2Interactor);
      //fp.addInteractor(triBand3Interactor);
      //fp.addInteractor(triBand4Interactor);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      ////Cp calculation
      //////////////////////////////////////////////////////////////////////////
      //UbSchedulerPtr pcpSch(new UbScheduler(1, pcpStart, pcpStop));

      //double planeCenter = g_minX2+(g_maxX2-g_minX2)/2.0;
      //double planeX2min = planeCenter-deltaXfine;
      //double planeX2max = planeCenter;//planeCenter+deltaXfine;
      //GbCuboid3DPtr plane(new GbCuboid3D(g_minX1,planeX2min,g_minX3,g_maxX1,planeX2max,g_maxX3));
      //if (myid==0) GbSystem3D::writeGeoObject(plane.get(), pathOut+"/geo/plane", WbWriterVtkXmlASCII::getInstance());

      //PressureCoefficientCoProcessor pcp(grid, pcpSch, plane, pathOut+"/cp/cp", comm);
      //pcp.addInteractor(fngIntrWhole2);
      //////////////////////////////////////////////////////////////////////////

      UbSchedulerPtr tavSch(new UbScheduler(1, timeAvStart, timeAvStop));
      TimeAveragedValuesCoProcessorPtr tav(new TimeAveragedValuesCoProcessor(grid, pathOut, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,
        TimeAveragedValuesCoProcessor::Density | TimeAveragedValuesCoProcessor::Velocity | TimeAveragedValuesCoProcessor::Fluctuations));
      tav->setWithGhostLayer(true);

      IntegrateValuesHelperPtr mic1(new IntegrateValuesHelper(grid, comm,300-deltaXcoarse,35,-600-deltaXcoarse,
         300,65,-600));
      if (myid==0) GbSystem3D::writeGeoObject(mic1->getBoundingBox().get(), pathOut+"/geo/mic1", WbWriterVtkXmlBinary::getInstance());
      UbSchedulerPtr stepMV(new UbScheduler(1));
      //TimeseriesCoProcessor tsp1(grid, stepMV, mic1, pathOut+"/mic/mic1", comm);

      //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
      //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, tavSch, CalculationManager::MPI));
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, tavSch, CalculationManager::PrePostBc));

      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
      if (myid==0) UBLOG(logINFO, "Simulation-end");
      
      ////////////////////////////////////////////////////////////////////////
      //MPIIORestart2CoProcessor 
      //grid->deleteBlockIDs();
      //RenumberBlockVisitor renumber;
      //grid->accept(renumber);
      //UbSchedulerPtr iSch(new UbScheduler(1));
      //MPIIORestart2CoProcessor rcpInit(grid, iSch, pathOut+"/mpiio2", comm);
      //rcpInit.process(0);
      ////////////////////////////////////////////////////////////////////////

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 10");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }
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

int main(int argc, char* argv[])
{

   if (argv!=NULL)
   {
      if (argv[1]!=NULL)
      {
         run(string(argv[1]));
      }
      else
      {
         cout<<"Configuration file must be set!: "<<argv[0]<<" <config file>"<<endl<<std::flush;
      }
   }

   return 0;
}


#include <iostream>
#include <string>

#include <PointerDefinitions.h>
#include "VirtualFluids.h"
#include <omp.h>
using namespace std;


//////////////////////////////////////////////////////////////////////////
void reinitGrid(SPtr<Grid3D> grid)
{
   int minInitLevel = grid->getCoarsestInitializedLevel();
   int maxInitLevel = grid->getFinestInitializedLevel();
   int gridRank = grid->getRank();

   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      std::vector<SPtr<Block3D>> blockVector;
      grid->getBlocks(level, gridRank, blockVector);
      for (SPtr<Block3D> block : blockVector)
      {
         if (block)
         {
            SPtr<ILBMKernel> kernel = block->getKernel();
            if (!kernel)
               throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: "+block->toString());

            SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

            LBMReal f[D3Q27System::ENDF+1];
            LBMReal vx1,vx2,vx3,rho;

            size_t nx1 = distributions->getNX1();
            size_t nx2 = distributions->getNX2();
            size_t nx3 = distributions->getNX3();

            for (int ix3=0; ix3<bcArray->getNX3(); ix3++)
               for (int ix2=0; ix2<bcArray->getNX2(); ix2++)
                  for (int ix1=0; ix1<bcArray->getNX1(); ix1++)
                  {
                     distributions->getDistribution(f, ix1, ix2, ix3);
                     D3Q27System::calcCompMacroscopicValues(f,rho,vx1,vx2,vx3);
                     if (rho > 0.1 || rho < -0.1)
                     {
                        D3Q27System::calcCompFeq(f, 0, 0, 0, 0);
                        distributions->setDistribution(f, ix1, ix2, ix3);
                        distributions->setDistributionInv(f, ix1, ix2, ix3);
                     }
                  }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void voxelMatrixTransformation(vector<int> pmNX, double lthreshold, double uthreshold, vector<float> voxelDeltaX, string pathGeoTEvoxel, string pathOut, SPtr<Communicator> comm)
{  
   int myid = comm->getProcessID();

   vector<double> origin(3);           
   origin[0] = 0;                     
   origin[1] = 0;
   origin[2] = 0;

   SPtr<GbVoxelMatrix3D> voxelMatrix (new GbVoxelMatrix3D(pmNX[0], pmNX[1], pmNX[2], 0, lthreshold, uthreshold));
   voxelMatrix->readMatrixFromRawFile<unsigned short>(pathGeoTEvoxel, GbVoxelMatrix3D::BigEndian);
   voxelMatrix->setVoxelMatrixDelta(voxelDeltaX[0], voxelDeltaX[1], voxelDeltaX[2]);
   voxelMatrix->setVoxelMatrixMininum(origin[0], origin[1], origin[2]);

   voxelMatrix->rotate90aroundZ();
   voxelMatrix->rotate90aroundZ();
   voxelMatrix->rotate90aroundZ();
   voxelMatrix->rotate90aroundX();
   voxelMatrix->translate(0.2813, 0, 0.0042);
   //voxelMatrix->translate(0.3, 0, 0.006);

   if (myid==0) voxelMatrix->writeToVTKImageDataAppended(pathOut+"/geo/fngTEvoxel");

}
//////////////////////////////////////////////////////////////////////////
void initPteBlock(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   int gridRank = grid->getRank();
   int blockRank = block->getRank();

   if (blockRank == gridRank)
   {
      SPtr<ILBMKernel> kernel = block->getKernel();
      if (!kernel)
         throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: "+block->toString());

      SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();
      SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

      LBMReal f[D3Q27System::ENDF+1];

      size_t nx1 = distributions->getNX1();
      size_t nx2 = distributions->getNX2();
      size_t nx3 = distributions->getNX3();

      for (int ix3=0; ix3<bcArray->getNX3(); ix3++)
         for (int ix2=0; ix2<bcArray->getNX2(); ix2++)
            for (int ix1=0; ix1<bcArray->getNX1(); ix1++)
            {
               D3Q27System::calcCompFeq(f, 0, 0, 0, 0);
               distributions->setDistribution(f, ix1, ix2, ix3);
               distributions->setDistributionInv(f, ix1, ix2, ix3);
            }
      block->setActive(true);
   }
}
//////////////////////////////////////////////////////////////////////////
void initPteFs(SPtr<Grid3D> grid, vector<SPtr<Block3D>>& vectorTE)
{
   for (SPtr<Block3D> block : vectorTE)
   {
      initPteBlock(grid, block);
   }
}
//////////////////////////////////////////////////////////////////////////
void initPteBC(SPtr<Grid3D> grid, vector<SPtr<Block3D>>& vectorTE, string fngFileTEmesh, int parts, string pathGeoTE, string pathOut, SPtr<Communicator> comm)
{
   int myid = comm->getProcessID();

   if (myid==0) UBLOG(logINFO, "vectorTE.size="<<vectorTE.size());

   for (int i = 0; i < parts; i++)
   {
      SPtr<GbTriFaceMesh3D> fngMeshTEmesh;
      if (myid==0) UBLOG(logINFO, "Read fngMeshTEmesh "<<fngFileTEmesh<<"_"<<UbSystem::toString(i)<<".stl"<<":start");
      fngMeshTEmesh = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeoTE+"/"+fngFileTEmesh+"_"+UbSystem::toString(i)+".stl",
         "fngMeshTEmesh", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
      if (myid==0) UBLOG(logINFO, "Read fngMeshTEmesh:end");
      //if (myid==0) GbSystem3D::writeGeoObject(fngMeshTEmesh.get(), pathOut+"/geo/fngMeshTEmesh", WbWriterVtkXmlASCII::getInstance());

      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
      SPtr<Interactor3D> fngIntrTEmesh = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshTEmesh, grid, noSlipBCAdapter, Interactor3D::SOLID));
      //fngMeshTEmesh.reset();

      SPtr<GbCuboid3D> boundingBox(new GbCuboid3D(fngMeshTEmesh->getX1Minimum(), fngMeshTEmesh->getX2Minimum(), fngMeshTEmesh->getX3Minimum(),
         fngMeshTEmesh->getX1Maximum(), fngMeshTEmesh->getX2Maximum(), fngMeshTEmesh->getX3Maximum()));

      for (SPtr<Block3D> block : vectorTE)
      {
         UbTupleDouble3 coords = grid->getBlockWorldCoordinates(block);
         UbTupleDouble3 deltas = grid->getBlockLengths(block);
         if (boundingBox->isCellInsideOrCuttingGbObject3D(val<1>(coords)
            , val<2>(coords)
            , val<3>(coords)
            , val<1>(coords)+val<1>(deltas)
            , val<2>(coords)+val<2>(deltas)
            , val<3>(coords)+val<3>(deltas)))
         {
            block->setActive(true);
            fngIntrTEmesh->setDifferencesToGbObject3D(block);
         }
      }

      //for (SPtr<Block3D> block : bb)
      //{
      //   block->setActive(true);
      //   fngIntrTEmesh->setDifferencesToGbObject3D(block);
      //   initPteBlock(grid, block);
      //}

      dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrTEmesh)->setQs(0);
   }

   //fngIntrTE->removeSolidBlocks();
   //fngIntrTE->removeBcBlocks();
}
//////////////////////////////////////////////////////////////////////////
void initPte(SPtr<Grid3D> grid, SPtr<Interactor3D> fngIntrTE, SPtr<Interactor3D> fngIntrTEmesh, string pathOut, SPtr<Communicator> comm)
{
   SetSolidOrBoundaryBlockVisitor v1(fngIntrTE, SetSolidOrBoundaryBlockVisitor::SOLID);
   grid->accept(v1);
   SetSolidOrBoundaryBlockVisitor v2(fngIntrTE, SetSolidOrBoundaryBlockVisitor::BC);
   grid->accept(v2);
   std::vector<SPtr<Block3D>>& sb = fngIntrTE->getSolidBlockSet();
   std::vector<SPtr<Block3D>>& bb = fngIntrTE->getBcBlocks();

   for (SPtr<Block3D> block : sb)
   {
      block->setActive(true);
      fngIntrTEmesh->setDifferencesToGbObject3D(block);
      initPteBlock(grid, block);
   }

   fngIntrTE->removeSolidBlocks();

   for (SPtr<Block3D> block : bb)
   {
      block->setActive(true);
      fngIntrTEmesh->setDifferencesToGbObject3D(block);
      initPteBlock(grid, block);;
   }

   dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrTEmesh)->setQs(0);


   fngIntrTE->removeBcBlocks();
}
//////////////////////////////////////////////////////////////////////////
void setPointsTE(vector<SPtr<Block3D>>& vectorTE)
{
   for (SPtr<Block3D> block : vectorTE)
   {
      block->setActive(true);

      SPtr<ILBMKernel> kernel = block->getKernel();
      if (!kernel)
         throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: "+block->toString());

      SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();
      SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

      LBMReal f[D3Q27System::ENDF+1];

      size_t nx1 = distributions->getNX1();
      size_t nx2 = distributions->getNX2();
      size_t nx3 = distributions->getNX3();

      int minX1 = 0;
      int minX2 = 0;
      int minX3 = 0;

      int maxX1 = (int)bcArray->getNX1();
      int maxX2 = (int)bcArray->getNX2();
      int maxX3 = (int)bcArray->getNX3();

      for (int ix3=minX3; ix3<maxX3; ix3++)
         for (int ix2=minX2; ix2<maxX2; ix2++)
            for (int ix1=minX1; ix1<maxX1; ix1++)
            {
               if (bcArray->isSolid(ix1, ix2, ix3) || bcArray->isUndefined(ix1,ix2,ix3))
               {
                  bool fFlag = true;
                  for (int k=-1; k<=1; k++)
                     for (int j=-1; j<=1; j++)
                        for (int i=-1; i<=1; i++)
                        {
                           if (ix1+i>minX1 && ix1+i<maxX1 && ix2+j>minX2 && ix2+j<maxX2 && ix3+k>minX3 && ix3+k<maxX3)
                           {
                              if (bcArray->isFluid(ix1+i, ix2+j, ix3+k))
                              {
                                 fFlag = fFlag && true;
                              }
                              else
                              {
                                 fFlag = fFlag && false;
                              }
                           }
                        }
                  if (fFlag)
                  {
                     bcArray->setFluid(ix1, ix2, ix3);
                     D3Q27System::calcCompFeq(f, 0, 0, 0, 0);
                     distributions->setDistribution(f, ix1, ix2, ix3);
                     distributions->setDistributionInv(f, ix1, ix2, ix3);
                  }
               }
               /*else if (bcArray->isFluid(ix1, ix2, ix3))
               {
                  bool sFlag = false;
                  for (int k=-1; k<=1; k++)
                     for (int j=-1; j<=1; j++)
                        for (int i=-1; i<=1; i++)
                        {
                           if (ix1+i>minX1 && ix1+i<maxX1 && ix2+j>minX2 && ix2+j<maxX2 && ix3+k>minX3 && ix3+k<maxX3)
                           {
                              if (bcArray->isSolid(ix1+i, ix2+j, ix3+k))
                              {
                                 sFlag = sFlag && true;
                              }
                              else
                              {
                                 sFlag = sFlag && false;
                              }
                           }
                        }
                  if (sFlag)
                  {
                     bcArray->setSolid(ix1, ix2, ix3);
                  }
               }*/
            }
   }
}


void run(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathOut = config.getValue<string>("pathOut");
      string          pathGeo = config.getValue<string>("pathGeo");
      //string          pathGeoTE = config.getValue<string>("pathGeoTE");
      string          fngFileNoTapeFull = config.getValue<string>("fngFileNoTapeFull");
      string          fngFileFull = config.getValue<string>("fngFileFull");
      string          fngFileNoTapeBody = config.getValue<string>("fngFileNoTapeBody");
      string          fngFileBody = config.getValue<string>("fngFileBody");
      string          fngFileTEmesh = config.getValue<string>("fngFileTEmesh");
      string          fngFileTE = config.getValue<string>("fngFileTE");
      //string          tapeFile = config.getValue<string>("tapeFile");
      int             accuracy = config.getValue<int>("accuracy");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blockNx = config.getVector<int>("blockNx");
      vector<double>  boundingBox = config.getVector<double>("boundingBox");
      double          restartStep = config.getValue<double>("restartStep");
      double          cpStart = config.getValue<double>("cpStart");
      double          cpStep = config.getValue<double>("cpStep");
      int             endTime = config.getValue<int>("endTime");
      double          outTimeStep = config.getValue<double>("outTimeStep");
      double          outTimeStart = config.getValue<double>("outTimeStart");
      double          availMem = config.getValue<double>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      bool            logToFile = config.getValue<bool>("logToFile");
      double          deltaXfine = config.getValue<double>("deltaXfine");
      double          refineDistance = config.getValue<double>("refineDistance");
      double          startDistance = config.getValue<double>("startDistance");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");
      bool            newStart = config.getValue<bool>("newStart");
      bool            writeBlocks = config.getValue<bool>("writeBlocks");
      string          pathReInit = config.getValue<string>("pathReInit");
      int             stepReInit = config.getValue<int>("stepReInit");
      bool            reinit = config.getValue<bool>("reinit");

      double          pcpStart = config.getValue<double>("pcpStart");
      double          pcpStop  = config.getValue<double>("pcpStop");

      double          timeAvStart       = config.getValue<double>("timeAvStart");
      double          timeAvStop        = config.getValue<double>("timeAvStop");

      vector<int>     pmNX              = config.getVector<int>("pmNX");
      double          lthreshold        = config.getValue<double>("lthreshold");
      double          uthreshold        = config.getValue<double>("uthreshold");
      vector<float>   voxelDeltaX       = config.getVector<float>("voxelDeltaX");
      string          pathGeoTEvoxel    = config.getValue<string>("pathGeoTEvoxel");
      


      SPtr<Communicator> comm = MPICommunicator::getInstance();
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

      double g_minX1 = boundingBox[0];//*1000.0;
      double g_minX2 = boundingBox[2];//*1000.0;
      double g_minX3 = boundingBox[4];//*1000.0;
      double g_maxX1 = boundingBox[1];//*1000.0;
      double g_maxX2 = boundingBox[3];//*1000.0;
      double g_maxX3 = boundingBox[5];//*1000.0;
      //deltaXfine *=1000.0;

      //////////////////////////////////////////////////////////////////////////
      double deltaXcoarse = deltaXfine*(double)(1<<refineLevel);
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
      double csReal = uReal / Ma;
      double hLB = lReal / deltaXcoarse;

      LBMUnitConverter unitConverter(lReal, csReal, rhoReal, hLB);

      double uLB = uReal   * unitConverter.getFactorVelocityWToLb();
      double nuLB = nuReal * unitConverter.getFactorViscosityWToLb();
      double lLB = lReal/deltaXcoarse;
      //double nuLB = (uLB*lLB)/Re; //0.005;
      //double nuLB = 0.005;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      ////////////////////////////////////////////////////////////////////////
      //Grid
      //////////////////////////////////////////////////////////////////////////
      SPtr<Grid3D> grid(new Grid3D(comm));

      //BC adapters
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new ThinWallNoSlipBCAlgorithm()));

      //SPtr<BCAdapter> slipBCAdapter(new SlipBCAdapter());
      //slipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new SlipBCAlgorithm()));

      mu::Parser fct;
      fct.SetExpr("U");
      fct.DefineConst("U", uLB);
      SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

      //fct.SetExpr("U");
      //fct.DefineConst("U", 0.01);
      //SPtr<BCAdapter> velBCAdapterOut(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
      //velBCAdapterOut->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));

      SPtr<BCAdapter> outflowBCAdapter(new DensityBCAdapter(rhoLB));
      outflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(velBCAdapter);
      bcVisitor.addBC(outflowBCAdapter);

      SPtr<BCProcessor> bcProc;
      //bcProc = SPtr<BCProcessor>(new BCProcessor());
      bcProc = SPtr<BCProcessor>(new ThinWallBCProcessor());

      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
      //dynamicPointerCast<CompressibleCumulantLBMKernel>(kernel)->setRelaxationParameter(CompressibleCumulantLBMKernel::NORMAL);

      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      //dynamicPointerCast<CompressibleCumulant4thOrderViscosityLBMKernel>(kernel)->setBulkViscosity(nuLB*2.0e3);

      kernel->setBCProcessor(bcProc);

      SPtr<LBMKernel> spKernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
      spKernel->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIORestartCoProcessor> restartCoProcessor(new MPIIORestartCoProcessor(grid, rSch, pathOut, comm));
      //SPtr<MPIIOMigrationCoProcessor> restartCoProcessor(new MPIIOMigrationCoProcessor(grid, rSch, pathOut, comm));
      restartCoProcessor->setLBMKernel(kernel);
      restartCoProcessor->setBCProcessor(bcProc);

      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIOMigrationCoProcessor> migCoProcessor(new MPIIOMigrationCoProcessor(grid, mSch, pathOut+"/mig", comm));
      migCoProcessor->setLBMKernel(kernel);
      migCoProcessor->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 2");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }


      if (newStart)
      {
         ////////////////////////////////////////////////////////////////////////
         //define grid
         //////////////////////////////////////////////////////////////////////////
         grid->setDeltaX(deltaXcoarse);
         grid->setBlockNX(blockNx[0], blockNx[1], blockNx[2]);

         SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
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
            UBLOG(logINFO, "* dx_base             = "<<deltaXcoarse<<" m");
            UBLOG(logINFO, "* dx_refine           = "<<deltaXfine<<" m");
            UBLOG(logINFO, "* blocknx             = "<<blockNx[0]<<"x"<<blockNx[1]<<"x"<<blockNx[2]);
            UBLOG(logINFO, "* refineDistance      = "<<refineDistance);
            UBLOG(logINFO, "* number of levels    = "<<refineLevel+1);
            UBLOG(logINFO, "* number of threads   = "<<numOfThreads);
            UBLOG(logINFO, "* number of processes = "<<comm->getNumberOfProcesses());
            UBLOG(logINFO, "* path = "<<pathOut);
            UBLOG(logINFO, "Preprozess - start");
         }

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 3");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }
         //////////////////////////////////////////////////////////////////////////


         //voxelMatrixTransformation(pmNX, lthreshold, uthreshold, voxelDeltaX, pathGeoTEvoxel, pathOut, comm);

         //return;


         SPtr<GbTriFaceMesh3D> fngMeshTE;
         if (myid==0) UBLOG(logINFO, "Read fngMeshTE:start");
         fngMeshTE = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileTE, "fngMeshTE", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         if (myid==0) UBLOG(logINFO, "Read fngMeshTE:end");
         fngMeshTE->rotate(0.0, 0.5, 0.0);
         fngMeshTE->translate(0.0, 0.0, 0.0012);
         if (myid==0) GbSystem3D::writeGeoObject(fngMeshTE.get(), pathOut+"/geo/fngMeshTE", WbWriterVtkXmlBinary::getInstance());

         SPtr<Interactor3D> fngIntrTE = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshTE, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));



         if (refineLevel>0 && myid==0 && writeBlocks)
         {
            if (myid==0) UBLOG(logINFO, "Refinement - start");
            int rank = grid->getRank();
            grid->setRank(0);

            SPtr<GbTriFaceMesh3D> fngMeshNoTapeFull;
            if (myid==0) UBLOG(logINFO, "Read fngFileNoTapeFull:start");
            fngMeshNoTapeFull = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileNoTapeFull, "fngMeshNoTapeBody", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read fngFileNoTapeFull:end");
            fngMeshNoTapeFull->rotate(0.0, 0.5, 0.0);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshNoTapeFull.get(), pathOut+"/geo/fngMeshNoTapeFull", WbWriterVtkXmlBinary::getInstance());

            SPtr<Interactor3D> fngIntrNoTapeFull = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshNoTapeFull, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));

            int level;

            //level 1
            level = 1;
            if (refineLevel - level >= 0)
            {
               SPtr<GbObject3D> refBoxL1(new GbCuboid3D(0.12, -0.02625, -0.03, 1.0, 0.06, 0.03));
               if (myid==0) GbSystem3D::writeGeoObject(refBoxL1.get(), pathOut+"/geo/refBoxL1", WbWriterVtkXmlASCII::getInstance());
               RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL1(refBoxL1, level);
               grid->accept(refVisitorBoxL1);

               dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrNoTapeFull)->refineBlockGridToLevel(level, startDistance, 24.0*refineDistance);
            }

            //level 2
            level = 2;
            if (refineLevel - level >= 0)
            {
               SPtr<GbObject3D> refBoxL2(new GbCuboid3D(0.15, -0.03, -0.015, 0.7, 0.06, 0.015));
               if (myid==0) GbSystem3D::writeGeoObject(refBoxL2.get(), pathOut+"/geo/refBoxL2", WbWriterVtkXmlASCII::getInstance());
               RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL2(refBoxL2, level);
               grid->accept(refVisitorBoxL2);

               dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrNoTapeFull)->refineBlockGridToLevel(level, startDistance, 12.0*refineDistance);
            }

            //level 3
            level = 3;
            if (refineLevel - level >= 0)
            {
               SPtr<GbObject3D> refBoxL3(new GbCuboid3D(0.15, -0.03, -0.010, 0.4, 0.06, 0.012));
               if (myid==0) GbSystem3D::writeGeoObject(refBoxL3.get(), pathOut+"/geo/refBoxL3", WbWriterVtkXmlASCII::getInstance());
               RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL3(refBoxL3, level);
               grid->accept(refVisitorBoxL3);

               dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrNoTapeFull)->refineBlockGridToLevel(level, startDistance, 6.0*refineDistance);
            }

            //level 4
            level = 4;
            if (refineLevel - level >= 0)
            {
               SPtr<GbObject3D> refBoxL4(new GbCuboid3D(0.15, -0.03, -0.005, 0.32, 0.06, 0.01));
               if (myid==0) GbSystem3D::writeGeoObject(refBoxL4.get(), pathOut+"/geo/refBoxL4", WbWriterVtkXmlASCII::getInstance());
               RefineCrossAndInsideGbObjectBlockVisitor refVisitorBoxL4(refBoxL4, level);
               grid->accept(refVisitorBoxL4);

               dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrNoTapeFull)->refineBlockGridToLevel(level, startDistance, 3.0*refineDistance);
            }

            //level 5
            level = 5;
            if (refineLevel - level >= 0)
            {
               dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrNoTapeFull)->refineBlockGridToLevel(level, startDistance, refineDistance);


               //SPtr<GbTriFaceMesh3D> fngMeshTE;
               //if (myid==0) UBLOG(logINFO, "Read fngMeshTE:start");
               //fngMeshTE = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileTE, "fngMeshTE", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
               //if (myid==0) UBLOG(logINFO, "Read fngMeshTE:end");
               //fngMeshTE->rotate(0.0, 0.5, 0.0);
               //fngMeshTE->translate(0.0, 0.0, 0.0012);
               //if (myid==0) GbSystem3D::writeGeoObject(fngMeshTE.get(), pathOut+"/geo/fngMeshTE", WbWriterVtkXmlBinary::getInstance());

               ////SPtr<Interactor3D> fngIntrTE = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshTE, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));
               ////dynamicPointerCast<D3Q27TriFaceMeshInteractor>(fngIntrTE)->refineBlockGridToLevel(level, startDistance, refineDistance);
               RefineCrossAndInsideGbObjectBlockVisitor refVisitorTE(fngMeshTE, level);
               grid->accept(refVisitorTE);
            }

            ///////delete solid blocks
            if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");

            SPtr<GbTriFaceMesh3D> fngMeshNoTapeBody;
            if (myid==0) UBLOG(logINFO, "Read fngFileNoTapeBody:start");
            fngMeshNoTapeBody = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileNoTapeBody, "fngMeshNoTapeBody", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            if (myid==0) UBLOG(logINFO, "Read fngFileNoTapeBody:end");
            fngMeshNoTapeBody->rotate(0.0, 0.5, 0.0);
            fngMeshNoTapeBody->translate(0.0, 0.0, -0.00011);
            if (myid==0) GbSystem3D::writeGeoObject(fngMeshNoTapeBody.get(), pathOut+"/geo/fngMeshNoTapeBody", WbWriterVtkXmlBinary::getInstance());

            SPtr<Interactor3D> fngIntrNoTapeBody = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshNoTapeBody, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));//, Interactor3D::POINTS));

            SetSolidOrBoundaryBlockVisitor v(fngIntrNoTapeBody, SetSolidOrBoundaryBlockVisitor::SOLID);
            grid->accept(v);
            std::vector<SPtr<Block3D>>& sb = fngIntrNoTapeBody->getSolidBlockSet();
            for (SPtr<Block3D> block : sb)
            {
               grid->deleteBlock(block);
            }
            fngIntrNoTapeBody->removeSolidBlocks();
            fngIntrNoTapeBody->removeBcBlocks();


            if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - end");
            //////////////////////////////////////////

            {
               WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(0);
            }

            grid->setRank(rank);

            RatioBlockVisitor ratioVisitor(refineLevel);
            CheckRatioBlockVisitor checkRatio(refineLevel);
            int count = 0;

            do {
               if (myid==0) UBLOG(logINFO, "ratioVisitor - start");
               grid->accept(ratioVisitor);
               if (myid==0) UBLOG(logINFO, "ratioVisitor - end");
               if (myid==0) UBLOG(logINFO, "checkRatio - start");
               checkRatio.resetState();
               grid->accept(checkRatio);
               if (myid==0) UBLOG(logINFO, "checkRatio - end");
               if (myid==0) UBLOG(logINFO, "count = "<<count++<<" state = "<<checkRatio.getState());
            } while (!checkRatio.getState());


            {
               WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
               ppblocks.process(1);
            }

            OverlapBlockVisitor overlapVisitor(refineLevel, false);
            grid->accept(overlapVisitor);

            if (myid==0) UBLOG(logINFO, "Refinement - end");
         }
         else if (refineLevel>0 && !writeBlocks)
         {
            migCoProcessor->readBlocks(0);
         }
         grid->updateDistributedBlocks(comm);

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
         SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, velBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, velBCAdapter, Interactor3D::SOLID));

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathOut+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         //inflow
         SPtr<D3Q27Interactor> inflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //outflow
         SPtr<D3Q27Interactor> outflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowBCAdapter, Interactor3D::SOLID));

         //airfoil
         SPtr<GbTriFaceMesh3D> fngMeshBody;
         if (myid==0) UBLOG(logINFO, "Read fngFileBody:start");
         fngMeshBody = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileBody, "fngMeshBody", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         if (myid==0) UBLOG(logINFO, "Read fngFileBody:end");
         fngMeshBody->rotate(0.0, 0.5, 0.0);
         fngMeshBody->translate(0.0, 0.0, -0.00011);
         if (myid==0) GbSystem3D::writeGeoObject(fngMeshBody.get(), pathOut+"/geo/fngMeshBody", WbWriterVtkXmlBinary::getInstance());

         SPtr<Interactor3D> fngIntrBody = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshBody, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));
         fngMeshBody.reset();

         //SPtr<GbTriFaceMesh3D> fngMeshTEmesh;
         //if (myid==0) UBLOG(logINFO, "Read fngMeshTEmesh:start");
         //fngMeshTEmesh = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileTEmesh, "fngMeshTEmesh", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         //if (myid==0) UBLOG(logINFO, "Read fngMeshTEmesh:end");
         ////fngMeshTE->rotate(0.0, 0.5, 0.0);
         ////fngMeshTE->translate(0.0, 0.0, -0.0001);
         //if (myid==0) GbSystem3D::writeGeoObject(fngMeshTEmesh.get(), pathOut+"/geo/fngMeshTEmesh", WbWriterVtkXmlASCII::getInstance());


         //SPtr<Interactor3D> fngIntrTEmesh = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshTEmesh, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));
         //fngMeshTEmesh.reset();

         GbCuboid3DPtr geoAddWallP(new GbCuboid3D(0.269984375, g_minX2-blockLength, 0.0016, 0.2702421875, g_maxX2+blockLength, 0.0076));
         if (myid==0) GbSystem3D::writeGeoObject(geoAddWallP.get(), pathOut+"/geo/geoAddWallP", WbWriterVtkXmlASCII::getInstance());
         SPtr<D3Q27Interactor> addWallPIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoAddWallP, grid, noSlipBCAdapter, Interactor3D::SOLID));

         //////////////////////////////////////////////////////////////////////////
         vector<double> origin(3);
         origin[0] = 0;
         origin[1] = 0;
         origin[2] = 0;

         SPtr<GbVoxelMatrix3D> voxelMatrix1(new GbVoxelMatrix3D(pmNX[0], pmNX[1], pmNX[2], 0, lthreshold, uthreshold));
         voxelMatrix1->readMatrixFromRawFile<unsigned short>(pathGeoTEvoxel, GbVoxelMatrix3D::BigEndian);
         voxelMatrix1->setVoxelMatrixDelta(voxelDeltaX[0], voxelDeltaX[1], voxelDeltaX[2]);
         voxelMatrix1->setVoxelMatrixMininum(origin[0], origin[1], origin[2]);

         voxelMatrix1->rotate90aroundZ();
         voxelMatrix1->rotate90aroundZ();
         voxelMatrix1->rotate90aroundZ();
         voxelMatrix1->rotate90aroundX();
         voxelMatrix1->translate(0.2813, 0, 0.0042);
         double offset = ((g_maxX2-g_minX2)/2.0 - voxelMatrix1->getLengthX2())/2.0;
         voxelMatrix1->setVoxelMatrixMinX2(g_minX2+offset);

         if (myid==0) voxelMatrix1->writeToVTKImageDataAppended(pathOut+"/geo/fngTEvoxel1");

         SPtr<D3Q27Interactor> fngIntrTEvoxel1 = SPtr<D3Q27Interactor>(new D3Q27Interactor(voxelMatrix1, grid, noSlipBCAdapter, Interactor3D::SOLID));

         SPtr<GbVoxelMatrix3D> voxelMatrix2(new GbVoxelMatrix3D(pmNX[0], pmNX[1], pmNX[2], 0, lthreshold, uthreshold));
         voxelMatrix2->readMatrixFromRawFile<unsigned short>(pathGeoTEvoxel, GbVoxelMatrix3D::BigEndian);
         voxelMatrix2->setVoxelMatrixDelta(voxelDeltaX[0], voxelDeltaX[1], voxelDeltaX[2]);
         voxelMatrix2->setVoxelMatrixMininum(origin[0], origin[1], origin[2]);

         voxelMatrix2->rotate90aroundZ();
         voxelMatrix2->rotate90aroundZ();
         voxelMatrix2->rotate90aroundZ();
         voxelMatrix2->rotate90aroundX();
         voxelMatrix2->translate(0.2813, 0, 0.0042);
         voxelMatrix2->mirrorY();
         voxelMatrix2->setVoxelMatrixMinX2(voxelMatrix1->getX2Maximum());

         if (myid==0) voxelMatrix2->writeToVTKImageDataAppended(pathOut+"/geo/fngTEvoxel2");

         SPtr<D3Q27Interactor> fngIntrTEvoxel2 = SPtr<D3Q27Interactor>(new D3Q27Interactor(voxelMatrix2, grid, noSlipBCAdapter, Interactor3D::SOLID));
         //////////////////////////////////////////////////////////////////////////

         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(inflowIntr);
         intHelper.addInteractor(outflowIntr);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(fngIntrBody);
         intHelper.addInteractor(addWallPIntr);
         //intHelper.addInteractor(fngIntrTEmesh);
         intHelper.addInteractor(fngIntrTEvoxel1);
         intHelper.addInteractor(fngIntrTEvoxel2);
         intHelper.selectBlocks();

         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 4");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }
         //////////////////////////////////////

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

         if (writeBlocks)
         {
            migCoProcessor->writeBlocks(0);
         }

         {
            WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
            ppblocks.process(2);
         }

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
         if (myid==0) UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");

         if (myid==0) UBLOG(logINFO, "initPteBC:start");
         SetSolidOrBoundaryBlockVisitor v1(fngIntrTE, SetSolidOrBoundaryBlockVisitor::SOLID);
         grid->accept(v1);
         SetSolidOrBoundaryBlockVisitor v2(fngIntrTE, SetSolidOrBoundaryBlockVisitor::BC);
         grid->accept(v2);
         std::vector<SPtr<Block3D>>& vectorTE = fngIntrTE->getSolidBlockSet();
         std::vector<SPtr<Block3D>>& bb = fngIntrTE->getBcBlocks();
         vectorTE.insert(vectorTE.end(), bb.begin(), bb.end());
         //initPteBC(grid, vectorTE, fngFileTEmesh, 4, pathGeoTE, pathOut, comm);
         if (myid==0) UBLOG(logINFO, "initPteBC:end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 6");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         //initialization of distributions
         //InitDistributionsBlockVisitor initVisitor1;
         //initVisitor1.setVx1(fct);
         //grid->accept(initVisitor1);

         ////set connectors
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //initialization of distributions
         if (reinit)
         {
            InitDistributionsBlockVisitor initVisitor1;
            grid->accept(initVisitor1);
            SPtr<Grid3D> oldGrid(new Grid3D(comm));
            SPtr<UbScheduler> iSch(new UbScheduler());
            SPtr<MPIIORestartCoProcessor> rcp(new MPIIORestartCoProcessor(oldGrid, iSch, pathReInit, comm));
            rcp->setLBMKernel(kernel);
            rcp->setBCProcessor(bcProc);
            rcp->restart(stepReInit);
            InitDistributionsWithInterpolationGridVisitor initVisitor2(oldGrid, iProcessor, nuLB);
            grid->accept(initVisitor2);
         }
         else
         {
            InitDistributionsBlockVisitor initVisitor;
            initVisitor.setVx1(10.0/*fct*/);
            grid->accept(initVisitor);
         }

         initPteFs(grid, vectorTE);
         //initPte(grid, fngIntrTE, fngIntrTEmesh, pathOut, comm);

         //domain decomposition for threads
         //PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         //grid->accept(pqPartVisitor);

         //bcVisitor should be accept after initialization!!!!
         grid->accept(bcVisitor);
         if (myid==0) UBLOG(logINFO, "grid->accept(bcVisitor):end");

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 7");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         //Post process
         {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm);
            ppgeo.process(0);
         }

         if (myid==0)
         {
            UBLOG(logINFO, "PID = "<<myid<<" Point 8");
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
         }

         ////sponge layer
         ////////////////////////////////////////////////////////////////////////////

         GbCuboid3DPtr spongeLayerX1max(new GbCuboid3D(g_maxX1-0.4, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1max.get(), pathOut+"/geo/spongeLayerX1max", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX1max(spongeLayerX1max, spKernel, nuLB, D3Q27System::E);
         grid->accept(slVisitorX1max);

         GbCuboid3DPtr spongeLayerX1min(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1+0.2, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1min.get(), pathOut+"/geo/spongeLayerX1min", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX1min(spongeLayerX1min, spKernel, nuLB, D3Q27System::W);
         grid->accept(slVisitorX1min);

         GbCuboid3DPtr spongeLayerX3min(new GbCuboid3D(g_minX1+0.2, g_minX2-blockLength, g_minX3-blockLength, g_maxX1-0.4, g_maxX2+blockLength, g_minX3+0.2));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX3min.get(), pathOut+"/geo/spongeLayerX3min", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX3min(spongeLayerX3min, spKernel, nuLB, D3Q27System::B);
         grid->accept(slVisitorX3min);

         GbCuboid3DPtr spongeLayerX3max(new GbCuboid3D(g_minX1+0.2, g_minX2-blockLength, g_maxX3-0.2, g_maxX1-0.4, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX3max.get(), pathOut+"/geo/spongeLayerX3max", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX3max(spongeLayerX3max, spKernel, nuLB, D3Q27System::T);
         grid->accept(slVisitorX3max);

         /////////////////////////////////////////////////////////////////////////////
         if (myid==0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         //restartCoProcessor->restart((int)restartStep);
         migCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);

         WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
         ppblocks.process(3);
         ////////////////////////////////////////////////////////////////////////////
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //if (myid==0) UBLOG(logINFO, "setPointsTE:start");
         //SPtr<GbTriFaceMesh3D> fngMeshTE;
         //if (myid==0) UBLOG(logINFO, "Read fngMeshTE:start");
         //fngMeshTE = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileTE, "fngMeshTE", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         //if (myid==0) UBLOG(logINFO, "Read fngMeshTE:end");
         //fngMeshTE->rotate(0.0, 0.5, 0.0);
         //fngMeshTE->translate(0.0, 0.0, 0.0012);
         //if (myid==0) GbSystem3D::writeGeoObject(fngMeshTE.get(), pathOut+"/geo/fngMeshTE", WbWriterVtkXmlBinary::getInstance());
         //SPtr<Interactor3D> fngIntrTE = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshTE, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));
         //SetSolidOrBoundaryBlockVisitor v1(fngIntrTE, SetSolidOrBoundaryBlockVisitor::SOLID);
         //grid->accept(v1);
         //SetSolidOrBoundaryBlockVisitor v2(fngIntrTE, SetSolidOrBoundaryBlockVisitor::BC);
         //grid->accept(v2);
         //std::vector<SPtr<Block3D>>& vectorTE = fngIntrTE->getSolidBlockSet();
         //std::vector<SPtr<Block3D>>& bb = fngIntrTE->getBcBlocks();
         //vectorTE.insert(vectorTE.end(), bb.begin(), bb.end());
         //setPointsTE(vectorTE);
         //SPtr<UbScheduler> geoSch(new UbScheduler(1));
         //WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm);
         //ppgeo.process(1);
         //if (myid==0) UBLOG(logINFO, "setPointsTE:end");
//////////////////////////////////////////////////////////////////////////
         SPtr<GbTriFaceMesh3D> fngMeshBody;
         if (myid==0) UBLOG(logINFO, "Read fngFileBody:start");
         fngMeshBody = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo+"/"+fngFileBody, "fngMeshBody", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
         if (myid==0) UBLOG(logINFO, "Read fngFileBody:end");
         fngMeshBody->rotate(0.0, 0.5, 0.0);
         fngMeshBody->translate(0.0, 0.0, -0.00011);
         if (myid==0) GbSystem3D::writeGeoObject(fngMeshBody.get(), pathOut+"/geo/fngMeshBody2", WbWriterVtkXmlBinary::getInstance());
         SPtr<Interactor3D> fngIntrBody = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(fngMeshBody, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy));
         SetSolidOrBoundaryBlockVisitor v(fngIntrBody, SetSolidOrBoundaryBlockVisitor::BC);
         grid->accept(v);
         fngIntrBody->initInteractor();
//////////////////////////////////////////////////////////////////////////
         grid->accept(bcVisitor);

         ////sponge layer
         ////////////////////////////////////////////////////////////////////////////

         GbCuboid3DPtr spongeLayerX1max(new GbCuboid3D(g_maxX1-0.4, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1max.get(), pathOut+"/geo/spongeLayerX1max", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX1max(spongeLayerX1max, spKernel, nuLB, D3Q27System::E);
         grid->accept(slVisitorX1max);

         GbCuboid3DPtr spongeLayerX1min(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1+0.2, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX1min.get(), pathOut+"/geo/spongeLayerX1min", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX1min(spongeLayerX1min, spKernel, nuLB, D3Q27System::W);
         grid->accept(slVisitorX1min);

         GbCuboid3DPtr spongeLayerX3min(new GbCuboid3D(g_minX1+0.2, g_minX2-blockLength, g_minX3-blockLength, g_maxX1-0.4, g_maxX2+blockLength, g_minX3+0.2));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX3min.get(), pathOut+"/geo/spongeLayerX3min", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX3min(spongeLayerX3min, spKernel, nuLB, D3Q27System::B);
         grid->accept(slVisitorX3min);

         GbCuboid3DPtr spongeLayerX3max(new GbCuboid3D(g_minX1+0.2, g_minX2-blockLength, g_maxX3-0.2, g_maxX1-0.4, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(spongeLayerX3max.get(), pathOut+"/geo/spongeLayerX3max", WbWriterVtkXmlASCII::getInstance());
         SpongeLayerBlockVisitor slVisitorX3max(spongeLayerX3max, spKernel, nuLB, D3Q27System::T);
         grid->accept(slVisitorX3max);
      }

      SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      std::shared_ptr<NUPSCounterCoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> stepSch(new UbScheduler(outTimeStep, outTimeStart));

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Point 9");
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
      }

      SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      SPtr<GbObject3D> bbBox(new GbCuboid3D(g_minX1-blockLength, (g_maxX2-g_minX2)/2.0, g_minX3-blockLength, g_maxX1+blockLength, (g_maxX2-g_minX2)/2.0+deltaXcoarse, g_maxX3+blockLength));
      if (myid==0) GbSystem3D::writeGeoObject(bbBox.get(), pathOut+"/geo/bbBox", WbWriterVtkXmlASCII::getInstance());
      SPtr<WriteMQFromSelectionCoProcessor> writeMQSelectCoProcessor(new WriteMQFromSelectionCoProcessor(grid, stepSch, bbBox, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      SPtr<UbScheduler> tavSch(new UbScheduler(1, timeAvStart, timeAvStop));
      SPtr<TimeAveragedValuesCoProcessor> tav(new TimeAveragedValuesCoProcessor(grid, pathOut, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,
         TimeAveragedValuesCoProcessor::Density | TimeAveragedValuesCoProcessor::Velocity | TimeAveragedValuesCoProcessor::Fluctuations));
      tav->setWithGhostLayer(true);

      SPtr<IntegrateValuesHelper> mic1(new IntegrateValuesHelper(grid, comm, 0.3-deltaXfine, 0.015, 0.0005, 0.3, 0.015+deltaXfine, 0.0005+deltaXfine));
      if (myid==0) GbSystem3D::writeGeoObject(mic1->getBoundingBox().get(), pathOut+"/geo/mic1", WbWriterVtkXmlBinary::getInstance());
      SPtr<UbScheduler> stepMV(new UbScheduler(1, 0, 1000000));
      SPtr<TimeseriesCoProcessor> tsp1(new TimeseriesCoProcessor(grid, stepMV, mic1, pathOut+"/mic/mic1", comm));

      SPtr<IntegrateValuesHelper> mic2(new IntegrateValuesHelper(grid, comm, 0.3+deltaXfine, 0.015, 0.001685, 0.3, 0.015+deltaXfine, 0.001685+deltaXfine));
      if (myid==0) GbSystem3D::writeGeoObject(mic2->getBoundingBox().get(), pathOut+"/geo/mic2", WbWriterVtkXmlBinary::getInstance());
      SPtr<TimeseriesCoProcessor> tsp2(new TimeseriesCoProcessor(grid, stepMV, mic2, pathOut+"/mic/mic2", comm));

      SPtr<IntegrateValuesHelper> mic3(new IntegrateValuesHelper(grid, comm, 0.3-deltaXcoarse, 0.015, -0.46+4.25*deltaXcoarse, 0.3, 0.015+deltaXcoarse, -0.46+5.25*deltaXcoarse));
      if (myid==0) GbSystem3D::writeGeoObject(mic3->getBoundingBox().get(), pathOut+"/geo/mic3", WbWriterVtkXmlBinary::getInstance());
      SPtr<TimeseriesCoProcessor> tsp3(new TimeseriesCoProcessor(grid, stepMV, mic3, pathOut+"/mic/mic3", comm));

      SPtr<IntegrateValuesHelper> mic4(new IntegrateValuesHelper(grid, comm, 0.3-deltaXcoarse, 0.015, 0.46-5.25*deltaXcoarse, 0.3, 0.015+deltaXcoarse, 0.46-4.25*deltaXcoarse));
      if (myid==0) GbSystem3D::writeGeoObject(mic4->getBoundingBox().get(), pathOut+"/geo/mic4", WbWriterVtkXmlBinary::getInstance());
      SPtr<TimeseriesCoProcessor> tsp4(new TimeseriesCoProcessor(grid, stepMV, mic4, pathOut+"/mic/mic4", comm));

      SPtr<IntegrateValuesHelper> mic5(new IntegrateValuesHelper(grid, comm, 0.3+deltaXfine, 0.015, 0.000517+0.00037+7.0*deltaXfine, 0.3+2.0*deltaXfine, 0.015+deltaXfine, 0.000517+0.00037+8.0*deltaXfine));
      if (myid==0) GbSystem3D::writeGeoObject(mic5->getBoundingBox().get(), pathOut+"/geo/mic5", WbWriterVtkXmlBinary::getInstance());
      SPtr<TimeseriesCoProcessor> tsp5(new TimeseriesCoProcessor(grid, stepMV, mic5, pathOut+"/mic/mic5", comm));

      SPtr<IntegrateValuesHelper> mic6(new IntegrateValuesHelper(grid, comm, 0.3, 0.015, -0.46+4.25*deltaXcoarse, 0.3+deltaXcoarse, 0.015+deltaXcoarse, -0.46+5.25*deltaXcoarse));
      if (myid==0) GbSystem3D::writeGeoObject(mic6->getBoundingBox().get(), pathOut+"/geo/mic6", WbWriterVtkXmlBinary::getInstance());
      SPtr<TimeseriesCoProcessor> tsp6(new TimeseriesCoProcessor(grid, stepMV, mic6, pathOut+"/mic/mic6", comm));

      omp_set_num_threads(numOfThreads);
      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      //calculator->addCoProcessor(restartCoProcessor);
      calculator->addCoProcessor(migCoProcessor);
      calculator->addCoProcessor(writeMQSelectCoProcessor);
      calculator->addCoProcessor(writeMQCoProcessor);
      calculator->addCoProcessor(tsp1);
      calculator->addCoProcessor(tsp2);
      calculator->addCoProcessor(tsp3);
      calculator->addCoProcessor(tsp4);
      calculator->addCoProcessor(tsp5);
      calculator->addCoProcessor(tsp6);
      calculator->addCoProcessor(tav);


      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      if (myid==0) UBLOG(logINFO, "Simulation-end");

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

   //test_run();

   //SuperMUC
   //MPI_Finalize();

   return 0;
}


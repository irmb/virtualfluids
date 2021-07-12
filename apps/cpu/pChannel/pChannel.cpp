#include <iostream>
#include <string>
#include "VirtualFluids.h"
#include <omp.h>
double rangeRandom(double M, double N)
{
   return M + (rand() / (RAND_MAX / (N - M)));
}

double rangeRandom1()
{
   return (2.0*rand())/RAND_MAX - 1.0;
}

//double rangeRandom(double M, double N)
//{
//   return rand() % (int)N+(int)M;
//}



//#include <thread>

using namespace std;

std::vector<int> x1Nbr;
std::vector<int> x2Nbr;
std::vector<int> x3Nbr;

std::vector<int> x1NbrTemp;
std::vector<int> x2NbrTemp;
std::vector<int> x3NbrTemp;

//void findSolidNeighbor(SPtr<GbVoxelMatrix3D> voxelMatrix, int x1, int x2, int x3)
//{
//   for (int k3 = -1; k3<=1; k3++)
//   {
//      for (int k2 = -1; k2<=1; k2++)
//      {
//         for (int k1 = -1; k1<=1; k1++)
//         {
//            int j1 = x1+k1;
//            int j2 = x2+k2;
//            int j3 = x3+k3;
//            if (j1>=0&&j1<nodesX1 && j2>=0&&j2<nodesX2 && j3>=0&&j3<nodesX3)
//            {
//               if ((*voxelMatrix)(j1, j2, j3)==GbVoxelMatrix3D::FLUID)
//               {
//                  if (flagMatrix(j1, j2, j3)==0)
//                  {
//                     voxelMatrixTemp(j1, j2, j3) = GbVoxelMatrix3D::SOLID;
//                     flagMatrix(j1, j2, j3) = 1;
//                     x1NbrTemp.push_back(j1);
//                     x2NbrTemp.push_back(j2);
//                     x3NbrTemp.push_back(j3);
//                  }
//               }
//            }
//         }
//      }
//   }
//}

void changePorosity(SPtr<GbVoxelMatrix3D> sample, vector<int> pmNX)
{
   int minX1 = 0;
   int minX2 = 0;
   int minX3 = 0;
   int maxX1 = pmNX[0];
   int maxX2 = pmNX[1];
   int maxX3 = pmNX[2];
   sample->calculateNumberOfSolidAndFluid();
   double nSolid = sample->getNumberOfSolid();
   double nFluid = sample->getNumberOfFluid();
   double porosityStart = nFluid/(nSolid+nFluid);
   double porosityEnd = 0.5;
   double porosityStep = (porosityEnd-porosityStart)/(double)maxX3;
   double totallSliceVoxel = maxX1*maxX2;
   vector<int> fluidThreshold;

   SPtr<GbVoxelMatrix3D> sampleTemp = SPtr<GbVoxelMatrix3D>(sample->clone());

   int count=1;

   for (int ix3=minX3; ix3<maxX3; ix3++)
   {
      int cFluid = 0;
      for (int ix2=minX2; ix2<maxX2; ix2++)
      {
         for (int ix1=minX1; ix1<maxX1; ix1++)
         {
            if ((*sample)(ix1, ix2, ix3) == GbVoxelMatrix3D::FLUID)
            {
               cFluid++;
            }
         }
      }
      double slicePorosity = (double)cFluid/totallSliceVoxel;
      double porosityPercent = (porosityStep*(double)count)/slicePorosity;
      fluidThreshold.push_back((totallSliceVoxel-(double)cFluid)*porosityPercent);
      count++;
   }
   int solidCount = 0;
   count=0;

   for (int ix3=minX3; ix3<maxX3; ix3++)
   {
      //while (fluidThreshold[count] > 0)
      //{
     // int fTh = fluidThreshold[count];

         int solidCount = 0;
         for (int ix2=minX2; ix2<maxX2; ix2++)
         {
            for (int ix1=minX1; ix1<maxX1; ix1++)
            {
               if ((*sample)(ix1, ix2, ix3) == GbVoxelMatrix3D::SOLID)
               {
                  bool flagSolid = true;
                  for (int k3 = -1; k3<=1; k3++)
                  {
                     for (int k2 = -1; k2<=1; k2++)
                     {
                        for (int k1 = -1; k1<=1; k1++)
                        {
                           int j1 = ix1+k1;
                           int j2 = ix2+k2;
                           int j3 = ix3;//+k3;
                           if (j1>=0&&j1<maxX1 && j2>=0&&j2<maxX2 && j3>=0&&j3<maxX3)
                           {
                              if ((*sample)(j1, j2, j3) == GbVoxelMatrix3D::FLUID)
                              {
                                 flagSolid = flagSolid && false;
                              }
                           }
                        }
                     }
                  }
                  if (!flagSolid)
                  {
                     (*sample)(ix1, ix2, ix3)=GbVoxelMatrix3D::FLUID;
                      fluidThreshold[count]--;
                     solidCount++;
                  }
                  if ( fluidThreshold[count] == 0)
                  {
                     ix1=maxX1;
                     ix2=maxX2;
                  }
               }
            }
         }
         UBLOG(logINFO, "count = " << count);
         UBLOG(logINFO, "fTh = " <<  fluidThreshold[count]);
         UBLOG(logINFO, "solidCount = " << solidCount);
         //sample = sampleTemp;
         //sampleTemp = SPtr<GbVoxelMatrix3D>(sample->clone());
         
        count++;     
      
       }
      
   //}
   //sampleTemp->writeToLegacyVTKBinary("d:/temp/ChannelFlow/geo/sampleTemp.vtk");
}

//////////////////////////////////////////////////////////////////////////
void run(string configname)
{
   try
   {
      vf::basics::ConfigurationFile   config;
      config.load(configname);

      string          pathOut           = config.getValue<string>("pathOut");
      string          pathGeo           = config.getValue<string>("pathGeo");
      int             numOfThreads      = config.getValue<int>("numOfThreads");
      string          sampleFilename    = config.getValue<string>("sampleFilename");
      vector<int>     pmNX              = config.getVector<int>("pmNX");
      double          lthreshold        = config.getValue<double>("lthreshold");
      double          uthreshold        = config.getValue<double>("uthreshold");
      vector<float>   voxelDeltaX       = config.getVector<float>("voxelDeltaX");
      vector<int>     blocknx           = config.getVector<int>("blocknx");
      double          u_LB              = config.getValue<double>("u_LB");
      double          restartStep       = config.getValue<double>("restartStep");
      double          cpStep            = config.getValue<double>("cpStep");
      double          cpStart           = config.getValue<double>("cpStart");
      double          endTime           = config.getValue<double>("endTime");
      double          outTime           = config.getValue<double>("outTime");
      double          availMem          = config.getValue<double>("availMem");
      bool            rawFile           = config.getValue<bool>("rawFile");
      bool            logToFile         = config.getValue<bool>("logToFile");
      bool            writeSample       = config.getValue<bool>("writeSample");
      vector<double>  pmL               = config.getVector<double>("pmL");
      double          deltaXfine        = config.getValue<double>("deltaXfine");
      int             refineLevel       = config.getValue<int>("refineLevel");
      bool            thinWall          = config.getValue<bool>("thinWall");
      double          Re                = config.getValue<double>("Re");
      double          channelHigh       = config.getValue<double>("channelHigh");
      bool            changeQs          = config.getValue<bool>("changeQs");
      double          timeAvStart       = config.getValue<double>("timeAvStart");
      double          timeAvStop        = config.getValue<double>("timeAvStop");
      bool            averaging         = config.getValue<bool>("averaging");
      bool            averagingReset    = config.getValue<bool>("averagingReset");
      bool            newStart          = config.getValue<bool>("newStart");
      vector<double>  nupsStep          = config.getVector<double>("nupsStep");
      vector<double>  boundingBox       = config.getVector<double>("boundingBox");

      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

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

      //Sleep(30000);

      if (myid == 0) UBLOG(logINFO, "Testcase porous channel");

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;
      double deltaXcoarse = deltaXfine*(double)(1<<refineLevel);

      LBMReal rho_LB = 0.0;
      double rhoReal = 1.2041; //(kg/m3)
      double uReal = 48; //m/s
      double lReal = 0.008;//m
      double hLB = lReal / deltaXcoarse;
      double Ma = 0.13;//Ma-Real!
      double csReal = uReal / Ma;
      LBMUnitConverter unitConverter(lReal, csReal, rhoReal, hLB);
      if (myid==0) UBLOG(logINFO, unitConverter.toString());

      //double coord[6];

      vector<double> origin(3);
      origin[0] = 0;
      origin[1] = 0;
      origin[2] = 0;

      //real velocity is 49.63 m/s

      SPtr<Grid3D> grid(new Grid3D(comm));

      //BC adapters
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);

      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new BCProcessor());

      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel());

      mu::Parser fctForcingX1;
      fctForcingX1.SetExpr("Fx1");
      fctForcingX1.DefineConst("Fx1", 1.0e-6);
      kernel->setWithForcing(true);

      kernel->setBCProcessor(bcProc);

      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIORestartCoProcessor> restartCoProcessor(new MPIIORestartCoProcessor(grid, rSch, pathOut, comm));
      restartCoProcessor->setLBMKernel(kernel);
      restartCoProcessor->setBCProcessor(bcProc);

      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIOMigrationCoProcessor> migCoProcessor(new MPIIOMigrationCoProcessor(grid, mSch, pathOut+"/mig", comm));
      migCoProcessor->setLBMKernel(kernel);
      migCoProcessor->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////

      //bounding box
      double g_minX1 = boundingBox[0];
      double g_minX2 = boundingBox[1];
      double g_minX3 = boundingBox[2];

      double g_maxX1 = boundingBox[3];
      double g_maxX2 = boundingBox[4];
      double g_maxX3 = boundingBox[5];

      double blockLength = (double)blocknx[0]*deltaXcoarse;

      double channel_high = channelHigh; // g_maxX3-g_minX3;
      double channel_high_LB = channel_high/deltaXcoarse;
      //////////////////////////////////////////////////////////////////////////
      double nu_LB = (u_LB*channel_high_LB)/Re;
      //////////////////////////////////////////////////////////////////////////
      if (myid == 0)
      {
         UBLOG(logINFO, "Parameters:");
         UBLOG(logINFO, "Re                  = " << Re);
         UBLOG(logINFO, "u_LB                = " << u_LB);
         UBLOG(logINFO, "rho_LB              = " << rho_LB);
         UBLOG(logINFO, "nu_LB               = " << nu_LB);
         UBLOG(logINFO, "dx coarse           = " << deltaXcoarse << " m");
         UBLOG(logINFO, "dx fine             = " << deltaXfine << " m");
         UBLOG(logINFO, "channel_high        = " << channel_high << " m");
         UBLOG(logINFO, "channel_high_LB     = " << channel_high_LB);
         UBLOG(logINFO, "number of levels    = " << refineLevel + 1);
         UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
         UBLOG(logINFO, "number of threads   = " << numOfThreads);
         UBLOG(logINFO, "path = " << pathOut);
         UBLOG(logINFO, "Preprocess - start");
      }


      if (newStart)
      {
         if (myid == 0) UBLOG(logINFO, "new start...");



         grid->setPeriodicX1(true);
         grid->setPeriodicX2(true);
         grid->setPeriodicX3(false);
         grid->setDeltaX(deltaXcoarse);
         grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

         SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);


         //////////////////////////////////////////////////////////////////////////
         //refinement
         double blockLengthX3Fine = grid->getDeltaX(refineLevel) * blocknx[2];
         double refHight = 0.002;

         GbCuboid3DPtr refineBoxTop(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3-refHight, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxTop.get(), pathOut + "/geo/refineBoxTop", WbWriterVtkXmlASCII::getInstance());

         //GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+offsetMinX3+blockLengthX3Fine));
         GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+refHight));
         if (myid == 0) GbSystem3D::writeGeoObject(refineBoxBottom.get(), pathOut + "/geo/refineBoxBottom", WbWriterVtkXmlASCII::getInstance());

         if (refineLevel > 0)
         {
            if (myid == 0) UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
            refineHelper.addGbObject(refineBoxTop, refineLevel);
            refineHelper.addGbObject(refineBoxBottom, refineLevel);
            refineHelper.refine();
            if (myid == 0) UBLOG(logINFO, "Refinement - end");
         }
         //////////////////////////////////////////////////////////////////////////

         //walls
         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathOut+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathOut+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());


         //wall interactors
         SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
         ////////////////////////////////////////////
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

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

         //////////////////////////////////
         //undef nodes for refinement
         if (refineLevel > 0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }


         //BC
         intHelper.setBC();

         ////porous media
         if (true)
         {
            string samplePathname = pathGeo + "/" + sampleFilename;

            SPtr<GbVoxelMatrix3D> sample(new GbVoxelMatrix3D(pmNX[0], pmNX[1], pmNX[2], 0, lthreshold, uthreshold));
            if (rawFile)
            {
               sample->readBufferedMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
            }
            else
            {
               sample->readMatrixFromVtiASCIIFile(samplePathname);
            }

            sample->setVoxelMatrixDelta(voxelDeltaX[0], voxelDeltaX[1], voxelDeltaX[2]);
            sample->setVoxelMatrixMininum(origin[0], origin[1], origin[2]);

            changePorosity(sample, pmNX);

            sample->writeToLegacyVTKBinary(pathOut+"/geo/sample.vtk");
            return;

            int bounceBackOption = 1;
            bool vxFile = false;
            int i = 0;
            //for (int x = 0; x < lengthFactor; x+=2)
            int lenX = (int)((g_maxX1-g_minX1)/(pmL[0]));
            int lenY = (int)((g_maxX2-g_minX2)/(pmL[1]));

            for (int y = 0; y < lenY; y+=2)
               for (int x = 0; x < lenX; x+=2)
               {
                  double offsetX = pmL[0] * (double)x;
                  double offsetY = pmL[1] * (double)y;
                  //sample 0
                  if (myid == 0) UBLOG(logINFO, "sample # " << i);
                  sample->setVoxelMatrixMininum(origin[0]+offsetX, origin[1]+offsetY, origin[2]);
                  Utilities::voxelMatrixDiscretisation(sample, pathOut, myid, i, grid, bounceBackOption, vxFile);
                  i++;

                  if (myid == 0)
                  {
                     UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()/1073741824.0<<" GB");
                     UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed()/1073741824.0<<" GB");
                     UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
                  }

                  //sample 1
                  if (myid == 0) UBLOG(logINFO, "sample # " << i);
                  sample->setVoxelMatrixMininum(origin[0]+pmL[0]+offsetX, origin[1]+offsetY, origin[2]);
                  sample->mirrorX();
                  Utilities::voxelMatrixDiscretisation(sample, pathOut, myid, i, grid, bounceBackOption, vxFile);
                  i++;

                  if (myid == 0)
                  {
                     UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()/1073741824.0<<" GB");
                     UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed()/1073741824.0<<" GB");
                     UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
                  }

                  //sample 2
                  if (myid == 0) UBLOG(logINFO, "sample # " << i);
                  sample->setVoxelMatrixMininum(origin[0]+pmL[0]+offsetX, origin[1]+pmL[1]+offsetY, origin[2]);
                  sample->mirrorY();
                  Utilities::voxelMatrixDiscretisation(sample, pathOut, myid, i, grid, bounceBackOption, vxFile);
                  i++;

                  if (myid == 0)
                  {
                     UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()/1073741824.0<<" GB");
                     UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed()/1073741824.0<<" GB");
                     UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
                  }

                  //sample 3
                  if (myid == 0) UBLOG(logINFO, "sample # " << i);
                  sample->setVoxelMatrixMininum(origin[0]+offsetX, origin[1]+pmL[1]+offsetY, origin[2]);
                  sample->mirrorX();
                  Utilities::voxelMatrixDiscretisation(sample, pathOut, myid, i, grid, bounceBackOption, vxFile);
                  sample->mirrorY();
                  i++;
               }

         }

         grid->accept(bcVisitor);

         mu::Parser inflowProfileVx1, inflowProfileVx2, inflowProfileVx3, inflowProfileRho;
         inflowProfileVx1.SetExpr("x3 < h ? 0.0 : uLB+1*x1");
         inflowProfileVx1.DefineConst("uLB", u_LB);
         inflowProfileVx1.DefineConst("h", pmL[2]);

         InitDistributionsBlockVisitor initVisitor;
         initVisitor.setVx1(inflowProfileVx1);
         grid->accept(initVisitor);

         ////set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //Postrozess
         {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
            ppgeo.process(0);
         }

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         //restartCoProcessor->restart((int)restartStep);
         migCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);
         ////////////////////////////////////////////////////////////////////////////
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept(setConnsVisitor);

         grid->accept(bcVisitor);

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }

      SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      std::shared_ptr<NUPSCounterCoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> stepSch(new UbScheduler(outTime));

      SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      SPtr<GbObject3D> bbBox(new GbCuboid3D(g_minX1-blockLength, (g_maxX2-g_minX2)/2.0, g_minX3-blockLength, g_maxX1+blockLength, (g_maxX2-g_minX2)/2.0+deltaXcoarse, g_maxX3+blockLength));
      if (myid==0) GbSystem3D::writeGeoObject(bbBox.get(), pathOut+"/geo/bbBox", WbWriterVtkXmlASCII::getInstance());
      SPtr<WriteMQFromSelectionCoProcessor> writeMQSelectCoProcessor(new WriteMQFromSelectionCoProcessor(grid, stepSch, bbBox, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));


      SPtr<UbScheduler> AdjForcSch(new UbScheduler());
      AdjForcSch->addSchedule(10, 0, 10000000);
      SPtr<IntegrateValuesHelper> intValHelp(new IntegrateValuesHelper(grid, comm, g_minX1, g_minX2, g_minX3+pmL[2], g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(intValHelp->getBoundingBox().get(), pathOut + "/geo/IntValHelp", WbWriterVtkXmlBinary::getInstance());

      double vxTarget=u_LB;
      SPtr<AdjustForcingCoProcessor> AdjForcCoProcessor(new AdjustForcingCoProcessor(grid, AdjForcSch, pathOut, intValHelp, vxTarget, comm));

      //mu::Parser decrViscFunc;
      //decrViscFunc.SetExpr("nue0+c0/(t+1)/(t+1)");
      //decrViscFunc.DefineConst("nue0", nu_LB*4.0);
      //decrViscFunc.DefineConst("c0", 0.1);
      //SPtr<UbScheduler> DecrViscSch(new UbScheduler());
      //DecrViscSch->addSchedule(10, 0, 1000);
      //DecreaseViscosityCoProcessor decrViscPPPtr(grid, DecrViscSch, &decrViscFunc, comm);

     //if (changeQs)
     //{
       // double z1 = pmL[2];
       // SPtr<IntegrateValuesHelper> intValHelp2(new IntegrateValuesHelper(grid, comm,
         //  coord[0], coord[1], z1 - deltaXfine,
         //  coord[3], coord[4], z1 + deltaXfine));
       // if (myid == 0) GbSystem3D::writeGeoObject(intValHelp2->getBoundingBox().get(), pathOut + "/geo/intValHelp2", WbWriterVtkXmlBinary::getInstance());
       // Utilities::ChangeRandomQs(intValHelp2);
     //}

      std::vector<double> levelCoords;
      std::vector<int> levels;
      std::vector<double> bounds;
      //bounds.push_back(0);
      //bounds.push_back(0);
      //bounds.push_back(0);
      //bounds.push_back(0.004);
      //bounds.push_back(0.002);
      //bounds.push_back(0.003);
      //levels.push_back(1);
      //levels.push_back(0);
      //levels.push_back(1);
      //levelCoords.push_back(0);
      //levelCoords.push_back(0.0016);
      //levelCoords.push_back(0.0024);
      //levelCoords.push_back(0.003);
      bounds.push_back(0);
      bounds.push_back(0);
      bounds.push_back(0);
      bounds.push_back(0.004);
      bounds.push_back(0.002);
      bounds.push_back(0.002);
      levels.push_back(0);
      levelCoords.push_back(0);
      levelCoords.push_back(0.002);
      //SPtr<UbScheduler> tavSch(new UbScheduler(1, timeAvStart, timeAvStop));
      //SPtr<CoProcessor> tav(new TimeAveragedValuesCoProcessor(grid, pathOut, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,
      //   TimeAveragedValuesCoProcessor::Velocity | TimeAveragedValuesCoProcessor::Fluctuations | TimeAveragedValuesCoProcessor::Triplecorrelations,
      //   levels, levelCoords, bounds));


      //create line time series
      SPtr<UbScheduler> tpcSch(new UbScheduler(1, 1, 3));
      //GbPoint3DPtr p1(new GbPoint3D(0.0,0.005,0.01));
      //GbPoint3DPtr p2(new GbPoint3D(0.064,0.005,0.01));
      //SPtr<GbLine3D> line(new GbLine3D(p1.get(),p2.get()));
      SPtr<GbLine3D> line(new GbLine3D(new GbPoint3D(0.0, 0.005, 0.01), new GbPoint3D(0.064, 0.005, 0.01)));
      LineTimeSeriesCoProcessor lineTs(grid, tpcSch, pathOut+"/TimeSeries/line1.csv", line, 0, comm);
      if (myid==0) lineTs.writeLine(pathOut+"/geo/line1");

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      //omp_set_num_threads(numOfThreads);
       numOfThreads = 4;
      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, (int)endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(AdjForcCoProcessor);
      calculator->addCoProcessor(migCoProcessor);
      //calculator->addCoProcessor(restartCoProcessor);
      calculator->addCoProcessor(writeMQSelectCoProcessor);
      calculator->addCoProcessor(writeMQCoProcessor);

      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      if (myid == 0) UBLOG(logINFO, "Simulation-end");
   }
   catch (exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch (string& s)
   {
      cerr << s << endl;
   }
   catch (mu::Parser::exception_type &e)
   {
      std::cout << e.GetMsg() << std::endl;
   }
   catch (...)
   {
      cerr << "unknown exception" << endl;
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
         cout << "Configuration file is missing!" << endl;
      }
   }

   return 0;
}

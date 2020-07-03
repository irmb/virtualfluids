#include <iostream>
#include <string>

#include <VirtualFluids.h>


using namespace std;

void changeDP()
{
}
//////////////////////////////////////////////////////////////////////////
void run(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathname = config.getString("pathname");
      string          pathGeo = config.getString("pathGeo");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      string          sampleFilename = config.getString("sampleFilename");
      int             pmNX1 = config.getValue<int>("pmNX1");
      int             pmNX2 = config.getValue<int>("pmNX2");
      int             pmNX3 = config.getValue<int>("pmNX3");
      double          lthreshold = config.getValue<double>("lthreshold");
      double          uthreshold = config.getValue<double>("uthreshold");
      //double          pmL1 = config.getValue<double>("pmL1");
      //double          pmL2 = config.getValue<double>("pmL2");
      //double          pmL3 = config.getValue<double>("pmL3");
      int             blocknx = config.getValue<int>("blocknx");
      //double          nx3 = config.getValue<double>("nx3");
      double          dp_LB = config.getValue<double>("dp_LB");
      double          nu_LB = config.getValue<double>("nu_LB");
      string          timeSeriesFile = config.getString("timeSeriesFile");
      double          restartStep = config.getValue<double>("restartStep");
      //double          restartStepStart = config.getValue<double>("restartStepStart");
      int          endTime = config.getValue<int>("endTime");
      double          outTimeStep = config.getValue<double>("outTimeStep");
      double          outTimeStart = config.getValue<double>("outTimeStart");
      double          availMem = config.getValue<double>("availMem");
      bool            rawFile = config.getValue<bool>("rawFile");
      double          timeSeriesOutTime = config.getValue<double>("timeSeriesOutTime");
      bool            logToFile = config.getValue<bool>("logToFile");
      bool            spongeLayer = config.getValue<bool>("spongeLayer");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");
      double          deltax = config.getValue<double>("deltax");
      bool            newViscosity = config.getValue<bool>("newViscosity");
      bool            newPressure = config.getValue<bool>("newPressure");
      //bool            pmDeltas = config.getValue<bool>("pmDeltas");
      double          pmDeltaX1 = config.getValue<double>("pmDeltaX1");
      double          pmDeltaX2 = config.getValue<double>("pmDeltaX2");
      double          pmDeltaX3 = config.getValue<double>("pmDeltaX3");
      double          vx1 = config.getValue<double>("vx1");
      double          vx2 = config.getValue<double>("vx2");
      double          vx3 = config.getValue<double>("vx3");
      bool            yDir = config.getValue<bool>("yDir");
      bool            zDir = config.getValue<bool>("zDir");
      double          cpStep = config.getValue<double>("cpStep");
      double          cpStepStart = config.getValue<double>("cpStepStart");
      bool            newStart = config.getValue<bool>("newStart");
      int             chunk = config.getValue<int>("chunk");


      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (logToFile)
      {
#if defined(__unix__)
         if (myid==0)
         {
            const char* str = pathname.c_str();
            int status = mkdir(str, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
         }
#endif 

         if (myid==0)
         {
            stringstream logFilename;
            logFilename<<pathname+"/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }

      if (myid==0) UBLOG(logINFO, "Testcase permeability");

      //string machinename = UbSystem::getMachineName();
      //UBLOG(logINFO, "PID = " << myid << " Hostname: " << machinename);
      //UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
      //UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
      //UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());

      int blocknx1 = blocknx;
      int blocknx2 = blocknx;
      int blocknx3 = blocknx;

      LBMReal rho_LB = 0.0;
      double rhoLBinflow = dp_LB*3.0;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      double coord[6];
      //double deltax;

      SPtr<Grid3D> grid(new Grid3D(comm));

      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel());
      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new ThinWallBCProcessor());
      kernel->setBCProcessor(bcProc);

      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStepStart));
      SPtr<MPIIORestartCoProcessor> restartCoProcessor(new MPIIORestartCoProcessor(grid, rSch, pathname, comm));
      restartCoProcessor->setLBMKernel(kernel);
      restartCoProcessor->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////

      //BC Adapter
      //////////////////////////////////////////////////////////////////////////////
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new ThinWallNoSlipBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapterInflow(new DensityBCAdapter(rhoLBinflow));
      denBCAdapterInflow->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapterOutflow(new DensityBCAdapter(rho_LB));
      denBCAdapterOutflow->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));
      //////////////////////////////////////////////////////////////////////////////////
      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(denBCAdapterInflow);
      bcVisitor.addBC(denBCAdapterOutflow);

      if (newStart)
      {
         if (myid==0) UBLOG(logINFO, "new start..");
         if (myid==0) UBLOG(logINFO, "preprocess start..");

         //UBLOG(logINFO, "new start PID = " << myid << " Hostname: " << machinename);
         //UBLOG(logINFO, "new start PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         //UBLOG(logINFO, "new start PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         //UBLOG(logINFO, "new start PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());

         string samplePathname = pathGeo+sampleFilename;

         double deltaVoxelX1 = pmDeltaX1;
         double deltaVoxelX2 = pmDeltaX2;
         double deltaVoxelX3 = pmDeltaX3;

         if (myid==0) UBLOG(logINFO, "read voxel matrix: start");
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////         
         SPtr<GbVoxelMatrix3D> sample1(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0, lthreshold, uthreshold));
         if (rawFile)
         {
            sample1->readMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
         }
         else
         {
            sample1->readMatrixFromVtiASCIIFile(samplePathname);
         }
         sample1->setVoxelMatrixDelta((float)deltaVoxelX1, (float)deltaVoxelX2, (float)deltaVoxelX3);
         sample1->setVoxelMatrixMininum(0.0, 0.0, 0.0);

         if (myid==0) UBLOG(logINFO, "sample1: rotate voxel matrix: start");
         if (yDir)
         {
            sample1->rotate90aroundZ();
         }
         else if (zDir)
         {
            sample1->rotate90aroundY();
         }
         if (myid==0) UBLOG(logINFO, "sample1: rotate voxel matrix: end");
         
         if (myid==0) sample1->writeToVTKImageDataAppended(pathname+"/geo/sample1");
///////////////////////////////////////////////////////////////////////////////////////////////////////////         
         SPtr<GbVoxelMatrix3D> sample2(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0, lthreshold, uthreshold));
         if (rawFile)
         {
            sample2->readMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
         }
         else
         {
            sample2->readMatrixFromVtiASCIIFile(samplePathname);
         }
         sample2->setVoxelMatrixDelta((float)deltaVoxelX1, (float)deltaVoxelX2, (float)deltaVoxelX3);
         sample2->setVoxelMatrixMininum(0.0, sample1->getX2Maximum(), 0.0);
         if (myid==0) UBLOG(logINFO, "read voxel matrix: end");

         if (myid==0) UBLOG(logINFO, "sample2: rotate voxel matrix: start");
         if (yDir)
         {
            sample2->rotate90aroundZ();
            sample2->setVoxelMatrixMinX2(sample1->getX2Maximum());
         }
         else if (zDir)
         {
            sample2->rotate90aroundY();
         }

         sample2->mirrorY();

         if (myid==0) UBLOG(logINFO, "sample2: rotate voxel matrix: end");

         if (myid==0) sample2->writeToVTKImageDataAppended(pathname+"/geo/sample2");

////////////////////////////////////////////////////////////////////////////////////////////////////////////
         SPtr<GbVoxelMatrix3D> sample3(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0, lthreshold, uthreshold));
         if (rawFile)
         {
            sample3->readMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
         }
         else
         {
            sample3->readMatrixFromVtiASCIIFile(samplePathname);
         }
         sample3->setVoxelMatrixDelta((float)deltaVoxelX1, (float)deltaVoxelX2, (float)deltaVoxelX3);
         sample3->setVoxelMatrixMininum(0.0, 0.0, sample1->getX3Maximum());
         if (myid == 0) UBLOG(logINFO, "read voxel matrix: end");

         if (myid==0) UBLOG(logINFO, "sample3: rotate voxel matrix: start");
         if (yDir)
         {
            sample3->rotate90aroundZ();
         }
         else if (zDir)
         {
            sample3->rotate90aroundY();
         }
         sample3->mirrorZ();
         if (myid==0) UBLOG(logINFO, "sample3: rotate voxel matrix: end");

         if (myid==0) sample3->writeToVTKImageDataAppended(pathname+"/geo/sample3");

         ////////////////////////////////////////////////////////////////////////////////////////////////////////////         
         SPtr<GbVoxelMatrix3D> sample4(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0, lthreshold, uthreshold));
         if (rawFile)
         {
            sample4->readMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
         }
         else
         {
            sample4->readMatrixFromVtiASCIIFile(samplePathname);
         }
         sample4->setVoxelMatrixDelta((float)deltaVoxelX1, (float)deltaVoxelX2, (float)deltaVoxelX3);
         sample4->setVoxelMatrixMininum(0.0, sample1->getX1Maximum(), sample1->getX3Maximum());
         if (myid == 0) UBLOG(logINFO, "read voxel matrix: end");

         if (myid==0) UBLOG(logINFO, "sample4: rotate voxel matrix: start");

         if (yDir)
         {
            sample4->rotate90aroundZ();
            sample4->setVoxelMatrixMinX2(sample1->getX2Maximum());
            sample4->setVoxelMatrixMinX3(sample1->getX3Maximum());
         }
         else if (zDir)
         {
            sample4->rotate90aroundY();
         }
         sample4->mirrorY();
         sample4->mirrorZ();
         if (myid==0) UBLOG(logINFO, "sample4: rotate voxel matrix: end");

         if (myid==0) sample4->writeToVTKImageDataAppended(pathname+"/geo/sample4");

         ///////////////////////////////////////////////////////

         ////////////////////////////////////////////////////////////////////////

         double offset1 = sample1->getLengthX1()/10.0;
         double offset2 = 2.0*offset1;
         //double offset2 = offset1;
         //bounding box
         double g_minX1 = sample1->getX1Minimum()-offset1;
         double g_minX2 = sample1->getX2Minimum()-0.5*deltax;
         double g_minX3 = sample1->getX3Minimum()-0.5*deltax;

         double g_maxX1 = sample1->getX1Maximum()+offset2;
         double g_maxX2 = sample4->getX2Maximum()-0.5*deltax;
         double g_maxX3 = sample4->getX3Maximum()-0.5*deltax;

         if (myid==0)
         {
            UBLOG(logINFO, "g_minX1="<<g_minX1<<",g_minX2="<<g_minX2<<",g_minX3="<<g_minX3<<",g_maxX1="<<g_maxX1<<",g_maxX2="<<g_maxX2<<",g_maxX3="<<g_maxX3);
         }

         double blockLength = (double)blocknx1*deltax;

         grid->setPeriodicX1(false);
         grid->setPeriodicX2(true);
         grid->setPeriodicX3(true);
         grid->setDeltaX(deltax);
         grid->setBlockNX(blocknx1, blocknx2, blocknx3);

         SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid==0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         if (myid==0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "rho_LB = "<<rho_LB);
            UBLOG(logINFO, "nu_LB = "<<nu_LB);
            UBLOG(logINFO, "dp_LB = "<<dp_LB);
            UBLOG(logINFO, "dx = "<<deltax<<" m");
            UBLOG(logINFO, "numOfThreads = "<<numOfThreads);
            UBLOG(logINFO, "path = "<<pathname);
            UBLOG(logINFO, "Preprozess - start");
         }

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());
         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         SPtr<WriteBlocksCoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         //PM interactor
         SPtr<D3Q27Interactor> sample1Int(new D3Q27Interactor(sample1, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> sample2Int(new D3Q27Interactor(sample2, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> sample3Int(new D3Q27Interactor(sample3, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> sample4Int(new D3Q27Interactor(sample4, grid, noSlipBCAdapter, Interactor3D::SOLID));
         //inflow
         SPtr<D3Q27Interactor> inflowInt(new D3Q27Interactor(geoInflow, grid, denBCAdapterInflow, Interactor3D::SOLID));
         //outflow
         SPtr<D3Q27Interactor> outflowInt(new D3Q27Interactor(geoOutflow, grid, denBCAdapterOutflow, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::RECURSIVE));
         ////////////////////////////////////////////

         /////delete solid blocks
         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.addInteractor(sample1Int);
         intHelper.addInteractor(sample2Int);
         intHelper.addInteractor(sample3Int);
         intHelper.addInteractor(sample4Int);
         intHelper.selectBlocks();
         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         ppblocks->process(0);
         ppblocks.reset();


         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx1)* (blocknx2)* (blocknx3);
         unsigned long nod = nob * (blocknx1)* (blocknx2)* (blocknx3);
         unsigned long nodg = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);
         double needMemAll = double(nodg*(27*sizeof(double)+sizeof(int)+sizeof(float)*4));
         double needMem = needMemAll/double(comm->getNumberOfProcesses());

         if (myid==0)
         {
            UBLOG(logINFO, "Number of blocks = "<<nob);
            UBLOG(logINFO, "Number of nodes  = "<<nod);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for (int level = minInitLevel; level<=maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO, "Number of blocks for level "<<level<<" = "<<nobl);
               UBLOG(logINFO, "Number of nodes for level "<<level<<" = "<<nobl*nodb);
            }
            UBLOG(logINFO, "Necessary memory  = "<<needMemAll<<" bytes");
            UBLOG(logINFO, "Necessary memory per process = "<<needMem<<" bytes");
            UBLOG(logINFO, "Available memory per process = "<<availMem<<" bytes");
         }

         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nu_LB, availMem, needMem);
         grid->accept(kernelVisitor);

         //BC
         intHelper.setBC();

         //BS visitor
         grid->accept(bcVisitor);

         //Press*1.6e8+(14.76-coordsX)/3.5*5000
         //initialization of distributions
         mu::Parser fct;
         fct.SetExpr("(x1max-x1)/l*dp*3.0");
         fct.DefineConst("dp", dp_LB);
         fct.DefineConst("x1max", g_maxX1);
         fct.DefineConst("l", g_maxX1-g_minX1);

         InitDistributionsBlockVisitor initVisitor;
         initVisitor.setRho(fct);
         grid->accept(initVisitor);

         //Postrozess
         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         SPtr<WriteBoundaryConditionsCoProcessor> ppgeo(
            new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         ppgeo->process(0);
         ppgeo.reset();

         coord[0] = sample1->getX1Minimum();
         coord[1] = sample1->getX2Minimum();
         coord[2] = sample1->getX3Minimum();
         coord[3] = sample4->getX1Maximum();
         coord[4] = sample4->getX2Maximum();
         coord[5] = sample4->getX3Maximum();

         ////////////////////////////////////////////////////////
         UbFileOutputASCII outf(pathname+"/checkpoints/coord.txt");
         outf.writeDouble(deltax);
         outf.writeDouble(coord[0]);
         outf.writeDouble(coord[1]);
         outf.writeDouble(coord[2]);
         outf.writeDouble(coord[3]);
         outf.writeDouble(coord[4]);
         outf.writeDouble(coord[5]);
         outf.writeDouble(g_minX1);
         outf.writeDouble(g_maxX1);
         outf.writeDouble(availMem);
         outf.writeDouble(needMem);
         ////////////////////////////////////////////////////////

         grid->addInteractor(inflowInt);

         if (myid==0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         ////////////////////////////////////////////////////////
         UbFileInputASCII inf(pathname+"/checkpoints/coord.txt");
         deltax = inf.readDouble();
         coord[0] = inf.readDouble();
         coord[1] = inf.readDouble();
         coord[2] = inf.readDouble();
         coord[3] = inf.readDouble();
         coord[4] = inf.readDouble();
         coord[5] = inf.readDouble();
         double g_minX1 = inf.readDouble();
         double g_maxX1 = inf.readDouble();
         double availMem = inf.readDouble();
         double needMem = inf.readDouble();
         ////////////////////////////////////////////////////////

         restartCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);


         //new nu
         //if (newViscosity)
         //{
         //   ViscosityBlockVisitor nuVisitor(nu_LB);
         //   grid->accept(nuVisitor);
         //}

         ////new dp
         //if (newPressure)
         //{
         //   Grid3D::Interactor3DSet interactors = grid->getInteractors();
         //   interactors[0]->setGrid3D(grid);
         //   boost::dynamic_pointer_cast<D3Q27Interactor>(interactors[0])->deleteBCAdapter();
         //   BCAdapterPtr denBCAdapterFront(new DensityBCAdapter(rhoLBinflow));
         //   denBCAdapterFront->setBcAlgorithm(BCAlgorithmPtr(new EqDensityBCAlgorithm()));
         //   boost::dynamic_pointer_cast<D3Q27Interactor>(interactors[0])->addBCAdapter(denBCAdapterFront);
         //   interactors[0]->updateInteractor();
         //}

         if (myid==0)
         {
            UBLOG(logINFO, "Parameters:");         UBLOG(logINFO, "PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem());
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed());
            UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe());
            UBLOG(logINFO, "rho_LB = "<<rho_LB);
            UBLOG(logINFO, "nu_LB = "<<nu_LB);
            UBLOG(logINFO, "dp_LB = "<<dp_LB);
            UBLOG(logINFO, "dx = "<<deltax<<" m");
         }

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //BS visitor
         grid->accept(bcVisitor);

         if (myid==0) UBLOG(logINFO, "Restart - end");
      }



      SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      //nupsSch->addSchedule(nupsStep[0], nupsStep[1], nupsStep[2]);
      SPtr<NUPSCounterCoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> stepSch(new UbScheduler(outTimeStep, outTimeStart));

      SPtr<WriteMacroscopicQuantitiesCoProcessor> pp(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      deltax = grid->getDeltaX(baseLevel);
      double dxd2 = deltax/2.0;

      SPtr<IntegrateValuesHelper> ih1(new IntegrateValuesHelper(grid, comm, coord[0]-dxd2*10.0, coord[1]-dxd2, coord[2]-dxd2,
         coord[0]-dxd2*10.0-2.0*dxd2, coord[4]+dxd2, coord[5]+dxd2));

      SPtr<IntegrateValuesHelper> ih2(new IntegrateValuesHelper(grid, comm, coord[0], coord[1], coord[2], coord[3], coord[4], coord[5]));

      SPtr<IntegrateValuesHelper> ih3(new IntegrateValuesHelper(grid, comm, coord[3]+dxd2*10.0, coord[1]-dxd2, coord[2]-dxd2,
         coord[3]+dxd2*10.0+2.0*dxd2, coord[4]+dxd2, coord[5]+dxd2));

      if (myid==0) GbSystem3D::writeGeoObject(ih1->getBoundingBox().get(), pathname+"/geo/ih1", WbWriterVtkXmlBinary::getInstance());
      if (myid==0) GbSystem3D::writeGeoObject(ih2->getBoundingBox().get(), pathname+"/geo/ih2", WbWriterVtkXmlBinary::getInstance());
      if (myid==0) GbSystem3D::writeGeoObject(ih3->getBoundingBox().get(), pathname+"/geo/ih3", WbWriterVtkXmlBinary::getInstance());

      double factorp = 1; // dp_real / dp_LB;
      double factorv = 1;// dx / dt;
      SPtr<UbScheduler> stepMV(new UbScheduler(timeSeriesOutTime));

      SPtr<TimeseriesCoProcessor> tsp1(new TimeseriesCoProcessor(grid, stepMV, ih1, pathname+timeSeriesFile+"_1", comm));
      SPtr<TimeseriesCoProcessor> tsp2(new TimeseriesCoProcessor(grid, stepMV, ih2, pathname+timeSeriesFile+"_2", comm));
      SPtr<TimeseriesCoProcessor> tsp3(new TimeseriesCoProcessor(grid, stepMV, ih3, pathname+timeSeriesFile+"_3", comm));

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe());
      }

      omp_set_num_threads(numOfThreads);
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepSch, endTime));
      calculator->addCoProcessor(npr);
      calculator->addCoProcessor(restartCoProcessor);
      calculator->addCoProcessor(pp);
      calculator->addCoProcessor(tsp1);
      calculator->addCoProcessor(tsp2);
      calculator->addCoProcessor(tsp3);
      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      if (myid==0) UBLOG(logINFO, "Simulation-end");
   }
   catch (exception& e)
   {
      cerr<<e.what()<<endl<<flush;
   }
   catch (string& s)
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

   if (argv!=NULL)
   {
      if (argv[1]!=NULL)
      {
         run(string(argv[1]));
      }
      else
      {
         cout<<"Configuration file is missing!"<<endl;
      }
   }

   return 0;
}

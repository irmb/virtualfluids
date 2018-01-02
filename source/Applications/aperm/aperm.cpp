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
      int             numOfThreads = config.getInt("numOfThreads");
      string          sampleFilename = config.getString("sampleFilename");
      int             pmNX1 = config.getInt("pmNX1");
      int             pmNX2 = config.getInt("pmNX2");
      int             pmNX3 = config.getInt("pmNX3");
      double          lthreshold = config.getDouble("lthreshold");
      double          uthreshold = config.getDouble("uthreshold");
      double          pmL1 = config.getDouble("pmL1");
      double          pmL2 = config.getDouble("pmL2");
      double          pmL3 = config.getDouble("pmL3");
      int             blocknx = config.getInt("blocknx");
      //double          nx3 = config.getDouble("nx3");
      double          dpLB = config.getDouble("dp_LB");
      double          nu_LB = config.getDouble("nu_LB");
      string          timeSeriesFile = config.getString("timeSeriesFile");
      double          restartStep = config.getDouble("restartStep");
      //double          restartStepStart = config.getDouble("restartStepStart");
      double          endTime = config.getDouble("endTime");
      double          outTimeStep = config.getValue<double>("outTimeStep");
      double          outTimeStart = config.getValue<double>("outTimeStart");
      double          availMem = config.getDouble("availMem");
      bool            rawFile = config.getBool("rawFile");
      double          timeSeriesOutTime = config.getDouble("timeSeriesOutTime");
      bool            logToFile = config.getBool("logToFile");
      bool            spongeLayer = config.getBool("spongeLayer");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");
      double          deltax = config.getDouble("deltax");
      bool            newViscosity = config.getBool("newViscosity");
      bool            newPressure = config.getBool("newPressure");
      bool            pmDeltas = config.getBool("pmDeltas");
      double          pmDeltaX1 = config.getDouble("pmDeltaX1");
      double          pmDeltaX2 = config.getDouble("pmDeltaX2");
      double          pmDeltaX3 = config.getDouble("pmDeltaX3");
      double          vx1 = config.getDouble("vx1");
      double          vx2 = config.getDouble("vx2");
      double          vx3 = config.getDouble("vx3");
      bool            yDir = config.getBool("yDir");
      bool            zDir = config.getBool("zDir");
      double          cpStep = config.getDouble("cpStep");
      double          cpStepStart = config.getDouble("cpStepStart");
      bool            newStart = config.getValue<bool>("newStart");

      CommunicatorPtr comm = MPICommunicator::getInstance();
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

      //Sleep(30000);

      if (myid==0) UBLOG(logINFO, "Testcase permeability");

      string machinename = UbSystem::getMachineName();
      //UBLOG(logINFO, "PID = " << myid << " Hostname: " << machinename);
      //UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
      //UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
      //UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());

      int blocknx1 = blocknx;
      int blocknx2 = blocknx;
      int blocknx3 = blocknx;

      LBMReal rhoLB = 0.0;
      double rhoLBinflow = dpLB*3.0;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;

      double coord[6];
      //double deltax;

      Grid3DPtr grid(new Grid3D(comm));

      //////////////////////////////////////////////////////////////////////////
      //restart
      //UbSchedulerPtr rSch(new UbScheduler(cpStep, cpStepStart));
      //RestartCoProcessor rp(grid, rSch, comm, pathname, RestartCoProcessor::TXT);
      
      UbSchedulerPtr rSch2(new UbScheduler(cpStep, cpStepStart));
      MPIIORestart11CoProcessor rcp(grid, rSch2, pathname, comm);

      LBMKernelPtr kernel;
      kernel = LBMKernelPtr(new IncompressibleCumulantLBMKernel(blocknx1, blocknx2, blocknx3, IncompressibleCumulantLBMKernel::NORMAL));

      //BCProcessorPtr bcProc(new BCProcessor());
      BCProcessorPtr bcProc = BCProcessorPtr(new ThinWallBCProcessor());
      kernel->setBCProcessor(bcProc);

      rcp.setLBMKernel(kernel);
      rcp.setBCProcessor(bcProc);
      rcp.setChunk(1);
      //////////////////////////////////////////////////////////////////////////

      //BC Adapter
      //////////////////////////////////////////////////////////////////////////////
      BCAdapterPtr noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new ThinWallNoSlipBCAlgorithm()));
      //noSlipBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithm()));

      BCAdapterPtr denBCAdapterInflow(new DensityBCAdapter(rhoLBinflow));
      denBCAdapterInflow->setBcAlgorithm(BCAlgorithmPtr(new NonEqDensityBCAlgorithm()));

      BCAdapterPtr denBCAdapterOutflow(new DensityBCAdapter(rhoLB));
      denBCAdapterOutflow->setBcAlgorithm(BCAlgorithmPtr(new NonEqDensityBCAlgorithm()));
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

         double deltaVoxelX1 = pmL1/(double)pmNX1;
         double deltaVoxelX2 = pmL2/(double)pmNX2;
         double deltaVoxelX3 = pmL3/(double)pmNX3;

         if (pmDeltas)
         {
            deltaVoxelX1 = pmDeltaX1;
            deltaVoxelX2 = pmDeltaX2;
            deltaVoxelX3 = pmDeltaX3;
         }

         if (myid==0) UBLOG(logINFO, "read voxel matrix: start");
         GbVoxelMatrix3DPtr sample(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0, lthreshold, uthreshold));
         if (rawFile)
         {
            sample->readMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
         }
         else
         {
            sample->readMatrixFromVtiASCIIFile(samplePathname);
         }

         sample->setVoxelMatrixDelta((float)deltaVoxelX1, (float)deltaVoxelX2, (float)deltaVoxelX3);
         sample->setVoxelMatrixMininum(0.0, 0.0, 0.0);
         if (myid==0) UBLOG(logINFO, "read voxel matrix: end");


         if (myid==0) UBLOG(logINFO, "rotate voxel matrix: start");
         if (yDir)
         {
            sample->rotate90aroundZ();
            //sample->rotate90aroundZ();
            //sample->rotate90aroundZ();
         }
         if (zDir)
         {
            sample->rotate90aroundY();
         }
         if (myid==0) UBLOG(logINFO, "rotate voxel matrix: end");

         if (myid==0) sample->writeToVTKImageDataASCII(pathname+"/geo/sample");
        
         ///////////////////////////////////////////////////////

         ////////////////////////////////////////////////////////////////////////

         double offset1 = sample->getLengthX1()/10.0;
         double offset2 = 2.0*offset1;
         //double offset2 = offset1;
         //bounding box
         double g_minX1 = sample->getX1Minimum()-offset1;
         double g_minX2 = sample->getX2Minimum();
         double g_minX3 = sample->getX3Minimum();

         double g_maxX1 = sample->getX1Maximum()+offset2;
         double g_maxX2 = sample->getX2Maximum();
         double g_maxX3 = sample->getX3Maximum();

         ////////////////////////////////////////////////////////////////////////////
         //double nx1_temp = floor((g_maxX1-g_minX1)/(deltax*(double)blocknx));

         //deltax = (g_maxX1-g_minX1)/(nx1_temp*(double)blocknx);

         // g_maxX3 -= 0.5* deltax;
          ////////////////////////////////////////////////////////////////////////////

          ////deltax = (g_maxX3-g_minX3) /(nx3*blocknx3);

         double blockLength = (double)blocknx1*deltax;

         grid->setPeriodicX1(false);
         grid->setPeriodicX2(false);
         grid->setPeriodicX3(false);
         grid->setDeltaX(deltax);
         grid->setBlockNX(blocknx1, blocknx2, blocknx3);

         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid==0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         if (myid==0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "rho_LB = "<<rhoLB);
            UBLOG(logINFO, "nu_LB = "<<nu_LB);
            UBLOG(logINFO, "dp_LB = "<<dpLB);
            UBLOG(logINFO, "dx = "<<deltax<<" m");
            UBLOG(logINFO, "numOfThreads = "<<numOfThreads);
            UBLOG(logINFO, "path = "<<pathname);
            UBLOG(logINFO, "Preprozess - start");
         }

         //walls
         GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_minX2, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
         if (myid==0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1-blockLength, g_maxX2, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());


         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid==0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         WriteBlocksCoProcessorPtr ppblocks(new WriteBlocksCoProcessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         
         //PM interactor
         D3Q27InteractorPtr sampleInt(new D3Q27Interactor(sample, grid, noSlipBCAdapter, Interactor3D::SOLID));

         //wall interactors
         D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));


         D3Q27InteractorPtr inflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoInflow, grid, denBCAdapterInflow, Interactor3D::SOLID));

         //outflow
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, denBCAdapterOutflow, Interactor3D::SOLID));


         //UBLOG(logINFO, "PID = "<<myid<<" Hostname: "<<machinename);
         //UBLOG(logINFO, "PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem());
         //UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed());
         //UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe());


         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::RECURSIVE));
         ////////////////////////////////////////////
         //Zoltan
         //Grid3DVisitorPtr zoltanVisitor(new ZoltanPartitioningGridVisitor(comm, D3Q27System::BSW, 1));
         //grid->accept(zoltanVisitor);

         /////delete solid blocks
         if (myid==0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.addInteractor(sampleInt);
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

         //LBMKernelPtr kernel;
         //kernel = LBMKernelPtr(new IncompressibleCumulantLBMKernel(blocknx1, blocknx2, blocknx3, IncompressibleCumulantLBMKernel::NORMAL));

         ////BCProcessorPtr bcProc(new BCProcessor());
         //BCProcessorPtr bcProc = BCProcessorPtr(new ThinWallBCProcessor());
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
         fct.DefineConst("dp", dpLB);
         fct.DefineConst("x1max", g_maxX1);
         fct.DefineConst("l", g_maxX1-g_minX1);

         InitDistributionsBlockVisitor initVisitor(nu_LB, rhoLB);
         initVisitor.setRho(fct);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         WriteBoundaryConditionsCoProcessorPtr ppgeo(
            new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         ppgeo->process(0);
         ppgeo.reset();

         coord[0] = sample->getX1Minimum();
         coord[1] = sample->getX2Minimum();
         coord[2] = sample->getX3Minimum();
         coord[3] = sample->getX1Maximum();
         coord[4] = sample->getX2Maximum();
         coord[5] = sample->getX3Maximum();

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
         
         rcp.restart((int)restartStep);
         grid->setTimeStep(restartStep);

         //new nu
         if (newViscosity)
         {
            ViscosityBlockVisitor nuVisitor(nu_LB);
            grid->accept(nuVisitor);
         }

         //new dp
         if (newPressure)
         {
            Grid3D::Interactor3DSet interactors = grid->getInteractors();
            interactors[0]->setGrid3D(grid);
            boost::dynamic_pointer_cast<D3Q27Interactor>(interactors[0])->deleteBCAdapter();
            BCAdapterPtr denBCAdapterFront(new DensityBCAdapter(rhoLBinflow));
            denBCAdapterFront->setBcAlgorithm(BCAlgorithmPtr(new EqDensityBCAlgorithm()));
            boost::dynamic_pointer_cast<D3Q27Interactor>(interactors[0])->addBCAdapter(denBCAdapterFront);
            interactors[0]->updateInteractor();
         }

         if (myid==0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "rho_LB = "<<rhoLB);
            UBLOG(logINFO, "nu_LB = "<<nu_LB);
            UBLOG(logINFO, "dp_LB = "<<dpLB);
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

         UbSchedulerPtr geoSch(new UbScheduler(1));
         WriteBoundaryConditionsCoProcessor ppgeo = WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm);
         ppgeo.process(1);


         if (myid==0) UBLOG(logINFO, "Restart - end");
      }
      UbSchedulerPtr nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      //nupsSch->addSchedule(nupsStep[0], nupsStep[1], nupsStep[2]);
      NUPSCounterCoProcessor npr(grid, nupsSch, numOfThreads, comm);

      UbSchedulerPtr stepSch(new UbScheduler(outTimeStep,outTimeStart));

      WriteMacroscopicQuantitiesCoProcessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm);

      deltax = grid->getDeltaX(baseLevel);
      double dxd2 = deltax/2.0;

      IntegrateValuesHelperPtr ih1(new IntegrateValuesHelper(grid, comm, coord[0]-dxd2*10.0, coord[1]-dxd2, coord[2]-dxd2,
         coord[0]-dxd2*10.0-2.0*dxd2, coord[4]+dxd2, coord[5]+dxd2));

      //D3Q27IntegrateValuesHelperPtr ih2(new D3Q27IntegrateValuesHelper(grid, comm, coord[3]/2.0, coord[1] - dxd2, coord[2] - dxd2,
      //   coord[3]/2.0 + 2.0*dxd2, coord[4] + dxd2, coord[5] + dxd2));
      IntegrateValuesHelperPtr ih2(new IntegrateValuesHelper(grid, comm, coord[0], coord[1], coord[2], coord[3], coord[4], coord[5]));

      IntegrateValuesHelperPtr ih3(new IntegrateValuesHelper(grid, comm, coord[3]+dxd2*10.0, coord[1]-dxd2, coord[2]-dxd2,
         coord[3]+dxd2*10.0+2.0*dxd2, coord[4]+dxd2, coord[5]+dxd2));

      //D3Q27IntegrateValuesHelperPtr ih1(new D3Q27IntegrateValuesHelper(grid, comm, coord[0], coord[1], coord[2], coord[3], coord[4], coord[5]));
      if (myid==0) GbSystem3D::writeGeoObject(ih1->getBoundingBox().get(), pathname+"/geo/ih1", WbWriterVtkXmlBinary::getInstance());
      if (myid==0) GbSystem3D::writeGeoObject(ih2->getBoundingBox().get(), pathname+"/geo/ih2", WbWriterVtkXmlBinary::getInstance());
      if (myid==0) GbSystem3D::writeGeoObject(ih3->getBoundingBox().get(), pathname+"/geo/ih3", WbWriterVtkXmlBinary::getInstance());

      double factorp = 1; // dp_real / dp_LB;
      double factorv = 1;// dx / dt;
      UbSchedulerPtr stepMV(new UbScheduler(timeSeriesOutTime));

      TimeseriesCoProcessor tsp1(grid, stepMV, ih1, pathname+timeSeriesFile+"_1", comm);
      TimeseriesCoProcessor tsp2(grid, stepMV, ih2, pathname+timeSeriesFile+"_2", comm);
      TimeseriesCoProcessor tsp3(grid, stepMV, ih3, pathname+timeSeriesFile+"_3", comm);

      if (myid==0)
      {
         UBLOG(logINFO, "PID = "<<myid<<" Total Physical Memory (RAM): "<<Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used: "<<Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = "<<myid<<" Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe());
      }

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
      if (myid==0) UBLOG(logINFO, "Simulation-end");
      
      //////MPIIORestart2CoProcessor 
      //grid->deleteBlockIDs();
      //RenumberBlockVisitor renumber;
      //grid->accept(renumber);
      //UbSchedulerPtr iiSch(new UbScheduler(1));
      //MPIIORestart2CoProcessor rcpInit(grid, iiSch, pathname, comm);
      //rcpInit.process(300);
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

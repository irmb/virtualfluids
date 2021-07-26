#include <iostream>
#include <string>
#include <VirtualFluids.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////
void perm(string configname)
{
   try
   {
      vf::basics::ConfigurationFile   config;
      config.load(configname);

      string          pathname = config.getValue<string>("pathname");
      string          pathGeo = config.getValue<string>("pathGeo");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      string          sampleFilename = config.getValue<string>("sampleFilename");
      int             pmNX1 = config.getValue<int>("pmNX1");
      int             pmNX2 = config.getValue<int>("pmNX2");
      int             pmNX3 = config.getValue<int>("pmNX3");
      double          lthreshold = config.getValue<double>("lthreshold");
      double          uthreshold = config.getValue<double>("uthreshold");
      double          pmL1 = config.getValue<double>("pmL1");
      double          pmL2 = config.getValue<double>("pmL2");
      double          pmL3 = config.getValue<double>("pmL3");
      int             blocknx = config.getValue<int>("blocknx");
      double          dpLB = config.getValue<double>("dpLB");
      double          nuLB = config.getValue<double>("nuLB");
      string          timeSeriesFile = config.getValue<string>("timeSeriesFile");
      double          restartStep = config.getValue<double>("restartStep");
      double          restartStepStart = config.getValue<double>("restartStepStart");
      double          endTime = config.getValue<double>("endTime");
      double          availMem = config.getValue<double>("availMem");
      bool            rawFile = config.getValue<bool>("rawFile");
      double          timeSeriesOutTime = config.getValue<double>("timeSeriesOutTime");
      bool            logToFile = config.getValue<bool>("logToFile");
      bool            newStart = config.getValue<bool>("newStart");
      double          cpStart = config.getValue<double>("cpStart");
      double          cpStep = config.getValue<double>("cpStep");
      vector<double>  nupsStep = config.getVector<double>("nupsStep");
      double          outTimeStep = config.getValue<double>("outTimeStep");
      double          outTimeStart = config.getValue<double>("outTimeStart");
      double          deltax = config.getValue<double>("deltax");
      bool            writeSampleToFile = config.getValue<bool>("writeSampleToFile");

      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (logToFile)
      {
#if defined(__unix__)
         if (myid == 0)
         {
            const char* str = pathname.c_str();
            int status = mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
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

      if (myid == 0) UBLOG(logINFO, "Testcase permeability");

      if (myid == 0)
      {
         //string machinename = UbSystem::getMachineName();
         //UBLOG(logINFO, "PID = " << myid << " Hostname: " << machinename);
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      int blocknx1 = blocknx;
      int blocknx2 = blocknx;
      int blocknx3 = blocknx;

      LBMReal rhoLB = 0.0;

      double rhoLBinflow = dpLB*3.0;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

      double coord[6];



      ///close void space
      //////////////////////////////////////////////////////////////////////////
      //{
      //   string samplePathname = pathGeo + sampleFilename;

      //   double deltaVoxelX1 = pmL1/(double)pmNX1;
      //   double deltaVoxelX2 = pmL2/(double)pmNX2;
      //   double deltaVoxelX3 = pmL3/(double)pmNX3;

      //   GbVoxelMatrix3DPtr sample(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0, lthreshold, uthreshold));
      //   if (rawFile)
      //   {
      //      sample->readMatrixFromRawFile<unsigned short>(samplePathname, GbVoxelMatrix3D::BigEndian);
      //   }
      //   else
      //   {
      //      sample->readMatrixFromVtiASCIIFile(samplePathname);
      //   }

      //   sample->setVoxelMatrixDelta((float)deltaVoxelX1, (float)deltaVoxelX2, (float)deltaVoxelX3);
      //   sample->setVoxelMatrixMininum(0.0, 0.0, 0.0);

      //   if (myid == 0) sample->writeToVTKImageDataASCII(pathname + "/geo/sampleOpen");
      //   sample->calculateNumberOfSolidAndFluid();
      //   if (myid == 0)  UBLOG(logINFO, "number of solid = "<<sample->getNumberOfSolid());
      //   if (myid == 0)  UBLOG(logINFO, "number of fluid = "<<sample->getNumberOfFluid());

      //   sample->setClosedVoidSpaceToSolid();

      //   if (myid == 0) sample->writeToVTKImageDataASCII(pathname + "/geo/sampleClosed");

      //   sample->calculateNumberOfSolidAndFluid();
      //   if (myid == 0)  UBLOG(logINFO, "number of solid = "<<sample->getNumberOfSolid());
      //   if (myid == 0)  UBLOG(logINFO, "number of fluid = "<<sample->getNumberOfFluid());

      //   UBLOG(logINFO, "Finish!");
      //   return;
      //}
      //////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////
      //Grid
      //////////////////////////////////////////////////////////////////////////
      SPtr<Grid3D> grid(new Grid3D(comm));

      //BC adapters
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapterInflow(new DensityBCAdapter(rhoLBinflow));
      denBCAdapterInflow->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapterOutflow(new DensityBCAdapter(rhoLB));
      denBCAdapterOutflow->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

      //////////////////////////////////////////////////////////////////////////////////
      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(denBCAdapterInflow);
      bcVisitor.addBC(denBCAdapterOutflow);;

      SPtr<BCProcessor> bcProc;
      bcProc = SPtr<BCProcessor>(new BCProcessor());

      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());

      kernel->setBCProcessor(bcProc);


      //////////////////////////////////////////////////////////////////////////
      //restart
      SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
      SPtr<MPIIOMigrationCoProcessor> migCoProcessor(new MPIIOMigrationCoProcessor(grid, mSch, pathname+"/mig", comm));
      migCoProcessor->setLBMKernel(kernel);
      migCoProcessor->setBCProcessor(bcProc);
      //////////////////////////////////////////////////////////////////////////
      
	 if (myid == 0)
	 {
		UBLOG(logINFO, "Parameters:");
		UBLOG(logINFO, "rhoLB = " << rhoLB);
		UBLOG(logINFO, "nuLB = " << nuLB);
		UBLOG(logINFO, "dpLB = " << dpLB);
		UBLOG(logINFO, "dx = " << deltax << " m");

		UBLOG(logINFO, "numOfThreads = " << numOfThreads);
		UBLOG(logINFO, "path = " << pathname);
		UBLOG(logINFO, "Preprozess - start");
	 }      

      if (newStart)
      {
         if (myid == 0) UBLOG(logINFO, "new start..");

         if (myid == 0)
         {
            //UBLOG(logINFO, "new start PID = " << myid << " Hostname: " << machinename);
            UBLOG(logINFO, "new start PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
            UBLOG(logINFO, "new start PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
            UBLOG(logINFO, "new start PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
         }

         string samplePathname = pathGeo + sampleFilename;

         double deltaVoxelX1 = pmL1/(double)pmNX1;
         double deltaVoxelX2 = pmL2/(double)pmNX2;
         double deltaVoxelX3 = pmL3/(double)pmNX3;

         SPtr<GbVoxelMatrix3D> sample(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0, lthreshold, uthreshold));
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

         if (myid == 0 && writeSampleToFile) sample->writeToVTKImageDataASCII(pathname + "/geo/sample");

         ///////////////////////////////////////////////////////

         ////////////////////////////////////////////////////////////////////////

         double offset1 = sample->getLengthX1()/10.0;
         double offset2 = 2.0*offset1;
         //bounding box
         double g_minX1 = sample->getX1Minimum() - offset1;
         double g_minX2 = sample->getX2Minimum();
         double g_minX3 = sample->getX3Minimum();

         double g_maxX1 = sample->getX1Maximum() + offset2;
         double g_maxX2 = sample->getX2Maximum();
         double g_maxX3 = sample->getX3Maximum();

         double blockLength = (double)blocknx1*deltax;

         grid->setPeriodicX1(false);
         grid->setPeriodicX2(false);
         grid->setPeriodicX3(false);
         grid->setDeltaX(deltax);
         grid->setBlockNX(blocknx1, blocknx2, blocknx3);

         SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         //walls
         GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_minX2, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1-blockLength, g_maxX2, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());


         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         //PM interactor
         SPtr<D3Q27Interactor> sampleInt(new D3Q27Interactor(sample, grid, noSlipBCAdapter, Interactor3D::SOLID));

         //wall interactors
         SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));

		 //inflow
         SPtr<D3Q27Interactor> inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, denBCAdapterInflow, Interactor3D::SOLID));

         //outflow
         SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, denBCAdapterOutflow, Interactor3D::SOLID));;

         if (myid == 0)
         {
            //UBLOG(logINFO, "PID = " << myid << " Hostname: " << machinename);
            UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
            UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
            UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
         }

         ////////////////////////////////////////////
         //METIS
          SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B, MetisPartitioner::RECURSIVE));
         ////////////////////////////////////////////

         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.addInteractor(sampleInt);
         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         {
            WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm);
            ppblocks.process(1);
         }

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx1)* (blocknx2)* (blocknx3);
         unsigned long nod = nob * (blocknx1)* (blocknx2)* (blocknx3);
         unsigned long nodg = nob * (blocknx1 + gl) * (blocknx2 + gl) * (blocknx3 + gl);
         double needMemAll = double(nodg*(27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
         double needMem = needMemAll / double(comm->getNumberOfProcesses());

         if (myid == 0)
         {
            UBLOG(logINFO, "Number of blocks = " << nob);
            UBLOG(logINFO, "Number of nodes  = " << nod);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for (int level = minInitLevel; level <= maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO, "Number of blocks for level " << level << " = " << nobl);
               UBLOG(logINFO, "Number of nodes for level " << level << " = " << nobl*nodb);
            }
            UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
            UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
            UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
         }


         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);


         //BC
         intHelper.setBC();


         //Press*1.6e8+(14.76-coordsX)/3.5*5000
         //initialization of distributions
         mu::Parser fct;
         fct.SetExpr("(x1max-x1)/l*dp*3.0");
         fct.DefineConst("dp", dpLB);
         fct.DefineConst("x1max", g_maxX1);
         fct.DefineConst("l", g_maxX1-g_minX1);

         InitDistributionsBlockVisitor initVisitor;
         initVisitor.setRho(fct);
         grid->accept(initVisitor);

         //Post process
         {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm);
            ppgeo.process(0);
         }

         coord[0] = sample->getX1Minimum();
         coord[1] = sample->getX2Minimum();
         coord[2] = sample->getX3Minimum();
         coord[3] = sample->getX1Maximum();
         coord[4] = sample->getX2Maximum();
         coord[5] = sample->getX3Maximum();

         ////////////////////////////////////////////////////////
         FILE * pFile;
         string str = pathname + "/checkpoints/coord.txt";
         pFile = fopen(str.c_str(), "w");
         fprintf(pFile, "%g\n", deltax);
         fprintf(pFile, "%g\n", coord[0]);
         fprintf(pFile, "%g\n", coord[1]);
         fprintf(pFile, "%g\n", coord[2]);
         fprintf(pFile, "%g\n", coord[3]);
         fprintf(pFile, "%g\n", coord[4]);
         fprintf(pFile, "%g\n", coord[5]);
         fclose(pFile);
         ////////////////////////////////////////////////////////

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         ////////////////////////////////////////////////////////
         FILE * pFile;
         string str = pathname + "/checkpoints/coord.txt";
         pFile = fopen(str.c_str(), "r");
         fscanf(pFile, "%lg\n", &deltax);
         fscanf(pFile, "%lg\n", &coord[0]);
         fscanf(pFile, "%lg\n", &coord[1]);
         fscanf(pFile, "%lg\n", &coord[2]);
         fscanf(pFile, "%lg\n", &coord[3]);
         fscanf(pFile, "%lg\n", &coord[4]);
         fscanf(pFile, "%lg\n", &coord[5]);
         fclose(pFile);
         ////////////////////////////////////////////////////////

         migCoProcessor->restart((int)restartStep);
         grid->setTimeStep(restartStep);
         
         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      
      ////set connectors
      SPtr<InterpolationProcessor> iProcessor(new CompressibleOffsetMomentsInterpolationProcessor());
      SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
      grid->accept(setConnsVisitor);

      //bcVisitor should be accept after initialization!!!!
      grid->accept(bcVisitor);
      if (myid == 0) UBLOG(logINFO, "grid->accept(bcVisitor):end");
      
      SPtr<UbScheduler> nupsSch(new UbScheduler(nupsStep[0], nupsStep[1], nupsStep[2]));
      std::shared_ptr<CoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> stepSch(new UbScheduler(outTimeStep, outTimeStart));

      SPtr<CoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));

      deltax = grid->getDeltaX(baseLevel);
      double dxd2 = deltax / 2.0;

      SPtr<IntegrateValuesHelper> ih1(new IntegrateValuesHelper(grid, comm, coord[0] - dxd2*10.0, coord[1] - dxd2, coord[2] - dxd2,
         coord[0] - dxd2*10.0 - 2.0*dxd2, coord[4] + dxd2, coord[5] + dxd2));

      //D3Q27IntegrateValuesHelperPtr ih2(new D3Q27IntegrateValuesHelper(grid, comm, coord[3]/2.0, coord[1] - dxd2, coord[2] - dxd2,
      //   coord[3]/2.0 + 2.0*dxd2, coord[4] + dxd2, coord[5] + dxd2));
      SPtr<IntegrateValuesHelper> ih2(new IntegrateValuesHelper(grid, comm, coord[0], coord[1], coord[2], coord[3], coord[4], coord[5]));

      SPtr<IntegrateValuesHelper> ih3(new IntegrateValuesHelper(grid, comm, coord[3] + dxd2*10.0, coord[1] - dxd2, coord[2] - dxd2,
         coord[3] + dxd2*10.0 + 2.0*dxd2, coord[4] + dxd2, coord[5] + dxd2));

      //D3Q27IntegrateValuesHelperPtr ih1(new D3Q27IntegrateValuesHelper(grid, comm, coord[0], coord[1], coord[2], coord[3], coord[4], coord[5]));
      if (myid == 0) GbSystem3D::writeGeoObject(ih1->getBoundingBox().get(), pathname + "/geo/ih1", WbWriterVtkXmlBinary::getInstance());
      if (myid == 0) GbSystem3D::writeGeoObject(ih2->getBoundingBox().get(), pathname + "/geo/ih2", WbWriterVtkXmlBinary::getInstance());
      if (myid == 0) GbSystem3D::writeGeoObject(ih3->getBoundingBox().get(), pathname + "/geo/ih3", WbWriterVtkXmlBinary::getInstance());

      double factorp = 1; // dp_real / dpLB;
      double factorv = 1;// dx / dt;
      SPtr<UbScheduler> stepMV(new UbScheduler(timeSeriesOutTime));
      
      SPtr<CoProcessor> tsp1(new TimeseriesCoProcessor(grid, stepMV, ih1, pathname+timeSeriesFile+"_1", comm));
      SPtr<CoProcessor> tsp2(new TimeseriesCoProcessor(grid, stepMV, ih2, pathname+timeSeriesFile+"_2", comm));
      SPtr<CoProcessor> tsp3(new TimeseriesCoProcessor(grid, stepMV, ih3, pathname+timeSeriesFile+"_3", comm));

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      omp_set_num_threads(numOfThreads);
      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(tsp1);
      calculator->addCoProcessor(tsp2);
      calculator->addCoProcessor(tsp3);
      calculator->addCoProcessor(writeMQCoProcessor);
      calculator->addCoProcessor(migCoProcessor);
      


      if (myid==0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      if (myid==0) UBLOG(logINFO, "Simulation-end");
   }
   catch (exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch (string& s)
   {
      cerr << s << endl;
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
         perm(string(argv[1]));
      }
      else
      {
         cout<<"Configuration file must be set!: "<<argv[0]<<" <config file>"<<endl<<std::flush;
      }
   }

   return 0;
}

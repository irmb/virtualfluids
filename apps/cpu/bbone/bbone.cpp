#include <iostream>
#include <string>

#include <vfluids.h>

using namespace std;


void sbonepd(string configname)
{
   try
   {
      Configuration   config;
      config.load(configname);

      string          pathname          = config.getString("pathname");
      string          pathGeo           = config.getString("pathGeo");
      int             numOfThreads      = config.getInt("numOfThreads");
      string          sampleFilename    = config.getString("sampleFilename");
      int             pmNX1             = config.getInt("pmNX1");
      int             pmNX2             = config.getInt("pmNX2");
      int             pmNX3             = config.getInt("pmNX3");
      double          lthreshold        = config.getDouble("lthreshold");
      double          uthreshold        = config.getDouble("uthreshold");
      double          dp_real           = config.getDouble("dp_real");
      string          timeSeriesFile    = config.getString("timeSeriesFile");
      double          restartStep       = config.getDouble("restartStep");
      double          restartStepStart  = config.getDouble("restartStepStart");
      double          endTime           = config.getDouble("endTime");
      double          outTime           = config.getDouble("outTime");
      double          availMem          = config.getDouble("availMem");
      double          timeSeriesOutTime = config.getDouble("timeSeriesOutTime");
      bool            logToFile         = config.getBool("logToFile");
      double          deltaT            = config.getDouble("deltaT");
      
      CommunicatorPtr comm = vf::parallel::MPICommunicator::getInstance();
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



      if (myid == 0) UBLOG(logINFO, "Testcase big bone");

      Grid3DPtr grid(new Grid3D(comm));
      double deltaVoxel = 11.658e-6;

      double dx = deltaVoxel;

      const int blocknx1 = 64;
      const int blocknx2 = 64;
      const int blocknx3 = 64;

      LBMReal rho_LB = 0.0;
      //nueWasser = 1e-6 m^2/s
      double nu_real = 1e-6;
      LBMReal dt = deltaT; //1e-5; // s (frei gew�hlt)
      //dx - frei gew�hlt
      //
      LBMReal nu_LB = nu_real / (dx*dx / dt);


      //dp = 50000 Pa - 0 Pa = 50000 Pa
      //rho wasser = 1000 kg*m^-3
      double rho_real = 1000;
      //dp/rho = 50000/1000 = 50 m^2/s^2
      double dp_div_rho_real = dp_real / rho_real;

      double dp_LB = dp_div_rho_real / ((dx / dt)*(dx / dt));

      double rhoLBinflow;
      rhoLBinflow = dp_LB*3.0;
     
      double deltax = dx;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = 0;

      double coord[6];

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(restartStep, restartStepStart));
      RestartPostprocessor rp(grid, rSch, comm, pathname, RestartPostprocessor::BINARY);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep() == 0)
      {
         if (myid == 0) UBLOG(logINFO, "Neustart..");

         string boneFilename = pathGeo + sampleFilename;

         //int pmNX1 = 1164;  //abmessung einzelbild in x-richtung
         //int pmNX2 = 972; //abmessung einzelbild in y richtung
         //int pmNX3 = 900; //anzahl der bilder
         ////int pmNX3 = 10; //anzahl der bilder
         //float lthreshold = 109.0;
         //float uthreshold = 255.0;

         GbVoxelMatrix3DPtr bone(new GbVoxelMatrix3D(pmNX1, pmNX2, pmNX3, 0, lthreshold, uthreshold));
         bone->readMatrixFromRawFile<unsigned char>(boneFilename, GbVoxelMatrix3D::BigEndian);
         bone->setVoxelMatrixDelta(deltaVoxel, deltaVoxel, deltaVoxel);
         bone->setVoxelMatrixMininum(0.0, 0.0, 0.0);

         if (myid == 0) bone->writeToVTKImageDataASCII(pathname + "/geo/bone");

         ///////////////////////////////////////////////////////

         ////////////////////////////////////////////////////////////////////////

         double offset = 0.5e-3;
         //bounding box
         double g_minX1 = bone->getX1Minimum();
         double g_minX2 = bone->getX2Minimum();
         double g_minX3 = bone->getX3Minimum() - offset;

         double g_maxX1 = bone->getX1Maximum();
         double g_maxX2 = bone->getX2Maximum();
         double g_maxX3 = bone->getX3Maximum() + offset;

         double blockLength = (double)blocknx1*deltax;

         grid->setPeriodicX1(false);
         grid->setPeriodicX2(false);
         grid->setPeriodicX3(false);
         grid->setDeltaX(deltax);
         grid->setBlockNX(blocknx1, blocknx2, blocknx3);

         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "rho_LB = " << rho_LB);
            UBLOG(logINFO, "nu_LB = " << nu_LB);
            UBLOG(logINFO, "dp_LB = " << dp_LB);
            UBLOG(logINFO, "dx = " << dx << " m");
            UBLOG(logINFO, "dt = " << dt << " s");
            UBLOG(logINFO, "rho_real = " << rho_real << " kg*m^-3");
            UBLOG(logINFO, "nu_real = " << nu_real << " m^2/s");
            UBLOG(logINFO, "dp_real = " << dp_real << " Pa");

            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
            UBLOG(logINFO, "Preprozess - start");
         }

         //cylinder
         double radius = 0.0036;
         double cx1 = 0.007;
         double cx2 = 0.0046;

         GbObject3DPtr cylinder(new GbCylinder3D(cx1, cx2, g_minX3 - offset, cx1, cx2, g_maxX3 + offset, radius));
         GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength, g_maxX2 + blockLength, g_minX3));
         if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_maxX3, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         //bone interactor
         int bcOptionBone = 1; //0=simple Bounce Back, 1=quadr. BB, 2=thin wall
         D3Q27BoundaryConditionAdapterPtr bcBone(new D3Q27NoSlipBCAdapter(bcOptionBone));
         D3Q27InteractorPtr boneInt(new D3Q27Interactor(bone, grid, bcBone, Interactor3D::SOLID));

         //wall interactors
         int bcOptionWall = 1; //0=simple Bounce Back, 1=quadr. BB, 2=thin wall
         D3Q27BoundaryConditionAdapterPtr bcWall(new D3Q27NoSlipBCAdapter(bcOptionWall));
         D3Q27InteractorPtr cylInt(new D3Q27Interactor(cylinder, grid, bcWall, Interactor3D::INVERSESOLID));

         D3Q27BoundaryConditionAdapterPtr denBCAdapterInflow(new D3Q27DensityBCAdapter(rhoLBinflow));
         denBCAdapterInflow->setSecondaryBcOption(0);
         D3Q27InteractorPtr inflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoInflow, grid, denBCAdapterInflow, Interactor3D::SOLID));

         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapterOutflow(new D3Q27DensityBCAdapter(rho_LB));
         denBCAdapterOutflow->setSecondaryBcOption(0);
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoOutflow, grid, denBCAdapterOutflow, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW));
         ////////////////////////////////////////////
         //Zoltan
         //Grid3DVisitorPtr zoltanVisitor(new ZoltanPartitioningGridVisitor(comm, D3Q27System::BSW, 1));
         //////////////////////////////////////////////////////////////////////////
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(boneInt);
         intHelper.addInteractor(cylInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         ppblocks->update(0);
         ppblocks.reset();

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

         LBMKernel3DPtr kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));

         //BCProcessorPtr bcProc(new D3Q27ETForThinWallBCProcessor());
         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         BoundaryConditionPtr densityBC(new NonEqDensityBoundaryCondition());
         BoundaryConditionPtr noSlipBC(new NoSlipBoundaryCondition());
         bcProc->addBC(densityBC);
         bcProc->addBC(noSlipBC);
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nu_LB, availMem, needMem);
         grid->accept(kernelVisitor);

         //BC
         intHelper.setBC();

         BoundaryConditionBlockVisitor bcVisitor;
         grid->accept(bcVisitor);

         //Press*1.6e8+(14.76-coordsX)/3.5*5000
         //initialization of distributions
         mu::Parser fct;
         fct.SetExpr("(x3max-x3)/l*dp*3.0");
         fct.DefineConst("dp", dp_LB);
         fct.DefineConst("x3max", g_maxX3);
         fct.DefineConst("l", g_maxX3-g_minX3);

         D3Q27ETInitDistributionsBlockVisitor initVisitor(nu_LB, rho_LB);
         initVisitor.setRho(fct);
         //initVisitor.setVx1(fct);
         //initVisitor.setVx1(0.0);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
         ppgeo->update(0);
         ppgeo.reset();


         coord[0] = bone->getX1Minimum();
         coord[1] = bone->getX2Minimum();
         coord[2] = bone->getX3Minimum();//cylinder->getX3Centroid();
         coord[3] = bone->getX1Maximum();
         coord[4] = bone->getX2Maximum();
         coord[5] = bone->getX3Maximum(); //cylinder->getX3Centroid();

         ////////////////////////////////////////////////////////
         FILE * pFile;
         string str = pathname + "/checkpoints/coord.txt";
         pFile = fopen(str.c_str(), "w");
         fprintf(pFile, "%g\n", coord[0]);
         fprintf(pFile, "%g\n", coord[1]);
         fprintf(pFile, "%g\n", coord[2]);
         fprintf(pFile, "%g\n", coord[3]);
         fprintf(pFile, "%g\n", coord[4]);
         fprintf(pFile, "%g\n", coord[5]);
         fclose(pFile);
         ////////////////////////////////////////////////////////
         grid->addInteractor(inflowInt);
         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         Grid3D::Interactor3DSet interactors = grid->getInteractors();
         interactors[0]->setGrid3D(grid);
         boost::dynamic_pointer_cast<D3Q27Interactor>(interactors[0])->deleteBCAdapter();
         D3Q27BoundaryConditionAdapterPtr denBCAdapterFront(new D3Q27DensityBCAdapter(rhoLBinflow));
         boost::dynamic_pointer_cast<D3Q27Interactor>(interactors[0])->addBCAdapter(denBCAdapterFront);
         interactors[0]->updateInteractor();

         UBLOG(logINFO, "rhoLBinflow = "<<rhoLBinflow);

         BoundaryConditionBlockVisitor bcVisitor;
         grid->accept(bcVisitor);

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept(setConnsVisitor);

         ////////////////////////////////////////////////////////
         FILE * pFile;
         string str = pathname + "/checkpoints/coord.txt";
         pFile = fopen(str.c_str(), "r");
         fscanf(pFile, "%lg\n", &coord[0]);
         fscanf(pFile, "%lg\n", &coord[1]);
         fscanf(pFile, "%lg\n", &coord[2]);
         fscanf(pFile, "%lg\n", &coord[3]);
         fscanf(pFile, "%lg\n", &coord[4]);
         fscanf(pFile, "%lg\n", &coord[5]);
         fclose(pFile);
         ////////////////////////////////////////////////////////

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      //stepSch->addSchedule(10, 10, 10);
      //stepSch->addSchedule(100, 100, 100);
      //stepSch->addSchedule(1000, 1000, 1000);
      //stepSch->addSchedule(100, 1500, 2000);
      //stepSch->addSchedule(10000, 10000, 10000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv);

      double dxd2 = deltax / 2.0;
      D3Q27IntegrateValuesHelperPtr ih1(new D3Q27IntegrateValuesHelper(grid, comm, coord[0], coord[1], coord[2] - dxd2*10.0,
         coord[3], coord[4], coord[2] - dxd2*10.0 - 2.0*dxd2));
      D3Q27IntegrateValuesHelperPtr ih2(new D3Q27IntegrateValuesHelper(grid, comm, coord[0], coord[1], coord[2], coord[3], coord[4], coord[5]));
      D3Q27IntegrateValuesHelperPtr ih3(new D3Q27IntegrateValuesHelper(grid, comm, coord[0], coord[1], coord[5] + dxd2*10.0,
         coord[3], coord[4], coord[5] + dxd2*10.0 + 2.0*dxd2));

      if (myid == 0) GbSystem3D::writeGeoObject(ih1->getBoundingBox().get(), pathname + "/geo/ih1", WbWriterVtkXmlBinary::getInstance());
      if (myid == 0) GbSystem3D::writeGeoObject(ih2->getBoundingBox().get(), pathname + "/geo/ih2", WbWriterVtkXmlBinary::getInstance());
      if (myid == 0) GbSystem3D::writeGeoObject(ih3->getBoundingBox().get(), pathname + "/geo/ih3", WbWriterVtkXmlBinary::getInstance());

      UbSchedulerPtr stepMV(new UbScheduler(timeSeriesOutTime));

      TimeseriesPostprocessor tsp1(grid, stepMV, ih1, pathname+timeSeriesFile+"_1", comm);
      TimeseriesPostprocessor tsp2(grid, stepMV, ih2, pathname+timeSeriesFile+"_2", comm);
      TimeseriesPostprocessor tsp3(grid, stepMV, ih3, pathname+timeSeriesFile+"_3", comm);

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepMV));
      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculation->calculate();
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
         sbonepd(string(argv[1]));
   }

   return 0;
}

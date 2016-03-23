#include <iostream>
#include <string>

#include <vfluids.h>

using namespace std;


void sbonepd(const char *configname)
{
   try
   {

      string machine = QUOTEME(CAB_MACHINE);
      string pathname, pathGeo;
      int numOfThreads;
      double availMem;

      ConfigFileReader cf(configname);
      if (!cf.read())
      {
         std::string exceptionText = "Unable to read configuration file\n";
         throw exceptionText;
      }

      CommunicatorPtr comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if (machine == "BOMBADIL")
      {
         numOfThreads = 4;
         pathname = "d:/temp/bbone";
         pathGeo = "d:/Data/Bone/BigBone";
         availMem = 15.0e9;
      }
      else if (machine == "M01" || machine == "M02")
      {
         numOfThreads = 8;
         pathname = cf.getValue("pathname");
         pathGeo = cf.getValue("pathGeo");
         availMem = 12.0e9;

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
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");



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
      LBMReal dt = 5e-8; // s (frei gewählt)
      //dx - frei gewählt
      //
      LBMReal nu_LB = nu_real / (dx*dx / dt);


      //dp = 50000 Pa - 0 Pa = 50000 Pa
      double dp_real = UbSystem::stringTo<double>(cf.getValue("pressure")); //5000;
      //rho wasser = 1000 kg*m^-3
      double rho_real = 1000;
      //dp/rho = 50000/1000 = 50 m^2/s^2
      double dp_div_rho_real = dp_real / rho_real;

      double dp_LB = dp_div_rho_real / ((dx / dt)*(dx / dt));

      bool with_forcing = true;

      double rhoLBinflow;
      if (with_forcing)
      {
         rhoLBinflow = 0.0;
      }
      else
      {
         rhoLBinflow = dp_LB*3.0;
      }
      double deltax = dx;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = 0;

      double coord[6];

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr rSch(new UbScheduler(50000, 50000, 10000000));
      RestartPostprocessor rp(grid, rSch, comm, pathname, RestartPostprocessor::BINARY);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep() == 0)
      {
         if (myid == 0) UBLOG(logINFO, "Neustart..");

         string boneFilename = pathGeo + "/cyl_bone2.raw";

         int pmNX1 = 1164;  //abmessung einzelbild in x-richtung
         int pmNX2 = 972; //abmessung einzelbild in y richtung
         int pmNX3 = 900; //anzahl der bilder
         //int pmNX3 = 10; //anzahl der bilder
         float lthreshold = 109.0;
         float uthreshold = 255.0;

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
         grid->setPeriodicX3(true);
         grid->setDeltaX(deltax);
         grid->setBlockNX(blocknx1, blocknx2, blocknx3);

         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());


         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         double forcing = 0;
         if (with_forcing)
         {
            forcing = dp_LB / (blocknx3*grid->getNX3());
         }

         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "with forcing = " << with_forcing);
            UBLOG(logINFO, "rho_LB = " << rho_LB);
            UBLOG(logINFO, "nu_LB = " << nu_LB);
            UBLOG(logINFO, "dp_LB = " << dp_LB);
            UBLOG(logINFO, "forcing = " << forcing);
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
         int bcOptionBone = 0; //0=simple Bounce Back, 1=quadr. BB, 2=thin wall
         D3Q27BoundaryConditionAdapterPtr bcBone(new D3Q27NoSlipBCAdapter(bcOptionBone));
         D3Q27InteractorPtr boneInt(new D3Q27Interactor(bone, grid, bcBone, Interactor3D::SOLID));

         //wall interactors
         int bcOptionWall = 0; //0=simple Bounce Back, 1=quadr. BB, 2=thin wall
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

         mu::Parser fctForcingX3;
         fctForcingX3.SetExpr("Fx3");
         fctForcingX3.DefineConst("Fx3", forcing);

         kernel->setForcingX3(fctForcingX3);
         kernel->setWithForcing(true);

         //BCProcessorPtr bcProc(new D3Q27ETForThinWallBCProcessor());
         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nu_LB, availMem, needMem);
         grid->accept(kernelVisitor);


         //BC
         intHelper.setBC();

         //Press*1.6e8+(14.76-coordsX)/3.5*5000
         //initialization of distributions
         //mu::Parser fct;
         //fct.SetExpr("(x1max-x1)/l*dp*3.0");
         //fct.DefineConst("dp", dp_LB);
         //fct.DefineConst("x3max", g_maxX3);
         //fct.DefineConst("l", g_maxX3-g_minX3);

         D3Q27ETInitDistributionsBlockVisitor initVisitor(nu_LB, rho_LB);
         //initVisitor.setRho(fct);
         //initVisitor.setVx1(fct);
         initVisitor.setVx1(0.0);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
         ppgeo->update(0);
         ppgeo.reset();


         coord[0] = bone->getX1Minimum();
         coord[1] = bone->getX2Minimum();
         coord[2] = cylinder->getX3Centroid();
         coord[3] = bone->getX1Maximum();
         coord[4] = bone->getX2Maximum();
         coord[5] = cylinder->getX3Centroid();

         ////////////////////////////////////////////////////////
         FILE * pFile;
         string str = pathname + "/checkpoints/coord.txt";
         pFile = fopen(str.c_str(), "w");
         fprintf(pFile, "%f\n", coord[0]);
         fprintf(pFile, "%f\n", coord[1]);
         fprintf(pFile, "%f\n", coord[2]);
         fprintf(pFile, "%f\n", coord[3]);
         fprintf(pFile, "%f\n", coord[4]);
         fprintf(pFile, "%f\n", coord[5]);
         fclose(pFile);
         ////////////////////////////////////////////////////////

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
         grid->accept(setConnsVisitor);

         ////////////////////////////////////////////////////////
         FILE * pFile;
         string str = pathname + "/checkpoints/coord.txt";
         pFile = fopen(str.c_str(), "r");
         fscanf(pFile, "%f\n", &coord[0]);
         fscanf(pFile, "%f\n", &coord[1]);
         fscanf(pFile, "%f\n", &coord[2]);
         fscanf(pFile, "%f\n", &coord[3]);
         fscanf(pFile, "%f\n", &coord[4]);
         fscanf(pFile, "%f\n", &coord[5]);
         fclose(pFile);
         ////////////////////////////////////////////////////////

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      double outTime = 30000;
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      stepSch->addSchedule(10, 10, 10);
      stepSch->addSchedule(100, 100, 100);
      stepSch->addSchedule(1000, 1000, 1000);
      stepSch->addSchedule(100, 1500, 2000);
      stepSch->addSchedule(10000, 10000, 10000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv);

      double dxd2 = deltax / 2.0;
      D3Q27IntegrateValuesHelperPtr ih1(new D3Q27IntegrateValuesHelper(grid, comm, coord[0] - dxd2, coord[1] - dxd2, coord[2] - dxd2,
         coord[3] + dxd2, coord[4] + dxd2, coord[5] + dxd2));
      if (myid == 0) GbSystem3D::writeGeoObject(ih1->getBoundingBox().get(), pathname + "/geo/ih1", WbWriterVtkXmlBinary::getInstance());

      double factorp = dp_real / dp_LB;
      double factorv = dx / dt;
      UbSchedulerPtr stepMV(new UbScheduler(100));
      D3Q27MeanValuesPostprocessor mvp1(grid, stepMV, pathname + "/mv/mv1.txt", comm, ih1, factorp, factorv);


      //D3Q27IntegrateValuesHelperPtr ih2(new D3Q27IntegrateValuesHelper(grid, comm, g_maxX1-2.0*deltax, g_minX2, g_minX3,
      //   g_maxX1 - deltax, g_maxX2, g_maxX3));
      //if (myid == 0) GbSystem3D::writeGeoObject(ih2->getBoundingBox().get(), pathname + "/geo/ih2", WbWriterVtkXmlBinary::getInstance());

      //D3Q27MeanValuesPostprocessor mvp2(grid, stepSch, pathname + "/mv/mv2.txt", comm, ih2, factorp, factorv);

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      double endTime = UbSystem::stringTo<double>(cf.getValue("endTime")); //100001;//10001.0;

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
      sbonepd(argv[1]);
   }

   return 0;
}

#include <iostream>
#include <string>
#include "VirtualFluids.h"

using namespace std;

void generateCubes(double top_porosity, double bottom_porosity, std::array<double,3> dimensions, std::array<double,3> cubes_in_direction, string pathOut, SPtr<Grid3D> grid, SPtr<BCAdapter> noSlipBCAdapter)
{
   double x = dimensions[0];
   double y = dimensions[1];
   double z = dimensions[2];
   double num_x_cubes = cubes_in_direction[0];
   double num_y_cubes = cubes_in_direction[1];
   double num_z_cubes = cubes_in_direction[2];

   double H = z / 2.;

   double outer_cube_side_length = 20;
   double inner_cube_side_length = 10;

   double dx = (x - outer_cube_side_length) / (num_x_cubes - 1);
   double dy = (y - outer_cube_side_length) / (num_y_cubes - 1);
   double dz = 20; //(z - outer_cube_side_length) / (num_z_cubes - 1);

   double porosity_step_z = (top_porosity - bottom_porosity) / (num_z_cubes - 1);
   double porosity = float(bottom_porosity);
   double cube_side_length = pow(1.-porosity, 1. / 3) * outer_cube_side_length;
   double initial_x = 0. + outer_cube_side_length / 2.;
   double initial_y = 0. + outer_cube_side_length / 2.;
   double initial_z = 0. + outer_cube_side_length / 2.;

   double current_z = initial_z;
   
   double cube_side_length_temp = 0; 

   SPtr<GbObject3D> cube(new GbCuboid3D(0, 0, 0, inner_cube_side_length, inner_cube_side_length, inner_cube_side_length));
   double scaleFactor = cube_side_length / inner_cube_side_length;
   cube->scale(scaleFactor, scaleFactor, scaleFactor);

   std::vector< SPtr<Interactor3D> > interactors;

   SPtr<D3Q27Interactor> cubeInt(new D3Q27Interactor(cube, grid, noSlipBCAdapter, Interactor3D::SOLID));

   int num = 0;
   for(int k = 0; k < num_z_cubes; k++)
   {
      double current_y = initial_y;
      for(int j = 0; j < num_y_cubes; j++)
      {
         double current_x = initial_x;
         for(int i = 0; i < num_x_cubes; i++)
         {
            cube->setCenterCoordinates(current_x, current_y, current_z);
            //GbSystem3D::writeGeoObject(cube.get(), pathOut + "/cubes/cube"+ UbSystem::toString(num), WbWriterVtkXmlBinary::getInstance());
            std::vector< std::shared_ptr<Block3D> > blockVector;
            UbTupleInt3 blockNX=grid->getBlockNX();
            SPtr<GbObject3D> geoObject(cubeInt->getGbObject3D());
            double ext = 0.0;
            std::array<double, 6> AABB ={ geoObject->getX1Minimum(),geoObject->getX2Minimum(),geoObject->getX3Minimum(),geoObject->getX1Maximum(),geoObject->getX2Maximum(),geoObject->getX3Maximum() };
            grid->getBlocksByCuboid(AABB[0]-(double)val<1>(blockNX)*ext, AABB[1]-(double)val<2>(blockNX)*ext, AABB[2]-(double)val<3>(blockNX)*ext, AABB[3]+(double)val<1>(blockNX)*ext, AABB[4]+(double)val<2>(blockNX)*ext, AABB[5]+(double)val<3>(blockNX)*ext, blockVector);
            cubeInt->getBcNodeIndicesMap();
            dynamic_pointer_cast<Interactor3D>(cubeInt)->removeBcBlocks();
            for (std::shared_ptr<Block3D> block : blockVector)
            {
               if (block->getKernel())
               {
                  cubeInt->setBCBlock(block);
               }
            }
            cubeInt->initInteractor();
            num ++;
            current_x += dx;
         }
         current_y += dy;
      }
      current_z += dz;
      porosity += porosity_step_z;
      cube_side_length_temp = cube_side_length;
      cube_side_length = pow(1.-porosity, 1. / 3) * outer_cube_side_length;
      scaleFactor = cube_side_length / cube_side_length_temp;
      cube->scale(scaleFactor, scaleFactor, scaleFactor);
   }
}

//////////////////////////////////////////////////////////////////////////
void run(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathOut           = config.getValue<string>("pathOut");
      string          pathGeo           = config.getValue<string>("pathGeo");
      int             numOfThreads      = config.getValue<int>("numOfThreads");
      vector<int>     blocknx           = config.getVector<int>("blocknx");
      double          u_LB              = config.getValue<double>("u_LB");
      double          restartStep       = config.getValue<double>("restartStep");
      double          cpStep            = config.getValue<double>("cpStep");
      double          cpStart           = config.getValue<double>("cpStart");
      double          endTime           = config.getValue<double>("endTime");
      double          outTime           = config.getValue<double>("outTime");
      double          availMem          = config.getValue<double>("availMem");
      bool            logToFile         = config.getValue<bool>("logToFile");
      double          deltaXfine        = config.getValue<double>("deltaXfine");
      int             refineLevel       = config.getValue<int>("refineLevel");
      double          Re                = config.getValue<double>("Re");
      double          timeAvStart       = config.getValue<double>("timeAvStart");
      double          timeAvStop        = config.getValue<double>("timeAvStop");
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
      double lReal = 0.008/2.0;//m
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

      //SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel());
      SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulant4thOrderViscosityLBMKernel());
      
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

      double channel_hight = (g_maxX3-g_minX3)/2.0;
      double channel_hight_LB = channel_hight/deltaXcoarse;
      double d_p = channel_hight/20.0;
      //////////////////////////////////////////////////////////////////////////
      double nu_LB = (u_LB*channel_hight_LB)/Re;
      //////////////////////////////////////////////////////////////////////////
      if (myid == 0)
      {
         UBLOG(logINFO, "Parameters:");
         UBLOG(logINFO, "Re                  = " << Re);
         UBLOG(logINFO, "u_LB                = " << u_LB);
         UBLOG(logINFO, "rho_LB              = " << rho_LB);
         UBLOG(logINFO, "nu_LB               = " << nu_LB);
         UBLOG(logINFO, "dx coarse           = " << deltaXcoarse);
         UBLOG(logINFO, "dx fine             = " << deltaXfine);
         UBLOG(logINFO, "channel_high        = " << channel_hight);
         UBLOG(logINFO, "channel_high_LB     = " << channel_hight_LB);
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
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::RECURSIVE));
         
         //DEBUG METIS 
         //////////////////////////////////////////////////////////////////////////
         
         SimpleGeometricPartitioner sgp;
         UbTupleInt3 dim = sgp.createDimensions(30, 20, 20, 20);


         dynamic_pointer_cast<MetisPartitioningGridVisitor>(metisVisitor)->setNumberOfProcesses(1500);
         grid->accept(metisVisitor);

         //SPtr<Grid3DVisitor> zoltanVisitor(new ZoltanPartitioningGridVisitor(comm, D3Q27System::BSW));
         //grid->accept(zoltanVisitor);

         WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
         ppblocks.process(0);
         return;
         //////////////////////////////////////////////////////////////////////////

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

         double bottom_porosity = 0.875;
         double top_porosity = 1 - (1 - bottom_porosity) / 9;
         generateCubes(top_porosity, bottom_porosity, std::array<double, 3>{600., 400., 400.}, std::array<double, 3>{30., 20., 9.}, pathOut, grid, noSlipBCAdapter);

         grid->accept(bcVisitor);

         mu::Parser inflowProfileVx1, inflowProfileVx2, inflowProfileVx3, inflowProfileRho;
         inflowProfileVx1.SetExpr("x3 < h ? 0.0 : uLB-1e-5*(x1+x2+x3)");
         inflowProfileVx1.DefineConst("uLB", u_LB);
         inflowProfileVx1.DefineConst("h", channel_hight-d_p);

         InitDistributionsBlockVisitor initVisitor;
         initVisitor.setVx1(inflowProfileVx1);
         //initVisitor.setVx1(u_LB);
         //initVisitor.setVx2(u_LB);
         //initVisitor.setVx3(u_LB);

         //initVisitor.setVx1(inflowProfileVx1);
         initVisitor.setVx1(0);

         grid->accept(initVisitor);

         ////set connectors
         //InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
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
         restartCoProcessor->restart((int)restartStep);
         //migCoProcessor->restart((int)restartStep);
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
      writeMQCoProcessor->process(0);

      SPtr<GbObject3D> bbBox(new GbCuboid3D(g_minX1-blockLength, (g_maxX2-g_minX2)/2.0, g_minX3-blockLength, g_maxX1+blockLength, (g_maxX2-g_minX2)/2.0+deltaXcoarse, g_maxX3+blockLength));
      if (myid==0) GbSystem3D::writeGeoObject(bbBox.get(), pathOut+"/geo/bbBox", WbWriterVtkXmlASCII::getInstance());
      SPtr<WriteMQFromSelectionCoProcessor> writeMQSelectCoProcessor(new WriteMQFromSelectionCoProcessor(grid, stepSch, bbBox, pathOut, WbWriterVtkXmlBinary::getInstance(), conv, comm));


      SPtr<UbScheduler> AdjForcSch(new UbScheduler());
      AdjForcSch->addSchedule(10, 0, 10000000);
      SPtr<IntegrateValuesHelper> intValHelp(new IntegrateValuesHelper(grid, comm, g_minX1, g_minX2, channel_hight, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(intValHelp->getBoundingBox().get(), pathOut + "/geo/IntValHelp", WbWriterVtkXmlBinary::getInstance());

      double vxTarget=u_LB;
      SPtr<AdjustForcingCoProcessor> AdjForcCoProcessor(new AdjustForcingCoProcessor(grid, AdjForcSch, pathOut, intValHelp, vxTarget, comm));


      std::vector<double> levelCoords;
      std::vector<int> levels;
      std::vector<double> bounds;

      bounds.push_back(g_minX1);
      bounds.push_back(g_minX2);
      bounds.push_back(g_minX3);
      bounds.push_back(g_maxX1);
      bounds.push_back(g_maxX2);
      bounds.push_back(g_maxX3);
      levels.push_back(0);
      levelCoords.push_back(g_minX3);
      levelCoords.push_back(g_maxX3);
      SPtr<UbScheduler> tavSch(new UbScheduler(1, timeAvStart, timeAvStop));
      SPtr<CoProcessor> timeAveragingCoProcessor(new TimeAveragedValuesCoProcessor(grid, pathOut, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,TimeAveragedValuesCoProcessor::Density |  TimeAveragedValuesCoProcessor::Velocity | TimeAveragedValuesCoProcessor::Fluctuations | TimeAveragedValuesCoProcessor::Triplecorrelations, levels, levelCoords, bounds));
      
      
      //create line time series
      //SPtr<UbScheduler> tpcSch(new UbScheduler(1,1,3));
      //GbPoint3DPtr p1(new GbPoint3D(0.0,0.005,0.01));
      //GbPoint3DPtr p2(new GbPoint3D(0.064,0.005,0.01));
      //SPtr<GbLine3D> line(new GbLine3D(p1.get(),p2.get()));
      //SPtr<GbLine3D> line(new GbLine3D(new GbPoint3D(0.0,0.005,0.01),new GbPoint3D(0.064,0.005,0.01)));
      //LineTimeSeriesCoProcessor lineTs(grid, tpcSch,pathOut+"/TimeSeries/line1.csv",line, 0,comm);
      //if (myid==0) lineTs.writeLine(pathOut+"/geo/line1");

      if (myid == 0)
      {
         UBLOG(logINFO, "PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO, "PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      omp_set_num_threads(numOfThreads);
      SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, (int)endTime));
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(AdjForcCoProcessor);
      calculator->addCoProcessor(migCoProcessor);
      //calculator->addCoProcessor(restartCoProcessor);
      calculator->addCoProcessor(writeMQSelectCoProcessor);
      calculator->addCoProcessor(writeMQCoProcessor);
      calculator->addCoProcessor(timeAveragingCoProcessor);

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
   //Sleep(30000);

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

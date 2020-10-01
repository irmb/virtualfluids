#include <iostream>
#include <string>

#include <VirtualFluids.h>

using namespace std;


void pflowdp(string configname)
{
   try
   {
      ConfigurationFile   config;
      config.load(configname);

      string          pathname = config.getString("pathname");
      int             numOfThreads = config.getValue<int>("numOfThreads");
      vector<int>     blocknx = config.getVector<int>("blocknx");
      vector<double>  gridnx = config.getVector<double>("gridnx");
      double          nuLB = config.getValue<double>("nuLB");
      double          endTime = config.getValue<double>("endTime");
      double          outTime = config.getValue<double>("outTime");
      double          availMem = config.getValue<double>("availMem");
      int             refineLevel = config.getValue<int>("refineLevel");
      bool            logToFile = config.getValue<bool>("logToFile");
      double          restartStep = config.getValue<double>("restartStep");
      double          dpLB = config.getValue<double>("dpLB");
      bool            thinWall = config.getValue<bool>("thinWall");
      double          deltax = config.getValue<double>("deltax");


      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      LBMReal rhoLB = 0.0;
      double rhoLBinflow = dpLB*3.0;

      SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

      const int baseLevel = 0;

            //bounding box
            double g_minX1 = 0;
            double g_minX2 = 0;
            double g_minX3 = 0;
      
            double g_maxX1 = gridnx[0];
            double g_maxX2 = gridnx[1];
            double g_maxX3 = gridnx[2];

      double blockLength = (double)blocknx[0]*deltax;

      double h = (g_maxX2) / 2.0;
      double dex = g_maxX1;
      double Umax = (1.0 / (4.0*nuLB))*(dpLB / dex)*(h*h);
      double Re = (4 * h*Umax) / (3 * nuLB);

      //bc
      SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
      noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapterInflow(new DensityBCAdapter(rhoLBinflow));
      denBCAdapterInflow->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

      SPtr<BCAdapter> denBCAdapterOutflow(new DensityBCAdapter(rhoLB));
      denBCAdapterOutflow->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));

      //BS visitor
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(noSlipBCAdapter);
      bcVisitor.addBC(denBCAdapterInflow);
      bcVisitor.addBC(denBCAdapterOutflow);

      SPtr<Grid3D> grid(new Grid3D(comm));
      grid->setPeriodicX1(false);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(true);
      grid->setDeltaX(deltax);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

      double k1 = 4;
      double k2 = 8;

      SPtr<GbObject3D> refineCube1_1(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2 / k1 - 1.0, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube1_1.get(), pathname + "/geo/refineCube1_1", WbWriterVtkXmlBinary::getInstance());

      SPtr<GbObject3D> refineCube1_2(new GbCuboid3D(g_minX1, g_maxX2 - g_maxX2 / k1 + 1.0, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube1_2.get(), pathname + "/geo/refineCube1_2", WbWriterVtkXmlBinary::getInstance());

      SPtr<GbObject3D> refineCube2_1(new GbCuboid3D(g_minX1 + 2 * blockLength + 2 * deltax, g_minX2, g_minX3, g_maxX1 - 2 * blockLength - 2 * deltax, g_maxX2 / k2 - 1.0, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube2_1.get(), pathname + "/geo/refineCube2_1", WbWriterVtkXmlBinary::getInstance());

      SPtr<GbObject3D> refineCube2_2(new GbCuboid3D(g_minX1 + 2 * blockLength + 2 * deltax, g_maxX2 - g_maxX2 / k2 + 1.0, g_minX3, g_maxX1 - 2 * blockLength - 2 * deltax, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube2_2.get(), pathname + "/geo/refineCube2_2", WbWriterVtkXmlBinary::getInstance());

      SPtr<GbObject3D> refineCube2_3(new GbCuboid3D(g_minX1 + blockLength + 2 * deltax, g_minX2 + blockLength + 2 * deltax, g_minX3 + blockLength + 2 * deltax, g_maxX1 - blockLength - 2 * deltax, g_maxX2 - blockLength - 2 * deltax, g_maxX3 - blockLength - 2 * deltax));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube2_3.get(), pathname + "/geo/refineCube2_3", WbWriterVtkXmlBinary::getInstance());

      SPtr<GbObject3D> refineCube3(new GbCuboid3D(g_minX1 + 2.0*(g_maxX1/3.0), g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(refineCube3.get(), pathname + "/geo/refineCube3", WbWriterVtkXmlBinary::getInstance());

      //////////////////////////////////////////////////////////////////////////
      //restart
      //SPtr<UbScheduler> rSch(new UbScheduler(restartStep));
      //RestartCoProcessor rp(grid, rSch, comm, pathname, RestartCoProcessor::TXT);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep() == 0)
      {
         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         if (myid == 0)
         {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "h = " << h);
            UBLOG(logINFO, "rho = " << rhoLB);
            UBLOG(logINFO, "nue = " << nuLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << deltax);
            UBLOG(logINFO, "dpLB = " << dpLB);
            UBLOG(logINFO, "Umax = " << Umax);
            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
            UBLOG(logINFO, "Preprozess - start");
         }

         //walls
         GbCuboid3DPtr addWallYmin (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_minX3-4.0*blockLength, g_maxX1+4.0*blockLength, g_minX2, g_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax (new GbCuboid3D(g_minX1-4.0*blockLength, g_maxX2, g_minX3-4.0*blockLength, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1 - 4.0*blockLength, g_minX2 - 4.0*blockLength, g_minX3 - 4.0*blockLength, g_minX1, g_maxX2 + 4.0*blockLength, g_maxX3 + 4.0*blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - 4.0*blockLength, g_minX3 - 4.0*blockLength, g_maxX1 + 4.0*blockLength, g_maxX2 + 4.0*blockLength, g_maxX3 + 4.0*blockLength));
         if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         ////inflow
         //GbCuboid3DPtr geoInflow (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_minX3-4.0*blockLength, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_minX3));
         //if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         ////outflow
         //GbCuboid3DPtr geoOutflow (new GbCuboid3D(g_minX1-4.0*blockLength, g_minX2-4.0*blockLength, g_maxX3, g_maxX1+4.0*blockLength, g_maxX2+4.0*blockLength, g_maxX3+4.0*blockLength));
         //if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         SPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

         if (refineLevel > 0)
         {
            if (myid == 0) UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
            //refineHelper.addGbObject(refineCube1_1, 1);
            //refineHelper.addGbObject(refineCube1_2, 1);
            //refineHelper.addGbObject(refineCube2_1, 2);
            //refineHelper.addGbObject(refineCube2_2, 2);
            //refineHelper.addGbObject(refineCube2_3, refineLevel);
            refineHelper.addGbObject(refineCube3, refineLevel);
            
            refineHelper.refine();
            if (myid == 0) UBLOG(logINFO, "Refinement - end");
         }

         //walls
         SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, noSlipBCAdapter,Interactor3D::SOLID));
         SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, noSlipBCAdapter,Interactor3D::SOLID));

         SPtr<D3Q27Interactor> inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, denBCAdapterInflow, Interactor3D::SOLID));

         //outflow
         SPtr<D3Q27Interactor> outflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, denBCAdapterOutflow, Interactor3D::SOLID));

         ////////////////////////////////////////////
         //METIS
         SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         ////////////////////////////////////////////
         /////delete solid blocks
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - start");
         InteractorsHelper intHelper(grid, metisVisitor);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();
         if (myid == 0) UBLOG(logINFO, "deleteSolidBlocks - end");
         //////////////////////////////////////

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         ppblocks->process(0);
         ppblocks.reset();

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nod = nob * (blocknx[0]) * (blocknx[1]) * (blocknx[2]);
         unsigned long nodg = nob * (blocknx[0] + gl) * (blocknx[1] + gl) * (blocknx[1] + gl);
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

         int kernelType = 2;
         SPtr<LBMKernel> kernel;
         kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel());
         //}

         SPtr<BCProcessor> bcProc(new BCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //walls
         intHelper.setBC();

         grid->accept(bcVisitor);

         //initialization of distributions
         mu::Parser fct;
         fct.SetExpr("-(1.0/(2.0*nu))*(dp/dx)*((x2-h)^2 - h^2)");
         fct.DefineConst("dp", dpLB);
         fct.DefineConst("dx", dex);
         fct.DefineConst("h", h);
         fct.DefineConst("nu", nuLB);

         InitDistributionsBlockVisitor initVisitor;
         //initVisitor.setVx1(fct);
         initVisitor.setVx1(0.0);
         //initVisitor.setVx3(fct);
         grid->accept(initVisitor);

         //Postrozess
         SPtr<UbScheduler> geoSch(new UbScheduler(1));
         SPtr<CoProcessor> ppgeo(
            new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
         ppgeo->process(0);
         ppgeo.reset();

         if (myid == 0) UBLOG(logINFO, "Preprozess - end");
      }
      else
      {
         grid->accept(bcVisitor);

         //set connectors
         InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
         SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept(setConnsVisitor);

         if (myid == 0) UBLOG(logINFO, "Restart - end");
      }
      SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
      SPtr<NUPSCounterCoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      SPtr<UbScheduler> stepSch(new UbScheduler(outTime));

      SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));


      omp_set_num_threads(numOfThreads);
      SPtr<Calculator> calculator(new BasicCalculator(grid, stepSch, endTime));
      calculator->addCoProcessor(npr);
      calculator->addCoProcessor(writeMQCoProcessor);

      if (myid == 0) UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      if (myid == 0) UBLOG(logINFO, "Simulation-end");
   }
   catch (std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch (std::string& s)
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
         //pflowForcing(string(argv[1]));
         pflowdp(string(argv[1]));
      }
      else
      {
         cout << "Configuration file is missing!" << endl;
      }
   }

   return 0;
}

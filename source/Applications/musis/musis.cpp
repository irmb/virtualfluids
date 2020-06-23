#include <iostream>
#include <string>

#include <vfluids.h>

using namespace std;


void run(const char *cstr1, const char *cstr2)
{
   
   try
   {
      ConfigFileReader cf(cstr1);
      if ( !cf.read() )
      {
         std::string exceptionText = "Unable to read configuration file\n";
         throw exceptionText;
      }

      string machine = QUOTEME(CAB_MACHINE);
      string pathname = cf.getValue("path"); 
      int numOfThreads = UbSystem::stringTo<int>(cf.getValue("numOfThreads"));
      double availMem = 0;
      string geoFile;

      CommunicatorPtr comm(new MPICommunicator());
      int myid = comm->getProcessID();

      int d1, d2, d3; // for reading musis sample 

      if(machine == "BOMBADIL") 
      //if(machine == "YWANG")
      {
         //pathname = "J:/TBL/scratch/C100_DrySampleTest/";
         //geoFile  = "J:/TBL/TBL_Sw_Geos/CS518/TBL_CS518_MulSim02.geo.00000000.raw";
         pathname = cf.getValue("path"); 
         geoFile  = cf.getValue("geoFile"); 
         numOfThreads = UbSystem::stringTo<int>(cf.getValue("numOfThreads"));
         
         availMem = 3.0e9;
         d1 = UbSystem::stringTo<int>(cf.getValue("geoDimX1"));;
         d2 = UbSystem::stringTo<int>(cf.getValue("geoDimX2"));;
         d3 = UbSystem::stringTo<int>(cf.getValue("geoDimX3"));;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         //pathname = "/hpc3lustre/work/wang/TBL/scratch/CS518_DrySampleTest/";
         //geoFile = "/hpc3lustre/work/wang/TBL/TBL_Sw_Geos/CS518.X2Y2Z1/TBL_CS518_MulSim02.geo.00030000.raw";
         //numOfThreads = 1;
         pathname = cf.getValue("path"); 
         geoFile  = cf.getValue("geoFile"); 
         numOfThreads = UbSystem::stringTo<int>(cf.getValue("numOfThreads"));
         availMem = 12.0e9;

         if(myid ==0)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
         d1 = UbSystem::stringTo<int>(cf.getValue("geoDimX1"));;
         d2 = UbSystem::stringTo<int>(cf.getValue("geoDimX2"));;
         d3 = UbSystem::stringTo<int>(cf.getValue("geoDimX3"));;
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");
      
      const int baseLevel = 0;
      const int refineLevel = UbSystem::stringTo<int>(cf.getValue("refineLevel"));//2;
      const int blocknx1    = UbSystem::stringTo<int>(cf.getValue("blocknx1"));//12; 
      const int blocknx2    = UbSystem::stringTo<int>(cf.getValue("blocknx1"));//12;
      const int blocknx3    = UbSystem::stringTo<int>(cf.getValue("blocknx1"));//12;

      const int numBaseBlockL1 = UbSystem::stringTo<int>(cf.getValue("numBaseBlock_L1"));//1;
      ////////////////////////////////////////////////////////////////////////////
      //// Geometrie
      ////////////////////////////////////////////////////////////////////////////
      double L1 = UbSystem::stringTo<double>(cf.getValue("L1"));//0.07; //m
      double L2 = UbSystem::stringTo<double>(cf.getValue("L2"));//0.07; //m
      double L3 = UbSystem::stringTo<double>(cf.getValue("L3"));//0.0379 + 0.379; //m
      double dx = L1/(blocknx1*numBaseBlockL1)/(pow(2.0,refineLevel)); //0.0379/272.0; //m

      LBMReal rhoReal = 1.0; //kg/m^3
      LBMReal uReal = 5.0;  //m/s
      LBMReal uLB = 0.1;
      LBMReal nueLB = UbSystem::stringTo<double>(cf.getValue("nueLB"));//0.00166666666667;
      LBMReal Re = 0.0;
      LBMReal rhoLB = 0.0;

      //LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter(1.0, 1/sqrt(3.0)*(uReal/uLB), 1.0, 1.0/dx, dx*dx*dx));
      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());
      
      //bounding box
      double d_minX1 = 0.0;
      double d_minX2 = 0.0;
      double d_minX3 = 0.0;

      double d_maxX1 = L1;
      double d_maxX2 = L2;
      double d_maxX3 = L3;

      double offs = 0.0;

      double g_minX1 = d_minX1-offs;
      double g_minX2 = d_minX2-offs;;
      double g_minX3 = d_minX3-offs;

      double g_maxX1 = d_maxX1+offs;
      double g_maxX2 = d_maxX2+offs;
      double g_maxX3 = d_maxX3+offs;

      double blockLength = blocknx1 * dx;

      Grid3DPtr grid(new Grid3D(comm));
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(true);
      grid->setPeriodicX3(false);
      grid->setDeltaX(pow(2.0,refineLevel)*dx); // for coarse
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);

      double restartDump = UbSystem::stringTo<double>(cf.getValue("restartDump"));
      UbSchedulerPtr rSch(new UbScheduler(restartDump));
      RestartPostprocessorPtr rp(new RestartPostprocessor(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY));
      //UbSchedulerPtr emSch(new UbScheduler(1000, 1000));
      //EmergencyExitPostprocessor em(grid, emSch, pathname+"/checkpoints/emex.txt", rp, comm);

      std::string opt;
      if(cstr2!= NULL)
         opt = std::string(cstr2);

      LBMKernel3DPtr kernel;
      double ForcingX1 = UbSystem::stringTo<double>(cf.getValue("ForcingX1"));

      mu::Parser fctForcingX1;
      mu::Parser fctForcingX2;
      mu::Parser fctForcingX3;

      fctForcingX2.SetExpr("0.0");
      fctForcingX3.SetExpr("0.0");
      fctForcingX1.SetExpr("c3*(tanh(c1*(x3-c2))+c4)*Fx1*dx");
      //fctForcingX1.SetExpr("Fx1*dx");
      fctForcingX1.DefineConst("Fx1", ForcingX1);
      fctForcingX1.DefineConst("c1", 0.5);       // incline
      double ForcingLevel = 0.039/dx;
      fctForcingX1.DefineConst("c2", ForcingLevel); // forcing switch level
      fctForcingX1.DefineConst("c3", 0.5); // const always
      fctForcingX1.DefineConst("c4", 1.0); // const always
      if(myid == 0) UBLOG(logINFO,"Forcing Level = " << ForcingLevel );

      if(cstr2!= NULL)
      {
         if(myid==0) UBLOG(logINFO,"Restart step: " << opt);
         grid = rp->restart(UbSystem::stringTo<int>(opt));
         //rp->reconnect(grid);

         //Forcing setzen falls nötig
         SetForcingBlockVisitor forcingVisitor(fctForcingX1,fctForcingX2,fctForcingX3);
         grid->accept(forcingVisitor); 

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );

         ////domain decomposition //useful if pro mpi processor contains more than 1 thread 
         //PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         //grid->accept(pqPartVisitor);

         //int option = 0;
         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, option));


        
      }
      else
      {
         if(myid ==0)
         {
            UBLOG(logINFO,"L = " << L1/dx );
            UBLOG(logINFO,"v = " << uLB );
            UBLOG(logINFO,"rho = " << rhoLB );
            UBLOG(logINFO,"nue = " << nueLB );
            UBLOG(logINFO,"Re = " << Re );
            UBLOG(logINFO,"dx = " << dx );
            //UBLOG(logINFO,conv->toString() );
            UBLOG(logINFO,"number of levels = " << refineLevel+1 );
            UBLOG(logINFO,"numOfThreads = " << numOfThreads );
            UBLOG(logINFO,"Preprozess - start");
         }

         // read musis geometry 
         //if(myid ==0) UBLOG(logINFO,"Read geometry: start");
         //GbVoxelMatrix3DPtr vmatrix(new GbVoxelMatrix3D(d1, d2, d3, float(GbVoxelMatrix3D::FLUID),8.0,8.0)); 

         //vmatrix->readMatrixFromRawFile<char>(geoFile);
         //if(myid ==0) UBLOG(logINFO,"Read geometry: end");

         //vmatrix->setVoxelMatrixDelta(L1/(d1-1),L1/(d1-1),L1/(d1-1));

         //if(myid ==0) UBLOG(logINFO,"Write geometry: start");
         //if(myid == 0) vmatrix->writeToLegacyVTK(pathname+"/geo/geo");
         //if(myid ==0) UBLOG(logINFO,"Write geometry: end");
         
         // domain
         GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
         if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

         //refinement area
         GbObject3DPtr refineCube1(new  GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+8.0*0.0379));
         if(myid ==0) GbSystem3D::writeGeoObject(refineCube1.get(),pathname + "/geo/refineCube1", WbWriterVtkXmlBinary::getInstance());
         GbObject3DPtr refineCube2(new  GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+4.0*0.0379));
         if(myid ==0) GbSystem3D::writeGeoObject(refineCube2.get(),pathname + "/geo/refineCube2", WbWriterVtkXmlBinary::getInstance());

         // walls
         GbCuboid3DPtr addWallZmin (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());
         GbCuboid3DPtr addWallZmax (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());
         
         GenBlocksGridVisitor genBlocks;
         genBlocks.addGeoObject(gridCube);
         grid->accept(genBlocks);

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));

         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");   
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            refineHelper.addGbObject(refineCube1, refineLevel-1);
            refineHelper.addGbObject(refineCube2, refineLevel);
            refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");   
         }
         
         MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B);
         grid->accept( metisVisitor );

         SolidBlocksHelper sd(grid, comm);

         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bc_noslip(new D3Q27NoSlipBCAdapter(bbOption));
         D3Q27BoundaryConditionAdapterPtr bc_slip(new D3Q27SlipBCAdapter(bbOption));
         // porous geometry
         //D3Q27InteractorPtr geoInt = D3Q27InteractorPtr ( new D3Q27Interactor(vmatrix, grid, bc_noslip,Interactor3D::SOLID));

         //mu::Parser fct; 
         //fct.DefineConst("U", uLB);//Vx
         //fct.SetExpr("U"); 
         
         //D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST)); 

         //walls
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bc_noslip,Interactor3D::SOLID));
         ////up velocity
         //D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, velBCAdapter,Interactor3D::SOLID)); 
         //up slip
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bc_slip  ,Interactor3D::SOLID));
         
         //sd.addInteractor(geoInt);
         sd.addInteractor(addWallZminInt);
         sd.addInteractor(addWallZmaxInt);
      
         sd.deleteSolidBlocks();

         grid->accept( metisVisitor );


         ppblocks->update(0);
         ppblocks.reset();

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nodb = (blocknx1) * (blocknx2) * (blocknx3);
         unsigned long nod = nob * (blocknx1) * (blocknx2) * (blocknx3);
         unsigned long nodg = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);
         double needMemAll  = double(nod*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());

         if(myid == 0)
         {
            UBLOG(logINFO,"Number of blocks = " << nob);
            UBLOG(logINFO,"Number of nodes  = " << nod);
            int minInitLevel = grid->getCoarsestInitializedLevel();
            int maxInitLevel = grid->getFinestInitializedLevel();
            for(int level = minInitLevel; level<=maxInitLevel; level++)
            {
               int nobl = grid->getNumberOfBlocks(level);
               UBLOG(logINFO,"Number of blocks for level " << level <<" = " << nobl);
               UBLOG(logINFO,"Number of nodes for level " << level <<" = " << nobl*nodb);
            }
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");

         }            

         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3));
         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));
         //option = 0 - ohne param., option = 1 - mit param.
         //int option = 0;
         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, option));
         

         int kernelType = UbSystem::stringTo<int>(cf.getValue("kernel"));
         LBMKernel3DPtr kernel;
         if (kernelType == 0)
         {
            rhoLB = 1.0;
            kernel = LBMKernel3DPtr(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));
         }
         else if (kernelType == 1)
         {
            rhoLB = 1.0;
            kernel = LBMKernel3DPtr(new LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3));
         }
         else if (kernelType == 2)
         {
            rhoLB = 0.0;
            kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, 0));
            //kernel = LBMKernel3DPtr(new LBMKernelESD3Q27CCLB(blocknx1, blocknx2, blocknx3, grid));
            //kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLBex(blocknx1, blocknx2, blocknx3, 0, grid));
         }
         
         kernel->setForcingX1(fctForcingX1);
         kernel->setForcingX2(fctForcingX2);
         kernel->setForcingX3(fctForcingX3);
         kernel->setWithForcing(true);
         

         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);

         mu::Parser fctnueLB;

         SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

      //   //walls
         grid->addAndInitInteractor(addWallZminInt);
         grid->addAndInitInteractor(addWallZmaxInt);
         // porous geometry
         //grid->addAndInitInteractor(geoInt);

         //initialization of distributions
         D3Q27ETInitDistributionsBlockVisitor initVisitor(rhoLB);
         //initVisitor.setVx1(0.0);
         grid->accept(initVisitor);

         ////Postrozess - Measurement
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, comm, true));
         ppgeo->update(0);
         ppgeo.reset();

         //if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 
      }
      
      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      double calcBegin = UbSystem::stringTo<double>(cf.getValue("calcBegin"));//0.07; //m
      double calcEnd = UbSystem::stringTo<double>(cf.getValue("calcEnd"));//0.07; //m
      double calcIntervall = UbSystem::stringTo<double>(cf.getValue("calcIntervall"));//0.0379 + 0.379; //m
      
      /*UbSchedulerPtr TBL_Sch(new UbScheduler(calcIntervall,calcBegin,calcEnd));
      UbSchedulerPtr TBL_rSch(new UbScheduler(100000));
      TurbulentStrengthSurfaceRoughnessPostprocessor TBLpp(grid,pathname +"/results/TBL", TBL_Sch,TBL_rSch,comm);
      */

      double outTime = UbSystem::stringTo<double>(cf.getValue("outTime"));
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname + "/steps/step", WbWriterVtkXmlASCII::getInstance(), conv, comm);

      double fdx = grid->getDeltaX(grid->getFinestInitializedLevel());

      //D3Q27IntegrateValuesHelperPtr h1(new D3Q27IntegrateValuesHelper(grid, comm, 
      //   g_minX1, g_minX2, g_minX3, 
      //   g_minX1+1.0*fdx, g_maxX2, g_maxX3));
      ////if(myid ==0) GbSystem3D::writeGeoObject(h1->getBoundingBox().get(),pathname + "/geo/iv1", WbWriterVtkXmlBinary::getInstance());
      //D3Q27IntegrateValuesHelperPtr h2(new D3Q27IntegrateValuesHelper(grid, comm, 
      //   g_maxX1-1.0*fdx, g_minX2, g_minX3, 
      //   g_maxX1, g_maxX2, g_maxX3));
      ////if(myid ==0) GbSystem3D::writeGeoObject(h2->getBoundingBox().get(),pathname + "/geo/iv2", WbWriterVtkXmlBinary::getInstance());
      //LBMReal rhoReal = rhoLB;
      //LBMReal uReal = uLB; 
      //D3Q27PressureDifferencePostprocessor rhopp(grid, stepSch, pathname + "/results/rho_diff.txt", h1, h2, rhoReal, uReal, uLB, comm);

      //UbSchedulerPtr resSch(new UbScheduler(1000,10000,10000));
      //UbSchedulerPtr visSch(new UbScheduler(1,0,1));         
      //AverageValuesPostprocessor TBLpp(grid, pathname + "/results/AvVelocity", WbWriterVtkXmlBinary::getInstance(), visSch, resSch, comm); 


      double endTime = UbSystem::stringTo<double>(cf.getValue("endTime"));;//10001.0;

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
      if(myid == 0) UBLOG(logINFO,"Simulation-start");
      calculation->calculate();
      if(myid == 0) UBLOG(logINFO,"Simulation-end");
      //
      //double point1[3] = {0.45, 0.20, 0.205};
      //double point2[3] = {0.55, 0.20, 0.205};
      //D3Q27IntegrateValuesHelperPtr h1(new D3Q27IntegrateValuesHelper(grid, comm, 
      //   point1[0]-1.0*fdx, point1[1]-1.0*fdx, point1[2]-1.0*fdx, 
      //   point1[0], point1[1], point1[2]));
      //if(myid ==0) GbSystem3D::writeGeoObject(h1->getBoundingBox().get(),pathname + "/geo/iv1", WbWriterVtkXmlBinary::getInstance());
      //D3Q27IntegrateValuesHelperPtr h2(new D3Q27IntegrateValuesHelper(grid, comm, 
      //   point2[0], point2[1]-1.0*fdx, point2[2]-1.0*fdx, 
      //   point2[0]+1.0*fdx, point2[1], point2[2]));
      //if(myid ==0) GbSystem3D::writeGeoObject(h2->getBoundingBox().get(),pathname + "/geo/iv2", WbWriterVtkXmlBinary::getInstance());
      ////D3Q27PressureDifferencePostprocessor rhopp(grid, visSch, pathname + "/results/rho_diff.txt", h1, h2, conv, comm);
      //D3Q27PressureDifferencePostprocessor rhopp(grid, visSch, pathname + "/results/rho_diff.txt", h1, h2, rhoReal, uReal, uLB, comm);
      //
      //double area = 2.0*radius*H;
      //double v    = 4.0*uLB/9.0;
      //D3Q27ForcesPostprocessor fp(grid, visSch, pathname + "/results/forces.txt", comm, rhoLB, v, area, D3Q27ForcesPostprocessor::X, D3Q27ForcesPostprocessor::Y);
      //fp.addInteractor(cylinderInt);
      //
      //UbSchedulerPtr nupsSch(new UbScheduler(10, 10, 40));
      //NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      //double endTime = 40001.0;
      //CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      //if(myid == 0) UBLOG(logINFO,"Simulation-start");
      //calculation->calculate();
      //if(myid == 0) UBLOG(logINFO,"Simulation-end");
   }
   catch(std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch(std::string& s)
   {
      cerr << s << endl;
   }
   catch(...)
   {
      cerr << "unknown exception" << endl;
   }

}
int main(int argc, char* argv[])
{
   if ( argv != NULL )
   {
      if (argc > 1)
      {
         run(argv[1], argv[2]);
      }
      else
      {
         cout << "Configuration file must be set!: " <<  argv[0] << " <config file>" << endl << std::flush;
      }
   }

   return 0;
 
}





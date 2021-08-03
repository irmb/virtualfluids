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

      CommunicatorPtr comm = vf::mpi::MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      if(machine == "BOMBADIL") 
      {
         numOfThreads = 4;
         pathname = "d:/temp/sbone2";
         pathGeo = "d:/Data/Bone/SmallBone";
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         numOfThreads = 8;
         pathname = cf.getValue("pathname"); //"/work/koskuche/Bone/SmallBone";
         pathGeo = cf.getValue("pathGeo"); //"/home/koskuche/data/Bone/SmallBone/vti";
         availMem = 1.0e9;

#if defined(__unix__)
         if (myid == 0)
         {
            const char* str = pathname.c_str();
            int status = mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         }
#endif 

         if(myid ==0)
         {
            stringstream logFilename;
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+".txt";
            UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      if(myid==0) UBLOG(logINFO,"Testcase small bone");

      //string boneFileName = pathGeo + "/sbone.stl";
      string boneFileName = pathGeo + "/boneimage.vti";

      double dx = 3.5e-3/175.0;

      const int blocknx1 = 16;
      const int blocknx2 = 16;
      const int blocknx3 = 16;

      LBMReal rho_LB = 0.0;
      //nueWasser = 1e-6 m^2/s
      double nu_real = 1e-6;
      LBMReal dt = 5e-8; // s (frei gew�hlt)
      //dx - frei gew�hlt
      //
      LBMReal nu_LB = nu_real/(dx*dx/dt);


      //dp = 50000 Pa - 0 Pa = 50000 Pa
      double dp_real = UbSystem::stringTo<double>(cf.getValue("pressure")); //5000;
      //rho wasser = 1000 kg*m^-3
      double rho_real = 1000;
      //dp/rho = 50000/1000 = 50 m^2/s^2
      double dp_div_rho_real = dp_real/rho_real;

      double dp_LB = dp_div_rho_real/((dx/dt)*(dx/dt));

      bool with_forcing = false;

      double rhoLBinflow;
      if (with_forcing)
      {
         rhoLBinflow = 0.0;
      } 
      else
      {
         rhoLBinflow = dp_LB*3.0;
      }

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = 0;


      //////////////////////////////////////////////////////////////////////////
      //bone STL
      //GbTriFaceMesh3DPtr bone (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(boneFileName,"Netz"));
      //if(myid == 0) GbSystem3D::writeGeoObject( bone.get(), pathname+"/geo/bone", WbWriterVtkXmlBinary::getInstance() );

      string boneFilename = pathGeo + "/boneimage.vti";

      int pmNX1=151;  //abmessung einzelbild in x-richtung
      int pmNX2=101; //abmessung einzelbild in y richtung
      int pmNX3=101; //anzahl der bilder
      float lthreshold = 1.0;
      float uthreshold = 255.0;

      GbVoxelMatrix3DPtr bone(new GbVoxelMatrix3D(pmNX1,pmNX2,pmNX3,0,lthreshold,uthreshold));
      bone->readMatrixFromVtiASCIIFile(boneFilename);
      bone->setVoxelMatrixMininum(11.5, 8.01, 5.01);

      double deltax = dx*1e3;
      double deltaVoxel = 11e-3;
      bone->setVoxelMatrixDelta(deltaVoxel, deltaVoxel, deltaVoxel);
      bone->setLbGridDx(deltax);

      if(myid == 0) bone->writeToLegacyVTKBinary(pathname+"/geo/bone");

      //bounding box
      double g_minX1 = bone->getX1Minimum()-0.25;
      double g_minX2 = bone->getX2Minimum()-0.25;
      double g_minX3 = bone->getX3Minimum()-0.25;

      double g_maxX1 = bone->getX1Maximum()+0.25;
      double g_maxX2 = bone->getX2Maximum()+0.25;
      double g_maxX3 = bone->getX3Maximum()+0.25;

      double blockLength = (double)blocknx1*deltax;

      //double h = g_maxX2/2.0;
      //double dpLB = (rhoLBinflow - rhoLB)/3.0;

      //
      //double dex = g_maxX1+1.0;
      //double Umax = (1.0/(2.0*nueLB))*(dpLB/dex)*(h*h);

      //double Re = (4*h*Umax)/(3*nueLB);

      Grid3DPtr grid(new Grid3D(comm));
      grid->setPeriodicX1(false);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(false);
      grid->setDeltaX(deltax);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);

      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());      


      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      double forcing = 0;
      if (with_forcing)
      {
         forcing = dp_LB/(blocknx1*grid->getNX1());
      }

      if(myid ==0)
      {
         UBLOG(logINFO,"Parameters:");
         UBLOG(logINFO,"with forcing = " << with_forcing );
         UBLOG(logINFO,"rho_LB = " << rho_LB );
         UBLOG(logINFO,"nu_LB = " << nu_LB );
         UBLOG(logINFO,"dp_LB = " << dp_LB );
         UBLOG(logINFO,"forcing = " << forcing );
         UBLOG(logINFO,"dx = " << dx << " m");
         UBLOG(logINFO,"dt = " << dt << " s");
         UBLOG(logINFO,"rho_real = " << rho_real << " kg*m^-3" );
         UBLOG(logINFO,"nu_real = " << nu_real << " m^2/s" );
         UBLOG(logINFO,"dp_real = " << dp_real << " Pa" );

         UBLOG(logINFO,"number of levels = " << refineLevel+1 );
         UBLOG(logINFO,"numOfThreads = " << numOfThreads );
         UBLOG(logINFO,"path = " << pathname );
         UBLOG(logINFO,"Preprozess - start");
      }

      //walls
      GbCuboid3DPtr addWallYmin (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_minX2, g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallYmax (new GbCuboid3D(g_minX1-blockLength, g_maxX2, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallZmin (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallZmax (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

      //inflow
      GbCuboid3DPtr geoInflow (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_minX1, g_maxX2+blockLength, g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

      //outflow
      GbCuboid3DPtr geoOutflow (new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));


      //   if (refineLevel > 0)
      //   {
      //      if(myid == 0) UBLOG(logINFO,"Refinement - start");	
      //      RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
      //      refineHelper.refine();
      //      if(myid == 0) UBLOG(logINFO,"Refinement - end");	
      //   }



      //bone interactor
      int bcOptionBone = 2; //0=simple Bounce Back, 1=quadr. BB, 2=thin wall
      D3Q27BoundaryConditionAdapterPtr bcBone(new D3Q27NoSlipBCAdapter(bcOptionBone));

      D3Q27InteractorPtr boneInt(new D3Q27Interactor(bone, grid, bcBone,Interactor3D::SOLID));

      //wall interactors
      int bcOptionWall = 1; //0=simple Bounce Back, 1=quadr. BB, 2=thin wall
      D3Q27BoundaryConditionAdapterPtr bcWall(new D3Q27NoSlipBCAdapter(bcOptionWall));
      D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcWall,Interactor3D::SOLID));
      D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcWall,Interactor3D::SOLID));
      D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcWall,Interactor3D::SOLID));
      D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcWall,Interactor3D::SOLID));

      //   //inflow
      //   //double dp_Ph=0.1*10000.0;//dp in Bar
      //   //double dp_lb=dp_Ph*0.001*(nueLB*dx)*(nueLB*dx);//nue_ph=10e-6
      //   //if(myid == 0) UBLOG(logINFO,"dp_lb = " << dp_lb );
      //   //double rhoLBinflow = 3.0*(dp_lb-rhoLB);

      D3Q27BoundaryConditionAdapterPtr denBCAdapterInflow(new D3Q27DensityBCAdapter(rhoLBinflow));
      denBCAdapterInflow->setSecondaryBcOption(0);
      D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, denBCAdapterInflow, Interactor3D::SOLID));

      //outflow
      D3Q27BoundaryConditionAdapterPtr denBCAdapterOutflow(new D3Q27DensityBCAdapter(rho_LB));
      denBCAdapterOutflow->setSecondaryBcOption(0);
      D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapterOutflow,Interactor3D::SOLID));

      ////////////////////////////////////////////
      //METIS
      Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));   
      ////////////////////////////////////////////
      /////delete solid blocks
      if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
      InteractorsHelper intHelper(grid, metisVisitor);
      intHelper.addInteractor(boneInt);
      intHelper.addInteractor(addWallYminInt);
      intHelper.addInteractor(addWallYmaxInt);
      intHelper.addInteractor(addWallZminInt);
      intHelper.addInteractor(addWallZmaxInt);
      intHelper.addInteractor(inflowInt);
      intHelper.addInteractor(outflowInt);
      intHelper.selectBlocks();
      if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	 
      //////////////////////////////////////

      //set connectors
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu_LB, iProcessor);
      grid->accept( setConnsVisitor );

      //domain decomposition for threads
      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      grid->accept(pqPartVisitor);

      ppblocks->update(0);
      ppblocks.reset();

      unsigned long nob = grid->getNumberOfBlocks();
      int gl = 3;
      unsigned long nodb = (blocknx1) * (blocknx2) * (blocknx3);
      unsigned long nod = nob * (blocknx1) * (blocknx2) * (blocknx3);
      unsigned long nodg = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);
      double needMemAll  = double(nodg*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
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

      LBMKernel3DPtr kernel = LBMKernel3DPtr(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));

      //mu::Parser fctForcingX1;
      //fctForcingX1.SetExpr("Fx1");
      //fctForcingX1.DefineConst("Fx1", forcing);

      //kernel->setForcingX1(fctForcingX1);
      //kernel->setWithForcing(true);

      BCProcessorPtr bcProc(new D3Q27ETForThinWallBCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nu_LB, availMem, needMem);
      grid->accept(kernelVisitor);


      //if (refineLevel > 0)
      //{
      //   D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
      //   grid->accept(undefNodesVisitor);
      //}

      //BC
      intHelper.setBC();

      //Press*1.6e8+(14.76-coordsX)/3.5*5000
      //initialization of distributions
      mu::Parser fct;
      fct.SetExpr("(x1max-x1)/l*dp*3.0");
      fct.DefineConst("dp", dp_LB);
      fct.DefineConst("x1max", g_maxX1);
      fct.DefineConst("l", g_maxX1-g_minX1);

      D3Q27ETInitDistributionsBlockVisitor initVisitor(nu_LB, rho_LB);
      initVisitor.setRho(fct);
      //initVisitor.setVx1(fct);
      initVisitor.setVx1(0.0);
      grid->accept(initVisitor);

      //Postrozess
      UbSchedulerPtr geoSch(new UbScheduler(1));
      D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
         new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
      ppgeo->update(0);
      ppgeo.reset();

      if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);

      double outTime = 1000;
      UbSchedulerPtr stepSch(new UbScheduler(outTime));
      stepSch->addSchedule(10,10,100);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv);


      double dxd2 = deltax / 2.0;
      D3Q27IntegrateValuesHelperPtr ih1(new D3Q27IntegrateValuesHelper(grid, comm, bone->getX1Minimum() - dxd2, bone->getX2Minimum() - dxd2, bone->getX3Minimum() - dxd2,
         bone->getX1Maximum() + dxd2, bone->getX2Maximum() + dxd2, bone->getX3Maximum() + dxd2));
      if (myid == 0) GbSystem3D::writeGeoObject(ih1->getBoundingBox().get(), pathname + "/geo/ih1", WbWriterVtkXmlBinary::getInstance());

      double factorp = dp_real/dp_LB;
      double factorv = dx/dt;
      D3Q27MeanValuesPostprocessor mvp1(grid, stepSch, pathname + "/mv/mv1.txt", comm, ih1, factorp, factorv);


      D3Q27IntegrateValuesHelperPtr ih2(new D3Q27IntegrateValuesHelper(grid, comm, g_maxX1-2.0*deltax, g_minX2, g_minX3,
         g_maxX1 - deltax, g_maxX2, g_maxX3));
      if (myid == 0) GbSystem3D::writeGeoObject(ih2->getBoundingBox().get(), pathname + "/geo/ih2", WbWriterVtkXmlBinary::getInstance());

      D3Q27MeanValuesPostprocessor mvp2(grid, stepSch, pathname + "/mv/mv2.txt", comm, ih2, factorp, factorv);

      if(myid == 0)
      {
         UBLOG(logINFO,"PID = " << myid << " Total Physical Memory (RAM): " << Utilities::getTotalPhysMem());
         UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used: " << Utilities::getPhysMemUsed());
         UBLOG(logINFO,"PID = " << myid << " Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe());
      }

      double endTime = UbSystem::stringTo<double>(cf.getValue("endTime")); //100001;//10001.0;

      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, stepSch));
      if(myid == 0) UBLOG(logINFO,"Simulation-start");
      calculation->calculate();
      if(myid == 0) UBLOG(logINFO,"Simulation-end");
   }
   catch(exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch(string& s)
   {
      cerr << s << endl;
   }
   catch(...)
   {
      cerr << "unknown exception" << endl;
   }

}
//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

   if ( argv != NULL )
   {
      sbonepd(argv[1]);
   }

   return 0;
}

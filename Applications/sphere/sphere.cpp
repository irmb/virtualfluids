#include <vfluids.h>
#include <set>
#include <map>
using namespace std;


////////////////////////////////////////////////////////////////////////
void chanel(const char *cstr)
{
   try
   {

      //Sleep(30000);

      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
      double availMem = 0;

      CommunicatorPtr comm = MPICommunicator::getInstance();

      int myid = comm->getProcessID();
      int mybundle = comm->getBundleID();
      int root = comm->getRoot();

      ConfigFileReader cf(cstr);
      if ( !cf.read() )
      {
         std::string exceptionText = "Unable to read configuration file\n";
         throw exceptionText;
      }

      pathname = cf.getValue("path");
      availMem = UbSystem::stringTo<double>(cf.getValue("memory"));
      string metafile = cf.getValue("metafile");
      double outstep = UbSystem::stringTo<double>(cf.getValue("outstep"));
      double endstep = UbSystem::stringTo<double>(cf.getValue("endstep"));
      int numOfThreads = UbSystem::stringTo<int>(cf.getValue("threads"));

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      double dx = 1;

      const int blocknx1 = 8;
      const int blocknx2 = 8;
      const int blocknx3 = 8;

      const int gridNx1 = 8;// 18;
      const int gridNx2 = 4;// 11;
      const int gridNx3 = 4;// 11;

      //const int blocknx1 = 40;
      //const int blocknx2 = 40;
      //const int blocknx3 = 40;

      //const int gridNx1 = 2;
      //const int gridNx2 = 2;
      //const int gridNx3 = 2;

      double L1 = gridNx1*blocknx1;
      double L2, L3;
      L2 = gridNx2*blocknx1;
      L3 = gridNx3*blocknx1;

      LBMReal radius = 2.5;
      LBMReal uLB = 0.0000001;
      LBMReal Re = 0.0001;
      LBMReal rhoLB = 0.0;
      //LBMReal nuLB = (uLB*2.0*radius)/Re;
      //LBMReal nuLB = (uLB*L2)/Re;
      LBMReal nuLB = 0.168666666667/100;

      double dp_LB =1e-6;
      double rhoLBinflow = dp_LB*3.0;

      Grid3DPtr grid(new Grid3D(comm));
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);

      BoundaryConditionProcessorPtr bcProcessor(new BoundaryConditionProcessor());

      //////////////////////////////////////////////////////////////////////////
      //restart
      UbSchedulerPtr restartSch(new UbScheduler(100000, 100000, 100000));
      RestartPostprocessor rp(grid, restartSch, comm, pathname, RestartPostprocessor::BINARY);
      //////////////////////////////////////////////////////////////////////////

      if (grid->getTimeStep() == 0)
      {

         const int baseLevel = 0;
         const int refineLevel = UbSystem::stringTo<int>(cf.getValue("level"));

         //bounding box
         double d_minX1 = 0.0;
         double d_minX2 = 0.0;
         double d_minX3 = 0.0;

         double d_maxX1 = L1;
         double d_maxX2 = L2;
         double d_maxX3 = L3;

         double blockLength = blocknx1*dx;

         if(myid ==0)
         {
            UBLOG(logINFO,"Parameters:");
            UBLOG(logINFO,"uLB = " << uLB );
            UBLOG(logINFO,"rhoLB = " << rhoLB );
            UBLOG(logINFO,"nueLB = " << nuLB );
            UBLOG(logINFO,"Re = " << Re );
            UBLOG(logINFO,"dx = " << dx );
            UBLOG(logINFO,"number of levels = " << refineLevel+1 );
            UBLOG(logINFO,"numOfThreads = " << numOfThreads );
            UBLOG(logINFO,"Preprozess - start");
         }

         GbObject3DPtr gridCube(new GbCuboid3D(d_minX1, d_minX2, d_minX3, d_maxX1, d_maxX2, d_maxX3));
         if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance()); 

         GenBlocksGridVisitor genBlocks(gridCube);
         grid->accept(genBlocks);

         //sphere
         //GbObject3DPtr sphereRef(new GbSphere3D(L1/4.0, L2*0.5, L3*0.5, radius+1.0));
         //GbSystem3D::writeGeoObject(sphereRef.get(),pathname + "/geo/sphereRef", WbWriterVtkXmlBinary::getInstance());

         //sphere
         GbObject3DPtr sphere(new GbSphere3D(L1/4.0, L2*0.5, L3*0.5, radius));
         //GbObject3DPtr sphere(new GbCuboid3D(L1/4.0-radius, L2/2.0-radius, L3/2.0-radius, L1/4.0+radius, L2/2.0+radius, L3/2.0+radius));
         GbSystem3D::writeGeoObject(sphere.get(),pathname + "/geo/sphere", WbWriterVtkXmlBinary::getInstance());

         double off = 0.0;
         GbObject3DPtr refCube(new GbCuboid3D(sphere->getX1Minimum()-off, sphere->getX2Minimum()-off, sphere->getX3Minimum(), 
            sphere->getX1Maximum()+off, sphere->getX2Maximum()+off, sphere->getX3Maximum()));
         if(myid ==0) GbSystem3D::writeGeoObject(refCube.get(),pathname + "/geo/refCube", WbWriterVtkXmlBinary::getInstance()); 

         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel);
            //refineHelper.addGbObject(sphere, refineLevel);
            refineHelper.addGbObject(refCube, refineLevel);
            refineHelper.refine();
            if(myid == 0) UBLOG(logINFO,"Refinement - end");	
         }

         //walls
         GbCuboid3DPtr addWallYmin (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_minX2, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmin (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_minX3));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_maxX2, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_maxX3, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow (new GbCuboid3D(d_minX1-4.0*blockLength, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_minX1, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow (new GbCuboid3D(d_maxX1, d_minX2-4.0*blockLength, d_minX3-4.0*blockLength, d_maxX1+4.0*blockLength, d_maxX2+4.0*blockLength, d_maxX3+4.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));



         //sphere
         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcNoSlip(new D3Q27NoSlipBCAdapter(bbOption));
         D3Q27InteractorPtr sphereInt = D3Q27InteractorPtr ( new D3Q27Interactor(sphere, grid, bcNoSlip,Interactor3D::SOLID));

         D3Q27BoundaryConditionAdapterPtr bcSlip(new D3Q27SlipBCAdapter(bbOption));

         //walls
         D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcNoSlip,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcNoSlip,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcNoSlip,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcNoSlip,Interactor3D::SOLID));

         mu::Parser fct;
         fct.SetExpr("U");
         fct.DefineConst("U", uLB);
 
         //inflow
         //D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
         //velBCAdapter->setSecondaryBcOption(2);
         //D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         D3Q27BoundaryConditionAdapterPtr denBCAdapterInflow(new D3Q27DensityBCAdapter(rhoLBinflow));
         denBCAdapterInflow->setSecondaryBcOption(0);
         D3Q27InteractorPtr inflowInt = D3Q27InteractorPtr(new D3Q27Interactor(geoInflow, grid, denBCAdapterInflow, Interactor3D::SOLID));

         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));

         Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B));
         InteractorsHelper intHelper(grid, metisVisitor);
         //intHelper.addInteractor(sphereInt);
         intHelper.addInteractor(addWallYminInt);
         intHelper.addInteractor(addWallZminInt);
         intHelper.addInteractor(addWallYmaxInt);
         intHelper.addInteractor(addWallZmaxInt);
         intHelper.addInteractor(inflowInt);
         intHelper.addInteractor(outflowInt);
         intHelper.selectBlocks();

         //domain decomposition for threads
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);


         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         grid->accept( setConnsVisitor );

         Block3DConnectorFactoryPtr factory(new Block3DConnectorFactory());
         //ConnectorBlockVisitor setConnsVisitor(comm, nuLB, iProcessor, factory);
         //grid->accept(setConnsVisitor);

         ppblocks->update(0);
         ppblocks.reset();

         unsigned long nob = grid->getNumberOfBlocks();
         int gl = 3;
         unsigned long nod = nob * (blocknx1+gl) * (blocknx2+gl) * (blocknx3+gl);

         double needMemAll  = double(nod*(27*sizeof(double) + sizeof(int) + sizeof(float)*4));
         double needMem  = needMemAll / double(comm->getNumberOfProcesses());

         if(myid == 0)
         {
            UBLOG(logINFO,"Number of blocks = " << nob);
            UBLOG(logINFO,"Number of nodes  = " << nod);
            UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
            UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
            UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
         }            

         LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB(blocknx1, blocknx2, blocknx3, LBMKernelETD3Q27CCLB::NORMAL));

         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         //BoundaryConditionPtr velocityBC(new VelocityBoundaryCondition(refineLevel+1));
         BoundaryConditionPtr densityBC(new NonEqDensityBoundaryCondition());
         //BoundaryConditionPtr noSlipBC(new NoSlipBoundaryCondition());

         //BoundaryConditionPtr velocityBC(new NonReflectingVelocityBoundaryCondition());
         //BoundaryConditionPtr densityBC(new NonReflectingDensityBoundaryCondition());
         BoundaryConditionPtr noSlipBC(new HighViscosityNoSlipBoundaryCondition());

         //bcProcessor->addBC(velocityBC);
         bcProcessor->addBC(densityBC);
         bcProcessor->addBC(noSlipBC);

         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nuLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         intHelper.setBC();


         BoundaryConditionBlockVisitor bcVisitor;
         grid->accept(bcVisitor);

         mu::Parser fctRoh;
         fctRoh.SetExpr("(x1max-x1)/l*dp*3.0");
         fctRoh.DefineConst("dp", dp_LB);
         fctRoh.DefineConst("x1max", d_maxX1);
         fctRoh.DefineConst("l", d_maxX1-d_minX1);

         //initialization of distributions
         D3Q27ETInitDistributionsBlockVisitor initVisitor(nuLB, rhoLB);
         //initVisitor.setVx1(fct);
         initVisitor.setRho(fctRoh);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, true));
         ppgeo->update(0);
         ppgeo.reset();

         if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 
      }
      else
      {
         UBLOG(logINFO,"SetConnectors - start, id="<<myid);

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27IncompressibleOffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
         //ConnectorFactoryPtr cFactory(new Block3DConnectorFactory());
         //ConnectorBlockVisitor setConnsVisitor(comm, nuLB, iProcessor, cFactory);
         grid->accept( setConnsVisitor );

         UBLOG(logINFO,"SetConnectors - end, id="<<myid); 
      }

      UbSchedulerPtr stepSch(new UbScheduler(outstep));
      //stepSch->addSchedule(10000, 0, 1000000);
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv);

      //D3Q27MacroscopicQuantitiesPostprocessor ppg(grid, stepSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv,true);

      //InSituVTKPostprocessor isp(grid, stepSch, metafile, conv);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 30, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, numOfThreads, comm);


      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endstep, stepSch, bcProcessor));

      if(myid == 0) 
         UBLOG(logINFO,"Simulation-start");

      calculation->calculate();

      if(myid == 0) 
         UBLOG(logINFO,"Simulation-end");

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

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   //struct iNode
   //{
   //   int x1, x2, x3;
   //   LBMReal x1off, x2off, x3off;
   //};

   //typedef boost::shared_ptr<iNode> iNodePtr;

  
   //

   //set <vector<int>>  iNodesSetSender00;
   //vector<double> offs;
   //map <vector<int>>  iNodesSetSender00;

   //{
   //   iNodePtr inode1(new iNode());
   //   inode1->x1 = 1;
   //   inode1->x3 = 2;
   //   inode1->x2 = 3;
   //   inode1->x1off = 0;
   //   inode1->x2off = 1;
   //   inode1->x3off = 2;

   //   iNodesSetSender00.insert(inode1);
   //}

   //{
   //   iNodePtr inode1(new iNode());
   //   inode1->x1 = 1;
   //   inode1->x3 = 2;
   //   inode1->x2 = 3;
   //   inode1->x1off = 0;
   //   inode1->x2off = 1;
   //   inode1->x3off = 2;

   //   iNodesSetSender00.insert(inode1);
   //}

   //std::pair<std::set<vector<int>>::iterator, bool> ret;

   //   {
   //      vector<int> inv;
   //      inv.push_back(1);
   //      inv.push_back(2);
   //      inv.push_back(3);
   //      inv.push_back(1);
   //      inv.push_back(2);
   //      inv.push_back(3);
   //      ret=iNodesSetSender00.insert(inv);
   //      if (ret.second==true)
   //      {
   //         offs.push_back(0);
   //         offs.push_back(1);
   //         offs.push_back(2);
   //      }
   //   }

   //         {
   //            vector<int> inv;
   //            inv.push_back(1);
   //            inv.push_back(1);
   //            inv.push_back(3);
   //            inv.push_back(1);
   //            inv.push_back(2);
   //            inv.push_back(3);
   //            ret=iNodesSetSender00.insert(inv);
   //            if (ret.second==true)
   //            {
   //               offs.push_back(0);
   //               offs.push_back(1);
   //               offs.push_back(2);
   //            }
   //         }

   //         return 0;


   if ( argv != NULL )
   {
      if (argc > 1)
      {
         chanel(argv[1]);
      }
      else
      {
         cout << "Configuration file must be set!: " <<  argv[0] << " <config file>" << endl << std::flush;
      }
   }
}


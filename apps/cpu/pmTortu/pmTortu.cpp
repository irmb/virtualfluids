#include <iostream>
#include <string>

#include <vfluids.h>

using namespace std;

//! \brief  Computes Flow thorugh a porous medium and writespathlines
//! \details Aim: determine tortuosity. pathlines are later integrated using python-script streamlinesMean.py (needs scipy,numpy)
//! \details If PM-data is large run on single visulalization node.
//! \details Created on: 01.07.2013
//! \author  Sonja Uphoff

void run(const char *cstr)
{
   try
   {
     string machine = QUOTEME(CAB_MACHINE);
      UBLOG(logINFO,"Testcase PMTortuosity");
      string pathname;
      string stlPath;
      int numOfThreads =1;
      bool logfile = false;
      stringstream logFilename;
      double availMem = 0;

      UbLog::reportingLevel() = logDEBUG5; //logINFO;

      CommunicatorPtr comm(new MPICommunicator());
      int myid = comm->getProcessID();

      if(machine == "PIPPINNEU")
      {

         pathname = "f:/temp/PMtortu";
         stlPath = "f:/GeometrienVliese";
         numOfThreads = 3;
         logfile = false;
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")
      {
         pathname = "/work/sonjaOutputs/PMTortu2metall450";
         stlPath = "/work/sonjaOutputs";
         numOfThreads = 4;
         availMem = 12.0e9;
         logfile = true;

         //if(myid ==0)
         //{
            logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(myid)+".txt";
         //}
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");



      //if(myid ==0 && logfile)
      //{
         UbLog::output_policy::setStream(logFilename.str());
      //}

      int baseLevel, refineLevel,nx[3],blocknx[3];
      double Re,velocity,rhoInit,vx1Init;//,vx2Init,vx3Init;

//////////////////////////////////////////////////////////////////////////
      //physik
//////////////////////////////////////////////////////////////////////////
      Re            = 1;// 13286;//13286;//gemessen 18.98 m/s...*5.0 zum  testen ob was passiert
      velocity      = 0.01;
      vx1Init       = 0.01;
      rhoInit       = 1.0;
      SimulationParametersPtr param = SimulationParameters::getInstanz();
param->setCollisionModelType(SimulationParameters::COMPRESSIBLE);

      ///////////////Knotenabmessungen:

    nx[0]=28;
    nx[1]=27;
       nx[2]=27;
    blocknx[0]=10;
    blocknx[1]=10;
    blocknx[2]=10;

   baseLevel   = 0;
   refineLevel = 0;

      bool periodicx1 = false;
      bool periodicx2 = false;
      bool periodicx3 = false;



   double minX1 = 0.0;
   double maxX1 = 280;
   double minX2 = 0.0;
   double maxX2 = 270;
   double minX3 = 0.0;
   double maxX3 = 270;
   double centerX1 = 0.5*(maxX1-minX1);
   double centerX2 = 0.5*(maxX2-minX2);
   //double scaleAsphalt = 0.0000625; //10/1600
   double scalepixeltomm=0.5;
   double scaleAsphalt = 1.0;
   minX1 = minX1*scaleAsphalt;
   minX2 = minX2*scaleAsphalt;
   minX3 = minX3*scaleAsphalt;
   maxX1 = maxX1*scaleAsphalt;
   maxX2 = maxX2*scaleAsphalt;
   maxX3 = maxX3*scaleAsphalt;

   //vorgabe geom. dx im feinsten = 1 -> abstand der voxel = 1
   double coarseNodeDx = (maxX2 - minX2) / (double)( blocknx[1]*nx[1] );
   double fineNodeDx   = coarseNodeDx / (double)(1<<refineLevel);

   double blockLengthx1 = blocknx[0]*coarseNodeDx;
   double blockLengthx2 = blocknx[1]*coarseNodeDx;
   double blockLengthx3 = blocknx[2]*coarseNodeDx;

   double originX1 = minX1;
   double originX2 = minX2;
   double originX3 = minX3;

   int nx1 = nx[0];
   int nx2 = nx[1];
   int nx3 = nx[2];
   int blocknx1      = blocknx[0];
   int blocknx2      = blocknx[1];
   int blocknx3      = blocknx[2];

   double gridOrigin[3] = { originX1, originX2, originX3 };

   //geom. GROBE Blocklaenge
   double coarseBlockLength[3];
   coarseBlockLength[0] = blockLengthx1;
   coarseBlockLength[1] = blockLengthx2;
   coarseBlockLength[2] = blockLengthx3;
   double geoLength[]   = {  nx[0]*blockLengthx1, nx[1]*blockLengthx2, nx[2]*blockLengthx3};

//////////////////////////////////////////////////////////////////////////
   // PM File
//////////////////////////////////////////////////////////////////////////
   string pmFilename;
   pmFilename = stlPath+"/metallrgbx271y271z270.vti";//
   int pmNX1=270;
   int pmNX2=271;
   int pmNX3=270;
   float threshold = 120.0;

         GbVoxelMatrix3DPtr pmMesh(GbVoxelMatrix3DCreator::getInstance()->createFromVtiASCIIFloatFile(pmFilename,pmNX1,pmNX2,pmNX3,threshold));

pmMesh->translate((maxX1-minX1)*0.05,-(maxX2-minX2)*0.01,-(maxX3-minX3)*0.01);
   pmMesh->setTransferViaFilename(true, pmFilename);

//##########################################################################
      //## physical parameters
//##########################################################################

      double rhoLB         = 1.0;
      double rhoReal       = 1.0;
      double nueReal  = 0.16;//0.015;

      double hReal         = maxX1;
      double uReal         = Re*nueReal/hReal;

      //##Machzahl:
      //#Ma     = uReal/csReal

      double csReal  = 1.0/sqrt(3.0);
      double cs_LB=1.0/sqrt(3.0);
      double Ma      = uReal/csReal;//0.0553;//Ma-Real!
      double hLB     = hReal;

      //LBMUnitConverter unitConverter(hReal, csReal, rhoReal, hLB);
      LBMUnitConverterPtr unitConverter = LBMUnitConverterPtr(new LBMUnitConverter(hReal, csReal, rhoReal, blocknx[0]*nx[0] ));

      double uLB           = uReal   * unitConverter->getFactorVelocityWToLb();
      double nueLB         = nueReal * unitConverter->getFactorViscosityWToLb();

      double realDeltaT     = (nueLB * hReal *hReal) / (nueReal * blocknx[0]*nx[0] *blocknx[0]*nx[0]);



      Grid3DPtr grid(new Grid3D());
      UbSchedulerPtr rSch(new UbScheduler(5000,5000,1000000));
      RestartPostprocessor rp(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY);

//////////////////////////////////////////////////////////////////////////

     std::string opt;

      if(cstr!= NULL)
         opt = std::string(cstr);

      //bounding box
      double g_minX1 = originX1;
      double g_minX2 = originX2;
      double g_minX3 = originX3;

      double g_maxX1 = originX1 + geoLength[0];
      double g_maxX2 = originX2 + geoLength[1];
      double g_maxX3 = originX3 + geoLength[2];

      //set grid
      grid->setDeltaX(coarseNodeDx);
      grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
      grid->setPeriodicX1(periodicx1);
      grid->setPeriodicX2(periodicx2);
      grid->setPeriodicX3(periodicx3);


      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if(myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/geo/gridCube", WbWriterVtkXmlASCII::getInstance());

      GenBlocksGridVisitor genBlocks;
      genBlocks.addGeoObject(gridCube);
      grid->accept(genBlocks);


//////////////////////////////////////////////////////////////////////////
      if(myid == 0)
      {
         UBLOG(logINFO, "*****************************************");
         UBLOG(logINFO, "* Parameters *");
         UBLOG(logINFO, "* Re            ="<<Re);
         UBLOG(logINFO, "* Ma            ="<<Ma);
         UBLOG(logINFO, "* uReal         ="<<uReal);
         UBLOG(logINFO, "* nueReal       ="<<nueReal);
         UBLOG(logINFO, "* nue           ="<<nueLB);
         UBLOG(logINFO, "* velocity      ="<<uLB);
         UBLOG(logINFO, "* LX1 (world/LB)="<<hReal<<"/"<<hReal/coarseNodeDx);
      //   UBLOG(logINFO, "* LX2 (world/LB)="<<kanalbreiteSI<<"/"<<kanalbreiteSI/coarseNodeDx);
      //   UBLOG(logINFO, "* LX3 (world/LB)="<<kanalhoeheSI<<"/"<<kanalhoeheSI/coarseNodeDx);
         UBLOG(logINFO, "* cdx           ="<<coarseNodeDx);
         UBLOG(logINFO, "* fdx           ="<<fineNodeDx);
         UBLOG(logINFO, "* dx_base ="<<coarseNodeDx<<" == "<<coarseNodeDx);
         UBLOG(logINFO, "* dx_refine ="<<fineNodeDx<<" == "<<fineNodeDx );
         UBLOG(logINFO, "* nx1/2/3 ="<<nx[0]<<"/"<<nx[1]<<"/"<<nx[2]);
         UBLOG(logINFO, "* blocknx1/2/3 ="<<blocknx[0]<<"/"<<blocknx[1]<<"/"<<blocknx[2]);
         UBLOG(logINFO, "* x2Periodic    ="<<periodicx2);
         UBLOG(logINFO, "* x3Periodic    ="<<periodicx3);
         UBLOG(logINFO, "*****************************************");
         UBLOGML(logINFO, "UnitConverter:"<<unitConverter->toString());
         UBLOG(logINFO, "*****************************************");
      }


      RatioBlockVisitor ratioVisitor(refineLevel);
      grid->accept(ratioVisitor);
      RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
      grid->accept(ratioSmoothVisitor);
      OverlapBlockVisitor overlapVisitor(refineLevel);
      grid->accept(overlapVisitor);
      std::vector<int> dirs;
      D3Q27System::getLBMDirections(dirs);
      SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
      grid->accept(interDirsVisitor);

      if(myid == 0) UBLOG(logINFO,"Refinement - end");

      MetisPartitioningGridVisitor metisVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::B);
      grid->accept( metisVisitor );

      if(myid == 0) UBLOG(logINFO,"Write blocks - start");
      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
      if(myid == 0) ppblocks->update(0);
      if(myid == 0) UBLOG(logINFO,"Write blocks - end");

      if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - start");
      SolidBlocksHelper sd(grid, comm);


      sd.deleteSolidBlocks();
      if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");

      if(myid == 0) UBLOG(logINFO,"Write blocks - start");
      grid->accept( metisVisitor );
      if(myid == 0) ppblocks->update(1);
      ppblocks.reset();
      if(myid == 0) UBLOG(logINFO,"Write blocks - end");

      unsigned long nob = grid->getNumberOfBlocks();
      unsigned long nod = nob * blocknx[0]*blocknx[1]*blocknx[2];
      unsigned long nod_real = nob * (blocknx[0]+3)*(blocknx[1]+3)*(blocknx[2]+3);

      double needMemAll  = double(nod_real*(27*sizeof(double) + sizeof(int)));
      double needMem  = needMemAll / double(comm->getNumberOfProcesses());

      if(myid == 0)
      {
         UBLOG(logINFO,"Number of blocks = " << nob);
         UBLOG(logINFO,"Number of nodes  = " << nod);
         UBLOG(logINFO,"Necessary memory  = " << needMemAll << " bytes");
         UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
         UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
      }

      LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(blocknx[0], blocknx[1], blocknx[2]));

  //    LBMKernel3DPtr kernel(new LBMKernelETD3Q27BGK(blocknx[0], blocknx[1], blocknx[2],1));
      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);


      SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);

      grid->accept(kernelVisitor);



//////////////////////////////////////////////////////////////////////////
    double geoOverlap = 5*coarseNodeDx;


//////////////////////////////////////////////////////////////////////////
   // Interactoren
//////////////////////////////////////////////////////////////////////////
//##########################################################################
   int noSlipSecOpt = 0; // #0=2nd order BB 1=simple BB
//##########################################################################
   int noSlipSecOptAsphalt = 1; // #0=2nd order BB 1=simple BB
//##########################################################################
     int bbOption1 = 0; //0=simple Bounce Back, 1=quadr. BB
     D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption1));
   D3Q27InteractorPtr PM1Interactor = D3Q27InteractorPtr ( new D3Q27Interactor(pmMesh, grid, bcObst,Interactor3D::SOLID)); //wo ist bc obst definiert?
 grid->addAndInitInteractor( PM1Interactor);
   //UBLOG(logINFO,"SpD3Q19Asphalt - send porous media to D3Q19InteractorService");
   //UBLOG(logINFO,"SpD3Q19Asphalt - send porous media = "<<pmInteractor->getName()<<" with "<<typeid(*pmInteractor->getGbObject3D()).name()<<" node("<<pmNX1<<"/"<<pmNX2<<"/"<<pmNX3<<")");
   UbTimer timer;
   timer.start();


   UBLOG(logINFO,"SpD3Q19Asphalt - send porous media to D3Q19InteractorService done in "<<timer.stop());


      if (refineLevel > 0)
      {
         D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
         grid->accept(undefNodesVisitor);
      }


      //set connectors
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
      grid->accept( setConnsVisitor );

      //domain decomposition
      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      grid->accept(pqPartVisitor);

      //initialization of decompositions
      D3Q27ETInitDistributionsBlockVisitor initVisitor(1.0);
      initVisitor.setVx1(0.0);
      grid->accept(initVisitor);


 //////////////////////////////////////////////////////////////////////////
   // BCs
//////////////////////////////////////////////////////////////////////////
      //Reparatur an den Seiten:
       UBLOG(logINFO,"inflow")
           double x3minMesh=0.000;  double x3maxMesh=originX3 + blockLengthx3*nx3 + geoOverlap;
   GbCuboid3DPtr leftCuboid(new GbCuboid3D( originX1 + blockLengthx1*nx1 - coarseNodeDx,
                                           originX2 - geoOverlap,
                                           x3minMesh,
                                           originX1 + blockLengthx1*nx1 + geoOverlap,
                                           originX2 + blockLengthx2*nx2 + geoOverlap,
                                           x3maxMesh));
   GbCuboid3DPtr rightCuboid(new GbCuboid3D( originX1 - geoOverlap,
                                           originX2 - geoOverlap,
                                           x3minMesh,
                                           originX1 + geoOverlap,
                                           originX2 + blockLengthx2*nx2 + geoOverlap,
                                           x3maxMesh));
   GbCuboid3DPtr northCuboid(new GbCuboid3D( originX1- geoOverlap,
                                           originX2 + blockLengthx2*nx2 - 0.5*coarseNodeDx,
                                           x3minMesh,
                                           originX1 + blockLengthx1*nx1 + geoOverlap,
                                           originX2 + blockLengthx2*nx2 + geoOverlap,
                                           x3maxMesh));
   GbCuboid3DPtr southCuboid(new GbCuboid3D( originX1 - geoOverlap,
                                           originX2 - geoOverlap,
                                           x3minMesh,
                                           originX1 + blockLengthx1*nx1 + geoOverlap,
                                           originX2 + geoOverlap,
                                           x3maxMesh));

//////////////////////////////////////////////////////////////////////////
   // inflow
//////////////////////////////////////////////////////////////////////////
   UBLOG(logINFO,"inflow")

   GbCuboid3DPtr densCuboid(new GbCuboid3D(
                                           originX1 - geoOverlap,
                                           originX2 - geoOverlap,
                                           originX3 + blockLengthx3*nx3 - coarseNodeDx,
                                           originX1 + blockLengthx1*nx1 + geoOverlap,
                                           originX2 + blockLengthx2*nx2 + geoOverlap,
                                           originX3 + blockLengthx3*nx3 + geoOverlap));


//////////////////////////////////////////////////////////////////////////
   // bottom/outflow
//////////////////////////////////////////////////////////////////////////
      double dRho=0.05;
      GbCuboid3DPtr densCuboid2(new GbCuboid3D(
                                                 originX1 - geoOverlap,
                                                 originX2 - geoOverlap,
                                                 originX3 - geoOverlap,
                                                 originX1 + blockLengthx1*nx1 + geoOverlap,
                                                 originX2 + blockLengthx2*nx2 + geoOverlap,
minX3+0.5*fineNodeDx   ));

      if(myid == 0) GbSystem3D::writeGeoObject(densCuboid2.get(), pathname+"/geo/densCuboid2", WbWriterVtkXmlASCII::getInstance());
      D3Q27BoundaryConditionAdapterPtr denBCAdapter2(new D3Q27DensityBCAdapter(rhoInit-dRho));
      D3Q27InteractorPtr densInteractor2( new D3Q27Interactor(leftCuboid,grid,denBCAdapter2,Interactor3D::SOLID) );
      grid->addAndInitInteractor( densInteractor2 );

            if(myid == 0) GbSystem3D::writeGeoObject(densCuboid.get(), pathname+"/geo/densCuboid", WbWriterVtkXmlASCII::getInstance());
      D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoInit+dRho));
      D3Q27InteractorPtr densInteractor( new D3Q27Interactor(rightCuboid,grid,denBCAdapter,Interactor3D::SOLID) );
      grid->addAndInitInteractor( densInteractor );

   D3Q27InteractorPtr leftInteractor = D3Q27InteractorPtr ( new D3Q27Interactor(densCuboid2, grid, bcObst,Interactor3D::SOLID));
   grid->addAndInitInteractor( leftInteractor);
   D3Q27InteractorPtr rightInteractor = D3Q27InteractorPtr ( new D3Q27Interactor(densCuboid, grid, bcObst,Interactor3D::SOLID));
  grid->addAndInitInteractor(rightInteractor);
   D3Q27InteractorPtr northInteractor = D3Q27InteractorPtr ( new D3Q27Interactor(northCuboid, grid, bcObst,Interactor3D::SOLID));
   grid->addAndInitInteractor(northInteractor);
   D3Q27InteractorPtr southInteractor = D3Q27InteractorPtr ( new D3Q27Interactor(southCuboid, grid, bcObst,Interactor3D::SOLID));
  grid->addAndInitInteractor(southInteractor);

  if(myid == 0) GbSystem3D::writeGeoObject(northCuboid.get(), pathname+"/geo/north", WbWriterVtkXmlASCII::getInstance());
if(myid == 0) GbSystem3D::writeGeoObject(southCuboid.get(), pathname+"/geo/south", WbWriterVtkXmlASCII::getInstance());
if(myid == 0) GbSystem3D::writeGeoObject(rightCuboid.get(), pathname+"/geo/right", WbWriterVtkXmlASCII::getInstance());
if(myid == 0) GbSystem3D::writeGeoObject(leftCuboid.get(), pathname+"/geo/left", WbWriterVtkXmlASCII::getInstance());


      UbSchedulerPtr geoSch(new UbScheduler(1));
      D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
           new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(),
unitConverter, comm, true));

              double raiseVelSteps = 0;

      grid->doPostProcess(0);
      ppgeo.reset();
      geoSch.reset();

     UbSchedulerPtr plSch(new UbScheduler(10, 2));
      vector<D3Q27PathLinePostprocessorPtr> pathlinepostPParray;

 for (int ppz=0; ppz<27; ppz++)
      {
      for (int ppy=0; ppy<27; ppy++)
      {
          char numstr[21];
          sprintf(numstr, "%d", ppy+20*ppz);
          std::string pathPL = pathname+"/pathline" + numstr+".dat";
         D3Q27PathLinePostprocessorPtr plptr1( new D3Q27PathLinePostprocessor(grid, pathPL, WbWriterVtkXmlASCII::getInstance(), unitConverter, plSch, comm, 8.0, 6.0+8.0*(double)ppy,5.0+8.0*(double)ppz, nueLB, iProcessor));
              pathlinepostPParray.push_back(plptr1);//new D3Q27PathLinePostprocessor(grid, pathname + "/pathLine", WbWriterVtkXmlASCII::getInstance(), conv, plSch, comm, 0.01+(double)ppx*0.0001, 0.00001,0.00001, nueLB, iProcessor));

          }
      }

 UbSchedulerPtr visSch(new UbScheduler());
      visSch->addSchedule(1,1,10);
      visSch->addSchedule(10,10,100);
      visSch->addSchedule(100,100,1000);
      visSch->addSchedule(1000,1000,100000);
      visSch->addSchedule(100000,100000,1000000);
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), unitConverter, comm);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 10, 100));
      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

//////////////////////////////////////////////////////////////////////////

      cout << "PID = " << myid << " Total Physical Memory (RAM): " << MemoryUtil::getTotalPhysMem()<<endl;
      cout << "PID = " << myid << " Physical Memory currently used: " << MemoryUtil::getPhysMemUsed()<<endl;
      cout << "PID = " << myid << " Physical Memory currently used by current process: " << MemoryUtil::getPhysMemUsedByMe()<<endl;

      double endTime = 40001;
      UbSchedulerPtr ghostLSch(new UbScheduler());
      ghostLSch->addSchedule(1,1,endTime);
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, ghostLSch));
      if(myid == 0) UBLOG(logINFO,"Simulation-start");
      calculation->calculate();
      if(myid == 0) UBLOG(logINFO,"Simulation-end");
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

   run(argv[1]);

   return 0;
} 

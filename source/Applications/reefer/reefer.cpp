#include <iostream>
#include <string>

#include "numerics/geometry3d/CoordinateTransformation3D.h"
#include "Grid3D.h"
#include "GenBlocksGridVisitor.h"
#include "numerics/geometry3d/GbSystem3D.h"
#include "numerics/geometry3d/GbCuboid3D.h"
#include "numerics/geometry3d/GbCylinder3D.h"
#include <numerics/geometry3d/GbSphere3D.h>
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "OverlapBlockVisitor.h"
#include "RefineInterGbObjectsVisitor.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "SetKernelBlockVisitor.h"
#include "LBMKernelETD3Q27Cascaded.h"
#include "D3Q27MacroscopicQuantitiesPostprocessor.h"
#include "MPICommunicator.h"
#include "D3Q27ETBCProcessor.h"
#include "SimulationParameters.h"
#include "D3Q27SetUndefinedNodesBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include "D3Q27SetConnectorsBlockVisitor.h"
#include "NullCommunicator.h"
#include "D3Q27ETInitDistributionsBlockVisitor.h"
#include "CalculationManager.h"
#include "PQueuePartitioningGridVisitor.h"
#include "MetisPartitioningGridVisitor.h"
#include "D3Q27Interactor.h"
#include "D3Q27NoSlipBCAdapter.h"
#include "D3Q27VelocityBCAdapter.h"
#include "D3Q27DensityBCAdapter.h"
#include "D3Q27BoundaryConditionAdapter.h"
#include "StringUtil.hpp"
#include "D3Q27OffsetInterpolationProcessor.h"
#include "D3Q27CompactInterpolationProcessor.h"
#include "SyncBcBlockVisitor.h"
#include "numerics/geometry3d/creator/GbTriFaceMesh3DCreator.h"
#include "numerics/geometry3d/GbTriFaceMesh3D.h"
#include "D3Q27TriFaceMeshInteractor.h"
#include "MathUtil.hpp"
#include "SolidBlocksHelper.h"
#include "LBMKernelETD3Q27CascadedTI.h"
#include "TurbulenceIntensityPostprocessor.h"
#include "RestartPostprocessor.h"

using namespace std;


void run(const char *cstr)
{
   try
   {
      string pathname = "c:/temp/reefer/out";
      string pathGeo = "c:/Data/reefer";

      //string pathname = "/work/koskuche/scratch/reefer2/out";
      //string pathGeo = "/home/koskuche/data/reefer/new";

      //string pathname = "/home/kucher/temp/reefer/out";
      //string pathGeo = "/home/kucher/data/reefer/new";

      int numOfThreads = 2;
      
      CommunicatorPtr comm(new MPICommunicator());
      int myid = comm->getProcessID();

      //if(myid ==0)
      //{
      //   stringstream logFilename;
      //   logFilename <<  "/work/koskuche/scratch/reefer2/logfile.log";
      //   UbLog::output_policy::setStream(logFilename.str());
      //}

      //const double dx = 13.6;
      const double dx = 2.0;
      double refLentgthWorld = dx/1000.0; //from mm to m
      double refLentgthLB = 1.0;
      LBMUnitConverterPtr uconv = LBMUnitConverterPtr(new LBMUnitConverter(refLentgthWorld, LBMUnitConverter::AIR_20C, refLentgthLB));
      LBMReal uSI = 10;//m/s
      LBMReal uLB = uSI * uconv->getFactorVelocityWToLb();
      LBMReal rhoLB = 1.0;
      LBMReal nueSI = 1.5e-5;
      LBMReal nueLB = nueSI * uconv->getFactorViscosityWToLb();//(uLB*l)/Re;

      Grid3DPtr grid(new Grid3D());
      UbSchedulerPtr rSch(new UbScheduler(1500,5000));
      RestartPostprocessor rp(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY);

      std::string opt;

      if(cstr!= NULL)
         opt = std::string(cstr);

      if(cstr!= NULL)
      {
         opt = std::string(cstr);

         if(myid==0) UBLOG(logINFO,"Restart step: " << opt);

         grid = rp.restart(UbSystem::stringTo<int>(opt));

         LBMReal nueLB = 1.5e-3;
         
         SimulationParametersPtr param = SimulationParameters::getInstanz();
         param->setCollisionModelType(SimulationParameters::COMPRESSIBLE);

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );

         //domain decomposition
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);
      }
      else
      {
      const int baseLevel = 0;
      const int refineLevel = 0;
      //////////////////////////////////////////////////////////////////////////
      // Geometries
      //////////////////////////////////////////////////////////////////////////
      //container
      GbTriFaceMesh3DPtr geoContainer (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo +"/Containerascii.stl","geoContainer"));
      if(myid == 0) GbSystem3D::writeGeoObject(geoContainer.get(), pathname+"/geo/geoContainer", WbWriterVtkXmlASCII::getInstance());
      //cargo
      //GbTriFaceMesh3DPtr geoCargo (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo + "/Kisten_fuer_Palettenascii.stl","geoCargo"));
      GbTriFaceMesh3DPtr geoCargo (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo + "/Kistenascii.stl","geoCargo"));
      if(myid == 0) GbSystem3D::writeGeoObject(geoCargo.get(), pathname+"/geo/geoCargo", WbWriterVtkXmlASCII::getInstance());
      //palette
      //GbTriFaceMesh3DPtr geoPalette (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo + "/Palettenascii.stl","geoPalette"));
      //if(myid == 0) GbSystem3D::writeGeoObject(geoPalette.get(), pathname+"/geoPalette", WbWriterVtkXmlASCII::getInstance());
      //reefer
      GbTriFaceMesh3DPtr geoBlower (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo + "/Solidblockascii.stl","geoReefer"));
      if(myid == 0) GbSystem3D::writeGeoObject(geoBlower.get(), pathname+"/geo/geoBlower", WbWriterVtkXmlASCII::getInstance());
      //T floor
      GbTriFaceMesh3DPtr geoTFloor (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathGeo + "/T-Floorascii.stl","geoTFloor"));
      if(myid == 0) GbSystem3D::writeGeoObject(geoTFloor.get(), pathname+"/geo/geoTFloor", WbWriterVtkXmlASCII::getInstance());

      //bounding box
      double g_minX1 = geoContainer->getX1Minimum();
      double g_minX2 = geoContainer->getX2Minimum();
      double g_minX3 = geoContainer->getX3Minimum();

      double g_maxX1 = geoContainer->getX1Maximum();
      double g_maxX2 = geoContainer->getX2Maximum();
      double g_maxX3 = geoContainer->getX3Maximum();

      const int nodesPerBlock = 10;
      //const double dx = 1.7;
      //const double dx = 13.6;
      const double blockLength = double(nodesPerBlock)*dx;

      const double gridOriginX1 = g_minX1;
      const double gridOriginX2 = g_minX2;
      const double gridOriginX3 = g_minX3;

      //add wall X
      GbCuboid3DPtr addWallXmax (new GbCuboid3D(g_maxX1, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      //GbCuboid3DPtr addWallXmax (new GbCuboid3D(geoBlower->getX1Maximum()+geoBlower->getLengthX1(), g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallXmax.get(), pathname+"/geo/addWallXmax", WbWriterVtkXmlASCII::getInstance());
      //add wall Y
      GbCuboid3DPtr addWallYmax (new GbCuboid3D(g_minX1-blockLength, geoBlower->getX2Maximum(), g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());
      //add wall Z
      GbCuboid3DPtr addWallZmax (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_maxX3, g_maxX1+blockLength, g_maxX2+blockLength, g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());
      //add wall X
      GbCuboid3DPtr addWallXmin (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, geoBlower->getX1Minimum(), g_maxX2+blockLength, g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallXmin.get(), pathname+"/geo/addWallXmin", WbWriterVtkXmlASCII::getInstance());
      //add wall Y
      GbCuboid3DPtr addWallYmin (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, geoBlower->getX2Minimum(), g_maxX3+blockLength));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());
      //add wall Z
      GbCuboid3DPtr addWallZmin (new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3-blockLength, g_maxX1+blockLength, g_maxX2+blockLength, geoTFloor->getX3Minimum()));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());
      //add wall for blower
      GbCuboid3DPtr addWallBlower (new GbCuboid3D(geoBlower->getX1Minimum()-3.0*blockLength, geoBlower->getX2Minimum()-3.0*blockLength, geoBlower->getX3Minimum()+4.0*dx, 
                                                  geoBlower->getX1Maximum(), geoBlower->getX2Maximum()+3.0*blockLength, geoBlower->getX3Maximum()-4.0*dx));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallBlower.get(), pathname+"/geo/addWallBlower", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallBlowerXmin (new GbCuboid3D(geoBlower->getX1Minimum(), geoBlower->getX2Minimum(), geoBlower->getX3Minimum(), 
                                       geoBlower->getX1Minimum()+2.0*dx, geoBlower->getX2Maximum(), geoBlower->getX3Maximum()));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallBlowerXmin.get(), pathname+"/geo/addWallBlowerXmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallBlowerXmax (new GbCuboid3D(geoBlower->getX1Maximum()-2.0*dx, geoBlower->getX2Minimum(), geoBlower->getX3Minimum(), 
                                                      geoBlower->getX1Maximum(), geoBlower->getX2Maximum(), geoBlower->getX3Maximum()));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallBlowerXmax.get(), pathname+"/geo/addWallBlowerXmax", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallBlowerYmin (new GbCuboid3D(geoBlower->getX1Minimum(), geoBlower->getX2Minimum(), geoBlower->getX3Minimum(), 
                                                      geoBlower->getX1Maximum(), geoBlower->getX2Minimum()+2.0*dx, geoBlower->getX3Maximum()));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallBlowerYmin.get(), pathname+"/geo/addWallBlowerYmin", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWallBlowerYmax (new GbCuboid3D(geoBlower->getX1Minimum()-2.0*dx, geoBlower->getX2Maximum()-2.0*dx, geoBlower->getX3Minimum(), 
                                                      geoBlower->getX1Maximum(), geoBlower->getX2Maximum(), geoBlower->getX3Maximum()));
      if(myid == 0) GbSystem3D::writeGeoObject(addWallBlowerYmax.get(), pathname+"/geo/addWallBlowerYmax", WbWriterVtkXmlASCII::getInstance());

      //inflow
      GbCuboid3DPtr geoInflow (new GbCuboid3D(geoBlower->getX1Minimum()+dx, geoBlower->getX2Minimum()+dx, geoBlower->getX3Minimum()+2.0*dx, 
                                              geoBlower->getX1Maximum()-dx, geoBlower->getX2Maximum()-dx, geoBlower->getX3Minimum()+4.0*dx));
      if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

      //outflow
      GbCuboid3DPtr geoOutflow (new GbCuboid3D(geoBlower->getX1Minimum()+2.0*dx, geoBlower->getX2Minimum()+2.0*dx, geoBlower->getX3Maximum()-4.0*dx, 
                                               geoBlower->getX1Maximum()-2.0*dx, geoBlower->getX2Maximum()-2.0*dx, geoBlower->getX3Maximum()-2.0*dx));
      if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

      //simulation parameters
      double lSI = g_maxX2 - g_minX2;
      double lLB = lSI / dx;
      //double refLentgthWorld = blockLength/1000.0; //from mm to m
      //double refLentgthLB = double(nodesPerBlock);
      //LBMUnitConverterPtr uconv = LBMUnitConverterPtr(new LBMUnitConverter(refLentgthWorld, LBMUnitConverter::AIR_20C, refLentgthLB));
      //LBMReal uSI = 10;//m/s
      //LBMReal uLB = uSI * uconv->getFactorVelocityWToLb();
      //LBMReal rhoLB = 1.0;
      //LBMReal nueSI = 1.5e-5;
      //LBMReal nueLB = nueSI * uconv->getFactorViscosityWToLb();//(uLB*l)/Re;
      //LBMReal nueLB = 1.5e-3;
      LBMReal Re = (uLB*(420/dx))/nueLB;

      if(myid ==0)
      {
         UBLOG(logINFO,"grid = " <<int((g_maxX1 - g_minX1)/dx)<<"x"<<int((g_maxX2 - g_minX2)/dx) << "x"<<int((g_maxX3 - g_minX3)/dx));
         UBLOG(logINFO,"dx = " << dx);
         UBLOG(logINFO,"nodes per block = " << nodesPerBlock);
         UBLOG(logINFO,"block length = " << blockLength << "mm");
         UBLOG(logINFO,"v = " << uLB );
         UBLOG(logINFO,"rho = " << rhoLB );
         UBLOG(logINFO,"nue = " << nueLB );
         UBLOG(logINFO,"Re = " << Re );
         UBLOG(logINFO,"Preprozess - start");
      }

      SimulationParametersPtr param = SimulationParameters::getInstanz();
      param->setCollisionModelType(SimulationParameters::COMPRESSIBLE);
      param->setRho(rhoLB);
      param->setVelocityX(uLB);
      param->setViscosity(nueLB);
      
      //set grid
      //Grid3DPtr grid(new Grid3D());
      grid->setDeltaX(dx);
      grid->setBlockNX(nodesPerBlock, nodesPerBlock, nodesPerBlock);
      
      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      if(myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname+"/gridCube", WbWriterVtkXmlASCII::getInstance());

      GenBlocksGridVisitor genBlocks;
      genBlocks.addGeoObject(gridCube);
      grid->accept(genBlocks);

      MetisPartitioningGridVisitor metisVisitor(numOfThreads, D3Q27System::B, comm, false);
      grid->accept( metisVisitor );

      SolidBlocksHelper sd(grid, comm);

      //iteractors
      int bbOption1 = 0; //0=simple Bounce Back, 1=quadr. BB
      D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption1));

      D3Q27TriFaceMeshInteractorPtr cargoInt( new D3Q27TriFaceMeshInteractor(geoCargo, grid, bcObst,Interactor3D::SOLID));
      sd.addInteractor(cargoInt);

      D3Q27InteractorPtr addWallBlowerInt(new D3Q27Interactor(addWallBlower, grid, bcObst,Interactor3D::SOLID));

      sd.addInteractor(addWallBlowerInt);

      //D3Q27TriFaceMeshInteractorPtr paletteInt( new D3Q27TriFaceMeshInteractor(geoPalette, grid, bcObst,Interactor3D::SOLID));
      //sd.addInteractor(paletteInt);

      D3Q27InteractorPtr addWallXmaxInt(new D3Q27Interactor(addWallXmax, grid, bcObst,Interactor3D::SOLID));
      sd.addInteractor(addWallXmaxInt);

      D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcObst,Interactor3D::SOLID));
      sd.addInteractor(addWallYmaxInt);

      D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcObst,Interactor3D::SOLID));
      sd.addInteractor(addWallZmaxInt);

      D3Q27InteractorPtr addWallXminInt(new D3Q27Interactor(addWallXmin, grid, bcObst,Interactor3D::SOLID));
      sd.addInteractor(addWallXminInt);

      D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcObst,Interactor3D::SOLID));
      sd.addInteractor(addWallYminInt);

      D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcObst,Interactor3D::SOLID));
      sd.addInteractor(addWallZminInt);

      sd.deleteSolidBlocks();
      if(myid == 0) UBLOG(logINFO,"deleteSolidBlocks - end");	      
      
      if (refineLevel > 0)
      {
         GbObject3DPtr refineCube1(new  GbCuboid3D(geoTFloor->getX1Minimum(), geoTFloor->getX2Minimum(), geoTFloor->getX3Minimum(), 
            geoTFloor->getX1Maximum(), geoTFloor->getX2Maximum(), geoTFloor->getX3Maximum()));
         GbSystem3D::writeGeoObject(refineCube1.get(),pathname + "/refineCube", WbWriterVtkXmlASCII::getInstance());

         RefineCrossAndInsideGbObjectBlockVisitor refVisitor(refineCube1, baseLevel, refineLevel-1);
         grid->accept(refVisitor);

         RatioBlockVisitor ratioVisitor(refineLevel);
         grid->accept(ratioVisitor);

         RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
         grid->accept(ratioSmoothVisitor);

         OverlapBlockVisitor overlapVisitor(refineLevel);
         grid->accept(overlapVisitor);
      }

      if(myid == 0) UBLOG(logINFO,"Write blocks - start");
      grid->accept( metisVisitor );
      if(myid == 0) grid->writeBlocks(pathname + "/blocks" + StringUtil::toString(myid), 0, WbWriterVtkXmlASCII::getInstance(), false);
      if(myid == 0) UBLOG(logINFO,"Write blocks - end");

      unsigned long nob = grid->getNumberOfBlocks();
      unsigned long nod = nob * nodesPerBlock * nodesPerBlock *nodesPerBlock;
      double availMem = 6.0e9;
      double needMemAll  = double(nod*(27*sizeof(double) + sizeof(int))*2);
      double needMem  = needMemAll / double(comm->getNumberOfProcesses());

      if(myid == 0)
      {
         UBLOG(logINFO,"Number of blocks = " << nob);
         UBLOG(logINFO,"Number of nodes  = " << nod);
         UBLOG(logINFO,"Necessary memory  = " << needMemAll  << " bytes");
         UBLOG(logINFO,"Necessary memory per process = " << needMem  << " bytes");
         UBLOG(logINFO,"Available memory per process = " << availMem << " bytes");
      }

      //LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(nodesPerBlock, nodesPerBlock, nodesPerBlock));
      LBMKernel3DPtr kernel(new LBMKernelETD3Q27CascadedTI(nodesPerBlock, nodesPerBlock, nodesPerBlock));
      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
      //SetKernelBlockVisitor kernelVisitor(kernel, nueLB);
      grid->accept(kernelVisitor);

      if (refineLevel > 0)
      {
         std::vector<int> dirs;
         D3Q27System::getLBMDirections(dirs);
         SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
         grid->accept(interDirsVisitor);

         D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
         grid->accept(undefNodesVisitor);
      }

      //discretization
      //D3Q27TriFaceMeshInteractorPtr containerInt( new D3Q27TriFaceMeshInteractor(geoContainer, grid, bcObst,Interactor3D::SOLID));
      //grid->addAndInitInteractor(containerInt);

      D3Q27TriFaceMeshInteractorPtr tFloorInt( new D3Q27TriFaceMeshInteractor(geoTFloor, grid, bcObst,Interactor3D::SOLID));
      grid->addAndInitInteractor(tFloorInt);

      grid->addAndInitInteractor(addWallBlowerInt);
      //grid->addAndInitInteractor(blowerInt);
      grid->addAndInitInteractor(cargoInt);
      //grid->addAndInitInteractor(paletteInt);
      grid->addAndInitInteractor(addWallXmaxInt);
      grid->addAndInitInteractor(addWallYmaxInt);
      grid->addAndInitInteractor(addWallZmaxInt);
      grid->addAndInitInteractor(addWallXminInt);
      grid->addAndInitInteractor(addWallYminInt);
      grid->addAndInitInteractor(addWallZminInt);

      D3Q27InteractorPtr addWallBlowerXminInt(new D3Q27Interactor(addWallBlowerXmin, grid, bcObst,Interactor3D::SOLID));
      grid->addAndInitInteractor(addWallBlowerXminInt);
      
      D3Q27InteractorPtr addWallBlowerXmaxInt(new D3Q27Interactor(addWallBlowerXmax, grid, bcObst,Interactor3D::SOLID));
      grid->addAndInitInteractor(addWallBlowerXmaxInt);
      
      D3Q27InteractorPtr addWallBlowerYminInt(new D3Q27Interactor(addWallBlowerYmin, grid, bcObst,Interactor3D::SOLID));
      grid->addAndInitInteractor(addWallBlowerYminInt);
      
      D3Q27InteractorPtr addWallBlowerYmaxInt(new D3Q27Interactor(addWallBlowerYmax, grid, bcObst,Interactor3D::SOLID));
      grid->addAndInitInteractor(addWallBlowerYmaxInt);

      //outflow
      D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
      D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));
      grid->addAndInitInteractor(outflowInt);

      //inflow
      double Cx = geoInflow->getX1Centroid();
      double Hx = geoInflow->getLengthX1();
      double Cy = geoInflow->getX2Centroid();
      double Hy = geoInflow->getLengthX2();
      mu::Parser fct = MathUtil::getDuctParaboloidZ(Cx,Hx,Cy,Hy,-uLB);
      //mu::Parser fct;
      //fct.SetExpr("vx3");
      //fct.DefineConst("vx3", uLB);

      D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (false, false ,true ,fct, 0, D3Q27BCFunction::INFCONST));
      velBCAdapter->setSecondaryBcOption(2);
      D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));
      grid->addAndInitInteractor(inflowInt);

      //set connectors
      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
      grid->accept( setConnsVisitor );

      //domain decomposition
      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      grid->accept(pqPartVisitor);

      //initialization of decompositions
      D3Q27ETInitDistributionsBlockVisitor initVisitor(1.0);
      grid->accept(initVisitor);


      //Postrozess
      //if(myid == 0) grid->writeBlocks(pathname + "/blocks" + StringUtil::toString(myid), 0, WbWriterVtkXmlASCII::getInstance(), false);

      //std::vector< UbTupleFloat3 > nodes;
      //std::vector< UbTupleInt2 >   lines;
      //sphereInt->addQsLineSet(nodes, lines);
      //WbWriterVtkXmlBinary::getInstance()->writeLines(pathname+"/qs",nodes,lines);

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());
      
      UbSchedulerPtr geoSch(new UbScheduler(1));
      D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
           new D3Q27MacroscopicQuantitiesPostprocessor(grid, pathname + "/geo/nodes_geo", WbWriterVtkXmlBinary::getInstance(), 
                                                       conv, geoSch, comm, true));
      grid->doPostProcess(0);
      ppgeo.reset();
      geoSch.reset();
      

      if(myid == 0) UBLOG(logINFO,"Preprozess - end");      

}
      
      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());
      //double outTime = 50000;
      double outTime = 500;
      UbSchedulerPtr visSch(new UbScheduler());
      visSch->addSchedule(1000,1000,10000);
      visSch->addSchedule(10000,10000,100000);
      visSch->addSchedule(100000,100000,1000000);
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, pathname + "/mq/nodes", WbWriterVtkXmlBinary::getInstance(), conv, visSch, comm);

      //turbulence intensity postprocessor
      UbSchedulerPtr tiSch(new UbScheduler());
      tiSch->addSchedule(1000, 5000, 5000);
      tiSch->addSchedule(10000, 50000, 50000);
      tiSch->addSchedule(100000, 500000, 500000);
      TurbulenceIntensityPostprocessor vp(grid, pathname + "/ti/TI", WbWriterVtkXmlBinary::getInstance(), tiSch, comm);

      double endTime = 1000001;
      //double endTime = 1001.0;
      UbSchedulerPtr upSch(new UbScheduler(1));
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, upSch));
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


#include <iostream>
#include <string>
#include "VirtualFluids.h"


using namespace std;
void run()
{
   SPtr<Communicator> comm = MPICommunicator::getInstance();
   int myid = comm->getProcessID();

   if (myid == 0) UBLOG(logINFO, "Testcase organ pipe");

   SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

   double  deltaXcoarse = 0.064;
   const int baseLevel = 0;
   int refineLevel = 9;
   double deltaXfine = deltaXcoarse/(double)(1 << refineLevel);

   string pathOut = "e:/temp/OrganPipe";
   string pathGeo = "C:/Users/maini/Desktop/organflute/02_organf_scaled.stl";

   LBMReal rho_LB = 0.0;
   double rhoReal = 1.2041; //(kg/m3)
   double uReal = 3.0; //m/s
   double L = 0.195; //m
   double hLB = L / deltaXcoarse;
   double csReal = 343.3;
   double nuReal = 1.51e-5; //m^2/s
   LBMUnitConverter unitConverter(L, csReal, rhoReal, hLB);
   if (myid == 0) UBLOG(logINFO, unitConverter.toString());

   double nu_LB = nuReal * unitConverter.getFactorViscosityWToLb();
   double u_LB = uReal * unitConverter.getFactorVelocityWToLb();

   vector<int> blocknx = { 8, 8, 8 };

   SPtr<Grid3D> grid(new Grid3D(comm));

   if (myid == 0) UBLOG(logINFO, "Read organ pipe geometry:start");
   SPtr<GbTriFaceMesh3D> opipeGeo = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(pathGeo, "opipeGeo", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
   if (myid == 0) UBLOG(logINFO, "Read organ pipe geometry:end");
   opipeGeo->rotate(0, 270, 0);
   opipeGeo->rotate(0, 0, 90);
   opipeGeo->translate(0.0, 0.0, -(0.0952+0.0165));
   if (myid == 0) GbSystem3D::writeGeoObject(opipeGeo.get(), pathOut + "/geo/opipeGeo", WbWriterVtkXmlBinary::getInstance());

   double offs1 = opipeGeo->getX1Minimum();

   //bounding box
   double g_minX1 = offs1;
   double g_minX2 = -2.048;
   double g_minX3 = -2.048;

   double g_maxX1 = 4.096 + offs1;
   double g_maxX2 = 2.048;
   double g_maxX3 = 2.048;
   double radius_inlet = 0.006;

   if (myid == 0)
   {
      UBLOG(logINFO, "Parameters:");
      UBLOG(logINFO, "u_LB                = " << u_LB);
      UBLOG(logINFO, "rho_LB              = " << rho_LB);
      UBLOG(logINFO, "nu_LB               = " << nu_LB);
      UBLOG(logINFO, "dx coarse           = " << deltaXcoarse);
      UBLOG(logINFO, "dx fine             = " << deltaXfine);
      UBLOG(logINFO, "number of levels    = " << refineLevel + 1);
      UBLOG(logINFO, "number of processes = " << comm->getNumberOfProcesses());
      UBLOG(logINFO, "path = " << pathOut);
      UBLOG(logINFO, "Preprocess - start");
   }

   //BC adapters
   SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
   noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

   mu::Parser fct;
   fct.SetExpr("U");
   fct.DefineConst("U", u_LB);

   SPtr<BCAdapter> velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
   velBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

   SPtr<BCAdapter> outflowBCAdapter(new DensityBCAdapter(rho_LB));
   outflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithm()));

   BoundaryConditionsBlockVisitor bcVisitor;
   bcVisitor.addBC(noSlipBCAdapter);
   bcVisitor.addBC(velBCAdapter);
   bcVisitor.addBC(outflowBCAdapter);

   SPtr<BCProcessor> bcProc;
   bcProc = SPtr<BCProcessor>(new BCProcessor());

   SPtr<LBMKernel> kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());

   kernel->setBCProcessor(bcProc);

   grid->setPeriodicX1(false);
   grid->setPeriodicX2(false);
   grid->setPeriodicX3(false);
   grid->setDeltaX(deltaXcoarse);
   grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

   //generate block grid
   SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
   if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathOut + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
   GenBlocksGridVisitor genBlocks(gridCube);
  
   grid->accept(genBlocks);

   //////////////////////////////////////////////////////////////////////////
   //refinement
   SPtr<GbSphere3D> refineSphereL5(new GbSphere3D(g_minX1, 0.0, 0.0, 3.0*L));
   if (myid == 0) GbSystem3D::writeGeoObject(refineSphereL5.get(), pathOut + "/geo/refineSphereL5", WbWriterVtkXmlASCII::getInstance());

   //full pipe refine
   SPtr<GbObject3D> refineBoxL7(new GbCuboid3D(-0.1107, -0.0165, -0.0165, 0.1317, 0.0165, 0.0165));
   if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL7.get(), pathOut + "/geo/refineBoxL7", WbWriterVtkXmlASCII::getInstance());

   //refine languid and upperlip
   SPtr<GbObject3D> refineBoxL9(new GbCuboid3D(-0.1007, -0.0165, -0.0165, -0.0627, 0.0165, 0.0165));
   if (myid == 0) GbSystem3D::writeGeoObject(refineBoxL9.get(), pathOut + "/geo/refineBoxL9", WbWriterVtkXmlASCII::getInstance());
  
   SPtr<GbSphere3D> refineSphereL9(new GbSphere3D(g_minX1+75e-3, 0.0, 0.015, 24e-3));
   if (myid == 0) GbSystem3D::writeGeoObject(refineSphereL9.get(), pathOut + "/geo/refineSphereL9", WbWriterVtkXmlASCII::getInstance());

   if (refineLevel > 0)
   {
      if (myid == 0) UBLOG(logINFO, "Refinement - start");
      RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
      refineHelper.addGbObject(refineSphereL5, 5);
      refineHelper.addGbObject(refineBoxL7, 7);
      refineHelper.addGbObject(refineBoxL9, 9);
      refineHelper.addGbObject(refineSphereL9, 9);
      refineHelper.refine();
      if (myid == 0) UBLOG(logINFO, "Refinement - end");
   }
   //////////////////////////////////////////////////////////////////////////

   {
      WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
      ppblocks.process(0);
   }

   //interactors
   SPtr<Interactor3D> opipeInter = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(opipeGeo, grid, noSlipBCAdapter, Interactor3D::SOLID));
  
   //walls
   GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1-0.001 , g_minX2-0.001, g_minX3-0.001, g_maxX1+0.001, g_minX2, g_maxX3+0.001));
   if (myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathOut + "/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());
   GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1-0.001, g_maxX2, g_minX3-0.001, g_maxX1+0.001, g_maxX2+0.001, g_maxX3+0.001));
   if (myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathOut + "/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());
   GbCuboid3DPtr addWallZmin(new GbCuboid3D(g_minX1-0.001, g_minX2-0.001, g_minX3-0.001, g_maxX1+0.001, g_maxX2+0.001, g_minX3));
   if (myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathOut + "/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());
   GbCuboid3DPtr addWallZmax(new GbCuboid3D(g_minX1-0.001, g_minX2-0.001, g_maxX3, g_maxX1+0.001, g_maxX2+0.001, g_maxX3+0.001));
   if (myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathOut + "/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

   //wall interactors
   SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, velBCAdapter, Interactor3D::SOLID));
   SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, velBCAdapter, Interactor3D::SOLID));
   SPtr<D3Q27Interactor> addWallZminInt(new D3Q27Interactor(addWallZmin, grid, velBCAdapter, Interactor3D::SOLID));
   SPtr<D3Q27Interactor> addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, velBCAdapter, Interactor3D::SOLID));

   //inflow
   GbCylinder3DPtr geoInflow(new GbCylinder3D(g_minX1- deltaXcoarse, 0.0, 0.0, g_minX1 + deltaXcoarse, 0.0, 0.0, radius_inlet));
   if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathOut + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());
   //outflow
   GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - deltaXcoarse, g_minX3 - deltaXcoarse, g_maxX1 + deltaXcoarse, g_maxX2 + deltaXcoarse, g_maxX3 + deltaXcoarse));
   if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathOut + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());
   //inflow
   SPtr<D3Q27Interactor> inflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));
   //outflow
   SPtr<D3Q27Interactor> outflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowBCAdapter, Interactor3D::SOLID));

}

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   run();
}
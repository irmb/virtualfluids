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

   ////create pseudo pipe
   double offset1 = 30e-3;
   //GbCuboid3DPtr pipe(new GbCuboid3D(g_minX1, g_maxX2*0.5 - offset1, g_maxX3*0.5 - offset1, 270e-3, g_maxX2*0.5 + offset1, g_maxX3*0.5 + offset1));
   //if (myid == 0) GbSystem3D::writeGeoObject(pipe.get(), pathOut + "/geo/pipe", WbWriterVtkXmlASCII::getInstance());
   



   //////////////////////////////////////////////////////////////////////////
   //refinement
   double blockLengthX3Fine = grid->getDeltaX(refineLevel) * blocknx[2];
   double refHight = 0.002;

   //GbCuboid3DPtr refineBoxTop(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_maxX3 - refHight, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxTop.get(), pathOut + "/geo/refineBoxTop", WbWriterVtkXmlASCII::getInstance());

   ////GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+offsetMinX3+blockLengthX3Fine));
   //GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3, g_maxX1 + blockLength, g_maxX2 + blockLength, g_minX3 + refHight));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxBottom.get(), pathOut + "/geo/refineBoxBottom", WbWriterVtkXmlASCII::getInstance());

   SPtr<GbSphere3D> refineSphereL5(new GbSphere3D(g_minX1, 0.0, 0.0, 3.0*L));
   if (myid == 0) GbSystem3D::writeGeoObject(refineSphereL5.get(), pathOut + "/geo/refineSphereL5", WbWriterVtkXmlASCII::getInstance());
   
   SPtr<GbSphere3D> refineSphereL9(new GbSphere3D(g_minX1+48e-3, 0.0, 0.0, 24e-3));
   if (myid == 0) GbSystem3D::writeGeoObject(refineSphereL9.get(), pathOut + "/geo/refineSphereL9", WbWriterVtkXmlASCII::getInstance());

   if (refineLevel > 0)
   {
      if (myid == 0) UBLOG(logINFO, "Refinement - start");
      RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
      refineHelper.addGbObject(refineSphereL5, 5);
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
   SPtr<Interactor3D> opipeInter = SPtr<D3Q27TriFaceMeshInteractor>(new D3Q27TriFaceMeshInteractor(opipeGeo, grid, noSlipBCAdapter, Interactor3D::SOLID));//, Interactor3D::POINTS));

   //SPtr<D3Q27Interactor> inflowIntr = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));
}

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   run();
}
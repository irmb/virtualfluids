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
   int refineLevel = 3;
   double deltaXfine = deltaXcoarse/(double)(1 << refineLevel);

   string pathOut = "e:/temp/OrganPipe";

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

   //bounding box
   double g_minX1 = 0.0;
   double g_minX2 = 0.0;
   double g_minX3 = 0.0;

   double g_maxX1 = 4.096;
   double g_maxX2 = 4.096;
   double g_maxX3 = 4.096;

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

   BoundaryConditionsBlockVisitor bcVisitor;
   bcVisitor.addBC(noSlipBCAdapter);

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

   //create pseudo pipe
   double offset1 = 30e-3;
   GbCuboid3DPtr pipe(new GbCuboid3D(g_minX1, g_maxX2*0.5 - offset1, g_maxX3*0.5 - offset1, 270e-3, g_maxX2*0.5 + offset1, g_maxX3*0.5 + offset1));
   if (myid == 0) GbSystem3D::writeGeoObject(pipe.get(), pathOut + "/geo/pipe", WbWriterVtkXmlASCII::getInstance());


   //////////////////////////////////////////////////////////////////////////
   //refinement
   double blockLengthX3Fine = grid->getDeltaX(refineLevel) * blocknx[2];
   double refHight = 0.002;

   //GbCuboid3DPtr refineBoxTop(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_maxX3 - refHight, g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxTop.get(), pathOut + "/geo/refineBoxTop", WbWriterVtkXmlASCII::getInstance());

   ////GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1-blockLength, g_minX2-blockLength, g_minX3, g_maxX1+blockLength, g_maxX2+blockLength, g_minX3+offsetMinX3+blockLengthX3Fine));
   //GbCuboid3DPtr refineBoxBottom(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3, g_maxX1 + blockLength, g_maxX2 + blockLength, g_minX3 + refHight));
   //if (myid == 0) GbSystem3D::writeGeoObject(refineBoxBottom.get(), pathOut + "/geo/refineBoxBottom", WbWriterVtkXmlASCII::getInstance());

   SPtr<GbSphere3D> refineSphere(new GbSphere3D(g_minX1, g_maxX2*0.5 - offset1, g_maxX3*0.5 - offset1, 3.0*L));
   if (myid == 0) GbSystem3D::writeGeoObject(refineSphere.get(), pathOut + "/geo/refineSphere", WbWriterVtkXmlASCII::getInstance());

   if (refineLevel > 0)
   {
      if (myid == 0) UBLOG(logINFO, "Refinement - start");
      RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
      refineHelper.addGbObject(refineSphere, refineLevel);
      refineHelper.refine();
      if (myid == 0) UBLOG(logINFO, "Refinement - end");
   }
   //////////////////////////////////////////////////////////////////////////

   {
      WriteBlocksCoProcessor ppblocks(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathOut, WbWriterVtkXmlBinary::getInstance(), comm);
      ppblocks.process(0);
   }
}

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   run();
}
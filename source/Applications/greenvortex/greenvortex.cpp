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
#include "D3Q27BoundaryConditionAdapter.h"
#include "D3Q27PathLinePostprocessor.h"
#include "D3Q27OffsetInterpolationProcessor.h"
#include "BlocksPostprocessor.h"

using namespace std;


void run(const char *cstr)
{
   try
   {
      string pathname = "c:/temp/greenvortex/out";

      int numOfThreads = 3;

      const int blocknx1 = 5;
      const int blocknx2 = 5;
      const int blocknx3 = 5;

      const int baseLevel = 0;
      const int refineLevel = 1;

      const double blockLentghX1 = 1.0;
      const double blockLentghX2 = 1.0;
      const double blockLentghX3 = 1.0;

      const double gridOriginX1 = 0.0;
      const double gridOriginX2 = 0.0;
      const double gridOriginX3 = 0.0;

      double L1 = 5.0;
      double L2 = 5.0;
      double L3 = 5.0;

      const double dx = blockLentghX1/static_cast<double>(blocknx1);

      CommunicatorPtr comm(new MPICommunicator());

      LBMReal uLB = 0.01;
      LBMReal Re = 20.0;
      LBMReal rhoLB = 1.0;
      LBMReal l = blockLentghX2 / dx;
      LBMReal nueLB = (uLB*l)/Re;

      SimulationParametersPtr param = SimulationParameters::getInstanz();
      param->setCollisionModelType(SimulationParameters::COMPRESSIBLE);
      param->setRho(rhoLB);
      param->setVelocityX(uLB);
      param->setViscosity(nueLB);

      Grid3DPtr grid(new Grid3D());
      grid->setDeltaX(dx);
      grid->setBlockNX(blocknx1, blocknx2, blocknx3);
      grid->setPeriodicX1(true);
      grid->setPeriodicX2(false);
      grid->setPeriodicX3(false);


      GbObject3DPtr gridCube(new GbCuboid3D(0.0, 0.0, 0.0, L1, L2, L3));
      GenBlocksGridVisitor genBlocks;
      genBlocks.addGeoObject(gridCube);
      grid->accept(genBlocks);

      LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3));

      mu::Parser fctForcingX1, fctForcingX2;
      fctForcingX1.SetExpr("2.0*rho* (4.0*PI*PI/(L1/dx)/(L2/dx))*( nue*vlb*sin(x1*2.0*PI/(L1/dx))*cos(x2*2.0*PI/(L2/dx)))");
      fctForcingX2.SetExpr("-2.0*rho*(4.0*PI*PI/(L1/dx)/(L2/dx))*( nue*vlb*cos(x1*2.0*PI/(L1/dx))*sin(x2*2.0*PI/(L2/dx)))");

      fctForcingX1.DefineConst("L1"     , static_cast<double>(L1*blocknx1));
      fctForcingX1.DefineConst("L2"     , static_cast<double>(L2*blocknx2));
      fctForcingX1.DefineConst("PI"     , PI);
      fctForcingX1.DefineConst("rho"    , rhoLB);
      fctForcingX1.DefineConst("vlb"    , uLB);

      fctForcingX2.DefineConst("L1"     , static_cast<double>(L1*blocknx1));
      fctForcingX2.DefineConst("L2"     , static_cast<double>(L2*blocknx2));
      fctForcingX2.DefineConst("PI"     , PI);
      fctForcingX2.DefineConst("rho"    , rhoLB);
      fctForcingX2.DefineConst("vlb"    , uLB);

      kernel->setForcingX1(fctForcingX1);
      kernel->setForcingX2(fctForcingX2);

      BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
      kernel->setBCProcessor(bcProc);

      SetKernelBlockVisitor kernelVisitor(kernel, nueLB);
      grid->accept(kernelVisitor);

      D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
      D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
      grid->accept( setConnsVisitor );

      PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
      grid->accept(pqPartVisitor);

      D3Q27ETInitDistributionsBlockVisitor initVisitor(1.0);

      mu::Parser fct,fct2,fct3;
      fct.SetExpr(" vLB*sin( ( (x1)*2.0*PI/ L1))*cos( (x2)*2.0*PI/L2)");
      fct.DefineConst("L1"     , L1);
      fct.DefineConst("L2"     , L2);
      fct.DefineConst("vLB"  , uLB);
      fct.DefineConst("PI"  , PI);
      initVisitor.setVx1(fct);

      fct2.SetExpr(" -vLB*cos( ( (x1)*2.0*PI/ L1))*sin( (x2)*2.0*PI/L2)");
      fct2.DefineConst("L1"     , L1);
      fct2.DefineConst("L2"     , L2);
      fct2.DefineConst("vLB"  , uLB           );
      fct2.DefineConst("PI"  , PI);
      initVisitor.setVx2(fct2);

      initVisitor.setVx3(0.0);

      fct3.SetExpr(" 1.0+(vLB*vLB)*3.0/4.0*(cos((x1)*4.0*PI/L1)+cos((x2)*4.0*PI/L2))");
      fct3.DefineConst("L1"     , L1);
      fct3.DefineConst("L2"     , L2);
      fct3.DefineConst("vLB"  , uLB           );
      fct3.DefineConst("PI"  , PI);
      initVisitor.setRho(fct3);

      grid->accept(initVisitor);

      BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));
      ppblocks->update(0);
      ppblocks.reset();

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());
      {
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessor ppgeo(grid, geoSch, pathname + "/nodes_geo", WbWriterVtkXmlBinary::getInstance(), conv, comm, true);
         grid->doPostProcess(0);
      }
      double outTime = 1000.0;
      UbSchedulerPtr visSch(new UbScheduler(outTime));
      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname + "/nodes", WbWriterVtkXmlBinary::getInstance(), conv, comm);

      //////////////////////////////////////////////////////////////////////////
      //PathLine
      UbSchedulerPtr plSch(new UbScheduler(1000, 1000));
      D3Q27PathLinePostprocessor pathLine(grid, pathname + "/pathLine", WbWriterVtkXmlASCII::getInstance(), conv, plSch, comm, 4.2, 4.2, 4.2, nueLB, iProcessor);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //Simulation
      //////////////////////////////////////////////////////////////////////////
      double endTime = 10000.0;
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
      UBLOG(logINFO,"Simulation-start");
      calculation->calculate();

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


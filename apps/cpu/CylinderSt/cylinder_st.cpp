#include <iostream>
#include <string>

#include "geometry3d/CoordinateTransformation3D.h"
#include "Grid3D.h"
#include "GenBlocksGridVisitor.h"
#include "geometry3d/GbSystem3D.h"
#include "geometry3d/GbCuboid3D.h"
#include "geometry3d/GbCylinder3D.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "OverlapBlockVisitor.h"
#include "RefineInterGbObjectsVisitor.h"
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
#include "D3Q27PressureDifferencePostprocessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "RestartPostprocessor.h"
#include "SolidBlocksHelper.h"
#include "NUPSCounterPostprocessor.h"
#include "BlocksPostprocessor.h"
#include "LBMKernelETD3Q27BGK.h"
#include "EmergencyExitPostprocessor.h"
#include "D3Q27ForcesPostprocessor.h"

using namespace std;


void run(const char *cstr)
{
   try
   {
      string machine = QUOTEME(CAB_MACHINE);
      string pathname; 
      int numOfThreads = 1;
      double availMem = 0;

      CommunicatorPtr comm(new MPICommunicator());
      int myid = comm->getProcessID();

      if(machine == "BOMBADIL") 
      {
         pathname = "c:/temp/cylinder_st";
         numOfThreads = 3;
         availMem = 3.0e9;
      }
      else if(machine == "M01" || machine == "M02")      
      {
         pathname = "/work/koskuche/scratch/cylinder_st";
         numOfThreads = 8;
         availMem = 12.0e9;

         if(myid ==0)
         {
           stringstream logFilename;
           logFilename <<  pathname + "/logfile.txt";
           UbLog::output_policy::setStream(logFilename.str());
         }
      }
      else throw UbException(UB_EXARGS, "unknown CAB_MACHINE");

      double dx = 0.00207051;

      double L1 = 2.5;
      double L2, L3, H;
      L2 = L3 = H = 0.41;

      LBMReal radius = 0.05;
      LBMReal uLB = 0.05;
      LBMReal Re = 1000.0;
      LBMReal rhoLB = 1.0;
      LBMReal l = L2 / dx;
      //LBMReal nueLB = (uLB*l)/Re;
      LBMReal nueLB = (((4.0/9.0)*uLB)*2.0*(radius/dx))/Re;

      LBMUnitConverterPtr conv = LBMUnitConverterPtr(new LBMUnitConverter());

      const int baseLevel = 0;
      const int refineLevel = 2;

      //obstacle
      GbObject3DPtr cylinder(new GbCylinder3D(0.5, 0.2, -0.1, 0.5, 0.2, L3+0.1, radius));
      GbSystem3D::writeGeoObject(cylinder.get(),pathname + "/geo/cylinder", WbWriterVtkXmlBinary::getInstance());

      D3Q27InteractorPtr cylinderInt;

      //bounding box
      double d_minX1 = 0.0;
      double d_minX2 = 0.0;
      double d_minX3 = 0.0;

      double d_maxX1 = L1;
      double d_maxX2 = L2;
      double d_maxX3 = L3;

      double offs = dx;

      //double g_minX1 = d_minX1-offs-0.499999*dx;
      double g_minX1 = d_minX1-offs;
      double g_minX2 = d_minX2-offs;
      double g_minX3 = d_minX3-offs;

      double g_maxX1 = d_maxX1+offs;
      double g_maxX2 = d_maxX2+offs;
      double g_maxX3 = d_maxX3+offs;

      GbObject3DPtr gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));

      const int blocknx1 = 10;
      const int blocknx2 = 10;
      const int blocknx3 = 10;
      
      dx = (0.41+2*dx)/(20.0*(int)blocknx2);

      double blockLength = blocknx1*dx;

      //refinement area
      //double rf = cylinder->getLengthX1()/5;
      //GbObject3DPtr refineCube(new  GbCuboid3D(cylinder->getX1Minimum()-rf, cylinder->getX2Minimum()-rf, cylinder->getX3Minimum(), 
      //   cylinder->getX1Maximum()+6.0*rf, cylinder->getX2Maximum()+rf, cylinder->getX3Maximum()));
      GbObject3DPtr refineCube(new  GbCuboid3D(g_minX1 + 20.0*blockLength, g_minX2 + 6.0*blockLength, cylinder->getX3Minimum(), 
                                               g_minX1 + 33.0*blockLength, g_maxX2 - 6.0*blockLength, cylinder->getX3Maximum()));

      Grid3DPtr grid(new Grid3D());

      UbSchedulerPtr rSch(new UbScheduler(100000, 100000));
      RestartPostprocessorPtr rp(new RestartPostprocessor(grid, rSch, comm, pathname+"/checkpoints", RestartPostprocessor::BINARY));

      //UbSchedulerPtr emSch(new UbScheduler(1000, 1000));
      //EmergencyExitPostprocessor em(grid, emSch, pathname+"/checkpoints/emex.txt", rp, comm);

      std::string opt;

      if(cstr!= NULL)
         opt = std::string(cstr);

      if/*(cstr== NULL)*/(cstr!= NULL)
      {
         opt = std::string(cstr);

         if(myid==0) UBLOG(logINFO,"Restart step: " << opt);

         grid = rp->restart(UbSystem::stringTo<int>(opt));
         rp->reconnect();

         //cylinderInt = 

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
         if(myid ==0)
         {
            UBLOG(logINFO,"L = " << l );
            UBLOG(logINFO,"v = " << uLB );
            UBLOG(logINFO,"rho = " << rhoLB );
            UBLOG(logINFO,"nue = " << nueLB );
            UBLOG(logINFO,"Re = " << Re );
            UBLOG(logINFO,"dx = " << dx );
            UBLOG(logINFO,"Preprozess - start");
         }
        
         grid->setDeltaX(dx);
         grid->setBlockNX(blocknx1, blocknx2, blocknx3);
         
         // UbTupleDouble6 bouningBox(gridCube->getX1Minimum(),gridCube->getX2Minimum(),gridCube->getX3Minimum(),
                                   // gridCube->getX1Maximum(),gridCube->getX2Maximum(),gridCube->getX3Maximum());
         // UbTupleInt3 blockNx(blocknx1, blocknx2, blocknx3);
         // UbTupleInt3 gridNx(8, 16, 16);
         // grid = Grid3DPtr(new Grid3D(bouningBox, blockNx, gridNx));
         
         if(myid ==0) GbSystem3D::writeGeoObject(gridCube.get(),pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());
         if(myid ==0) GbSystem3D::writeGeoObject(refineCube.get(),pathname + "/geo/refineCube", WbWriterVtkXmlBinary::getInstance());
      
         GenBlocksGridVisitor genBlocks;
         genBlocks.addGeoObject(gridCube);
         grid->accept(genBlocks);

         //walls
         GbCuboid3DPtr addWallYmin (new GbCuboid3D(d_minX1-blockLength, d_minX2-blockLength, d_minX3-blockLength, d_maxX1+blockLength, d_minX2, d_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmin.get(), pathname+"/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());
      
         GbCuboid3DPtr addWallZmin (new GbCuboid3D(d_minX1-blockLength, d_minX2-blockLength, d_minX3-blockLength, d_maxX1+blockLength, d_maxX2+blockLength, d_minX3));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmin.get(), pathname+"/geo/addWallZmin", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallYmax (new GbCuboid3D(d_minX1-blockLength, d_maxX2, d_minX3-blockLength, d_maxX1+blockLength, d_maxX2+2.0*blockLength, d_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallYmax.get(), pathname+"/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

         GbCuboid3DPtr addWallZmax (new GbCuboid3D(d_minX1-blockLength, d_minX2-blockLength, d_maxX3, d_maxX1+blockLength, d_maxX2+blockLength, d_maxX3+2.0*blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(addWallZmax.get(), pathname+"/geo/addWallZmax", WbWriterVtkXmlASCII::getInstance());

         //inflow
         GbCuboid3DPtr geoInflow (new GbCuboid3D(d_minX1-blockLength, d_minX2-blockLength, d_minX3-blockLength, d_minX1, d_maxX2+blockLength, d_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname+"/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

         //outflow
         GbCuboid3DPtr geoOutflow (new GbCuboid3D(d_maxX1, d_minX2-blockLength, d_minX3-blockLength, d_maxX1+2.0*blockLength, d_maxX2+blockLength, d_maxX3+blockLength));
         if(myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname+"/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

         BlocksPostprocessorPtr ppblocks(new BlocksPostprocessor(grid, UbSchedulerPtr(new UbScheduler(1)), pathname + "/grid/blocks", WbWriterVtkXmlBinary::getInstance(), comm));

         if (refineLevel > 0)
         {
            if(myid == 0) UBLOG(logINFO,"Refinement - start");	
            RefineCrossAndInsideGbObjectBlockVisitor refVisitor(refineCube, baseLevel, refineLevel-1);
            grid->accept(refVisitor);

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
         }

         MetisPartitioningGridVisitor metisVisitor(numOfThreads, D3Q27System::B, comm, false);
         grid->accept( metisVisitor );

         SolidBlocksHelper sd(grid, comm);

         int bbOption = 1; //0=simple Bounce Back, 1=quadr. BB
         D3Q27BoundaryConditionAdapterPtr bcObst(new D3Q27NoSlipBCAdapter(bbOption));
         cylinderInt = D3Q27InteractorPtr ( new D3Q27Interactor(cylinder, grid, bcObst,Interactor3D::SOLID));

         //walls
         D3Q27InteractorPtr addWallYminInt(new D3Q27Interactor(addWallYmin, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZminInt(new D3Q27Interactor(addWallZmin, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, bcObst,Interactor3D::SOLID));
         D3Q27InteractorPtr addWallZmaxInt(new D3Q27Interactor(addWallZmax, grid, bcObst,Interactor3D::SOLID));

         mu::Parser fct;
         fct.SetExpr("16*U*x2*x3*(H-x2)*(H-x3)/H^4");
         fct.DefineConst("U", uLB);
         fct.DefineConst("H", H);

         //inflow
         D3Q27BoundaryConditionAdapterPtr velBCAdapter(new D3Q27VelocityBCAdapter (true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST));
         velBCAdapter->setSecondaryBcOption(2);
         D3Q27InteractorPtr inflowInt  = D3Q27InteractorPtr( new D3Q27Interactor(geoInflow, grid, velBCAdapter, Interactor3D::SOLID));

         //outflow
         D3Q27BoundaryConditionAdapterPtr denBCAdapter(new D3Q27DensityBCAdapter(rhoLB));
         D3Q27InteractorPtr outflowInt = D3Q27InteractorPtr( new D3Q27Interactor(geoOutflow, grid, denBCAdapter,Interactor3D::SOLID));

         sd.addInteractor(cylinderInt);
         sd.addInteractor(addWallYminInt);
         sd.addInteractor(addWallZminInt);
         sd.addInteractor(addWallYmaxInt);
         sd.addInteractor(addWallZmaxInt);
         sd.addInteractor(inflowInt);
         sd.addInteractor(outflowInt);
         
         sd.deleteSolidBlocks();
         
         grid->accept( metisVisitor );


         ppblocks->update(0);
         ppblocks.reset();

         //set connectors
         D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
         D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
         grid->accept( setConnsVisitor );

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

         LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(blocknx1, blocknx2, blocknx3));
         //LBMKernel3DPtr kernel(new LBMKernelETD3Q27BGK(blocknx1, blocknx2, blocknx3, true));

         BCProcessorPtr bcProc(new D3Q27ETBCProcessor());
         kernel->setBCProcessor(bcProc);

         SetKernelBlockVisitor kernelVisitor(kernel, nueLB, availMem, needMem);
         grid->accept(kernelVisitor);

         if (refineLevel > 0)
         {
            D3Q27SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
         }

         //walls
         grid->addAndInitInteractor(addWallYminInt);
         grid->addAndInitInteractor(addWallZminInt);
         grid->addAndInitInteractor(addWallYmaxInt);
         grid->addAndInitInteractor(addWallZmaxInt);

         //obstacle
         grid->addAndInitInteractor(cylinderInt);

         //inflow
         grid->addAndInitInteractor(inflowInt);

         //outflow
         grid->addAndInitInteractor(outflowInt);

         //domain decomposition
         PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
         grid->accept(pqPartVisitor);

         //initialization of distributions
         D3Q27ETInitDistributionsBlockVisitor initVisitor(1.0);
         initVisitor.setVx1(fct);
         grid->accept(initVisitor);

         //Postrozess
         UbSchedulerPtr geoSch(new UbScheduler(1));
         D3Q27MacroscopicQuantitiesPostprocessorPtr ppgeo(
            new D3Q27MacroscopicQuantitiesPostprocessor(grid, geoSch, pathname + "/grid/nodes", WbWriterVtkXmlBinary::getInstance(), conv, comm, true));
         ppgeo->update(0);
         ppgeo.reset();
         
         if(myid == 0) UBLOG(logINFO,"Preprozess - end"); 
      }

      double outTime = 50000.0;
      UbSchedulerPtr visSch(new UbScheduler(outTime));
      visSch->addSchedule(1000, 1000, 10000);
      visSch->addSchedule(10000, 10000, 50000);
      visSch->addSchedule(1, 1, 10000);

      D3Q27MacroscopicQuantitiesPostprocessor pp(grid, visSch, pathname + "/steps/step", WbWriterVtkXmlBinary::getInstance(), conv, comm);

      double fdx = grid->getDeltaX(grid->getFinestInitializedLevel());
      double point1[3] = {0.45, 0.20, 0.205};
      double point2[3] = {0.55, 0.20, 0.205};

      D3Q27IntegrateValuesHelperPtr h1(new D3Q27IntegrateValuesHelper(grid, comm, 
         point1[0]-1.0*fdx, point1[1]-1.0*fdx, point1[2]-1.0*fdx, 
         point1[0], point1[1], point1[2]));
      if(myid ==0) GbSystem3D::writeGeoObject(h1->getBoundingBox().get(),pathname + "/geo/iv1", WbWriterVtkXmlBinary::getInstance());
      D3Q27IntegrateValuesHelperPtr h2(new D3Q27IntegrateValuesHelper(grid, comm, 
         point2[0], point2[1]-1.0*fdx, point2[2]-1.0*fdx, 
         point2[0]+1.0*fdx, point2[1], point2[2]));
      if(myid ==0) GbSystem3D::writeGeoObject(h2->getBoundingBox().get(),pathname + "/geo/iv2", WbWriterVtkXmlBinary::getInstance());
      D3Q27PressureDifferencePostprocessor rhopp(grid, visSch, pathname + "/results/rho_diff.txt", h1, h2, conv, comm);

      UbSchedulerPtr nupsSch(new UbScheduler(10, 10, 10));
      double area = (radius*radius*H)/fdx;
      double v    = 4.0*uLB/9.0;
      D3Q27ForcesPostprocessor fp(grid, visSch, pathname + "/results/forces.txt", cylinderInt, comm, rhoLB, v, area, D3Q27ForcesPostprocessor::X, D3Q27ForcesPostprocessor::Y);

      NUPSCounterPostprocessor npr(grid, nupsSch, pathname + "/results/nups.txt", comm);

      double endTime = 1000001.0;
      CalculationManagerPtr calculation(new CalculationManager(grid, numOfThreads, endTime, visSch));
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


//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file LidDrivenCavity.cpp
//! \ingroup Applications
//! \author Konstantin Kutscher
//=======================================================================================

#include <string>

#include "VirtualFluids.h"

using namespace std;

int main(int  /*argc*/, char*  /*argv*/[])
{
   try
   {
      //////////////////////////////////////////////////////////////////////////
      // Simulation parameters
      //////////////////////////////////////////////////////////////////////////

      // set your output path here
      string path = "./output";

      const double L = 1.0;
      const double Re = 1000.0;
      const double velocity = 1.0;
      const double dt = 0.5e-3;
      const unsigned int nx = 64;

      const double timeStepOut = 1000;
      const double timeStepEnd = 25000;

      // Number of OpenMP threads
      int numOfThreads = 4;

      //////////////////////////////////////////////////////////////////////////

      double dx = L / double(nx);
      const double velocityLB = velocity * dt / dx; // LB units
      const double u = velocityLB / sqrt(2.0); // LB units
      const double viscosityLB = nx * velocityLB / Re; // LB unit

      //////////////////////////////////////////////////////////////////////////
      // create grid
      //////////////////////////////////////////////////////////////////////////
      // bounding box
      double g_minX1 = -0.5;
      double g_minX2 = -0.5;
      double g_minX3 = -0.5;

      double g_maxX1 = 0.5;
      double g_maxX2 = 0.5;
      double g_maxX3 = 0.5;

      // NullCommunicator is a place-holder for interprocess communication
      SPtr<vf::mpi::Communicator> comm = NullCommunicator::getInstance();
      // new grid object
      SPtr<Grid3D> grid(new Grid3D(comm));
      // set grid spacing
      grid->setDeltaX(dx);
      // set block size for three dimensions
      int blockSize = nx / 2;
      grid->setBlockNX(blockSize,blockSize,blockSize);
      
      // Create simulation bounding box
      SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
      GbSystem3D::writeGeoObject(gridCube.get(), path + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

      UBLOG(logINFO, "Lid Driven Cavity:");
      UBLOG(logINFO, "Domain size = " << nx << " x "<< nx << " x "<< nx);
      UBLOG(logINFO, "Block size = " << blockSize << " x "<< blockSize << " x "<< blockSize);
      UBLOG(logINFO, "velocity    = " << velocity << " m/s");
      UBLOG(logINFO, "velocityLB  = " << velocityLB);
      UBLOG(logINFO, "viscosityLB = " << viscosityLB);
      UBLOG(logINFO, "u  = " << u);
      UBLOG(logINFO, "Re = " << Re);
      UBLOG(logINFO, "dx = " << dx);
      UBLOG(logINFO, "dt = " << dt);
      UBLOG(logINFO, "Preprocess - start");

      // Generate block grid
      GenBlocksGridVisitor genBlocks(gridCube);
      grid->accept(genBlocks);

      // Write block grid to VTK-file
      auto ppblocks = std::make_shared<WriteBlocksCoProcessor>(grid, SPtr<UbScheduler>(new UbScheduler(1)), path, WbWriterVtkXmlBinary::getInstance(), comm);
      ppblocks->process(0);
      ppblocks.reset();

      // Create LBM kernel
      auto kernel = std::make_shared<LBMKernelETD3Q27BGK>();
      // auto kernel = std::make_shared<CumulantK17LBMKernel>();

      //////////////////////////////////////////////////////////////////////////
      // Create boundary conditions (BC)
      //////////////////////////////////////////////////////////////////////////     
      // Create no-slip BC
      auto noSlipBCAdapter = std::make_shared<NoSlipBCAdapter>();
      noSlipBCAdapter->setBcAlgorithm(std::make_shared<NoSlipBCAlgorithm>());
      
      // Velocity BC
      mu::Parser fct;
      fct.SetExpr("u");
      fct.DefineConst("u", u);
      // Set the same velocity in x and y-direction
      auto velBCAdapter = std::make_shared<VelocityBCAdapter>(true, true, false, fct, 0, BCFunction::INFCONST);
      velBCAdapter->setBcAlgorithm(std::make_shared<VelocityBCAlgorithm>());

      // Add velocity boundary condition to visitor. No-slip boundary   
      BoundaryConditionsBlockVisitor bcVisitor;
      bcVisitor.addBC(velBCAdapter);

      // Create boundary conditions processor
      SPtr<BCProcessor> bcProc;
      bcProc = std::make_shared<BCProcessor>();
      kernel->setBCProcessor(bcProc);

      // Create boundary conditions geometry
      GbCuboid3DPtr wallXmin(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_minX1, g_maxX2 + dx, g_maxX3));
      GbSystem3D::writeGeoObject(wallXmin.get(), path + "/geo/wallXmin", WbWriterVtkXmlASCII::getInstance());
      GbCuboid3DPtr wallXmax(new GbCuboid3D(g_maxX1, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_maxX3));
      GbSystem3D::writeGeoObject(wallXmax.get(), path + "/geo/wallXmax", WbWriterVtkXmlASCII::getInstance());
      GbCuboid3DPtr wallYmin(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_minX2, g_maxX3));
      GbSystem3D::writeGeoObject(wallYmin.get(), path + "/geo/wallYmin", WbWriterVtkXmlASCII::getInstance());
      GbCuboid3DPtr wallYmax(new GbCuboid3D(g_minX1 - dx, g_maxX2, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_maxX3));
      GbSystem3D::writeGeoObject(wallYmax.get(), path + "/geo/wallYmax", WbWriterVtkXmlASCII::getInstance());
      GbCuboid3DPtr wallZmin(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_minX3));
      GbSystem3D::writeGeoObject(wallZmin.get(), path + "/geo/wallZmin", WbWriterVtkXmlASCII::getInstance());
      GbCuboid3DPtr wallZmax(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_maxX3, g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx));
      GbSystem3D::writeGeoObject(wallZmax.get(), path + "/geo/wallZmax", WbWriterVtkXmlASCII::getInstance());

      // Add boundary conditions to grid generator
      SPtr<D3Q27Interactor> wallXminInt(new D3Q27Interactor(wallXmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
      SPtr<D3Q27Interactor> wallXmaxInt(new D3Q27Interactor(wallXmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
      SPtr<D3Q27Interactor> wallYminInt(new D3Q27Interactor(wallYmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
      SPtr<D3Q27Interactor> wallYmaxInt(new D3Q27Interactor(wallYmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
      SPtr<D3Q27Interactor> wallZminInt(new D3Q27Interactor(wallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
      SPtr<D3Q27Interactor> wallZmaxInt(new D3Q27Interactor(wallZmax, grid, velBCAdapter, Interactor3D::SOLID));

      InteractorsHelper intHelper(grid);
      intHelper.addInteractor(wallZmaxInt);
      intHelper.addInteractor(wallXminInt);
      intHelper.addInteractor(wallXmaxInt);
      intHelper.addInteractor(wallZminInt);
      intHelper.addInteractor(wallYminInt);
      intHelper.addInteractor(wallYmaxInt);

      intHelper.selectBlocks();

      // Generate grid
      SetKernelBlockVisitor kernelVisitor(kernel, viscosityLB);
      grid->accept(kernelVisitor);

      intHelper.setBC();

      // Initialization of distributions
      InitDistributionsBlockVisitor initVisitor;
      grid->accept(initVisitor);

      // Set connectivity between blocks
      SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, viscosityLB);
      grid->accept(setConnsVisitor);

      // Create lists of boundary nodes
      grid->accept(bcVisitor);

      // Write grid with boundary conditions information to VTK-file
      SPtr<UbScheduler> geoSch(new UbScheduler(1));
      WriteBoundaryConditionsCoProcessor ppgeo(grid, geoSch, path, WbWriterVtkXmlBinary::getInstance(), comm);
      ppgeo.process(0);

      UBLOG(logINFO, "Preprocess - end");
      
      UBLOG(logINFO, "Total Physical Memory (RAM): " << Utilities::getTotalPhysMem()/1e9 << " GB");
      UBLOG(logINFO, "Physical Memory currently used: " << Utilities::getPhysMemUsed()/1e9 << " GB");
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe()/1e9 << " GB");

      // Create coprocessor object for writing macroscopic quantities to VTK-file
      SPtr<UbScheduler> visSch(new UbScheduler(timeStepOut));
      SPtr<CoProcessor> mqCoProcessor(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, path, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter(L, velocity, 1.0, nx, velocityLB)), comm));
      mqCoProcessor->process(0);

      // Create coprocessor object for writing NUPS
      SPtr<UbScheduler> nupsSch(new UbScheduler(100, 100));
      SPtr<CoProcessor> nupsCoProcessor(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

      // OpenMP threads control
#ifdef _OPENMP
      omp_set_num_threads(numOfThreads);
#endif
      // Create simulation
      SPtr<Calculator> calculator(new BasicCalculator(grid, visSch, (int)timeStepEnd));
      // Add coprocessors objects to simulation
      calculator->addCoProcessor(nupsCoProcessor);
      calculator->addCoProcessor(mqCoProcessor);
    
      //////////////////////////////////////////////////////////////////////////
      // Run simulation
      //////////////////////////////////////////////////////////////////////////

      UBLOG(logINFO, "Simulation-start");
      calculator->calculate();
      UBLOG(logINFO, "Simulation-end");
   }
   catch (std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch (std::string& s)
   {
      cerr << s << endl;
   }
   catch (...)
   {
      cerr << "unknown exception" << endl;
   }
}


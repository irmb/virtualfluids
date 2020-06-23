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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file Simulation.cpp
//! \ingroup LBM
//! \author Martin Schoenherr
//=======================================================================================
#include "Simulation.h"
#include "LBM/LB.h"
#include "Parameter/Parameter.h"
#include "GPU/GPU_Interface.h"
#include "GPU/devCheck.h"
#include "GPU/CudaMemoryManager.h"	
#include "Init/InitLattice.h"
#include "DataStructureInitializer/GridProvider.h"
#include "Output/DataWriter.h"

#include "Core/Timer/Timer.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <stdio.h>
#include <memory>
//////////////////////////////////////////////////////////////////////////

Simulation::Simulation()
{}

Simulation::~Simulation()
{}

void Simulation::init(SPtr<Parameter> para, SPtr<GridProvider> gridProvider, SPtr<DataWriter> dataWriter, SPtr<CudaMemoryManager> cudaManager)
{
	devCheck(0);

	this->dataWriter = dataWriter;
	this->gridProvider = gridProvider;
	this->cudaManager = cudaManager;
	this->para = para;
	
	para->initParameter();
	para->setRe(para->getVelocityLB() * (real)1.0 / para->getViscosityLB());

	gridProvider->allocAndCopyForcing();

	//////////////////////////////////////////////////////////////////////////
	// create and use log file
	output.setName(para->getPathAndFilename() + ".log");
	output.setConsoleOut(true);
	output.clearLogFile();

	output << "Re:         " << para->getRe()              << "\n";
	output << "vis_ratio:  " << para->getViscosityRatio()  << "\n";
	output << "u0_ratio:   " << para->getVelocityRatio()   << "\n";
	output << "rho_ratio:  " << para->getDensityRatio()    << "\n";
	output << "velocityLB: " << para->getVelocityLB()      << "\n";
	output << "viscosityLB:" << para->getViscosityLB()     << "\n";

	/////////////////////////////////////////////////////////////////////////
	// set the used device memory to 0 before starting the calculation 
	cudaManager->setMemsizeGPU(0, true);

	//////////////////////////////////////////////////////////////////////////
	// allocate the memory for several arrays 
	gridProvider->allocArrays_CoordNeighborGeo();
	gridProvider->allocArrays_BoundaryValues();
	gridProvider->allocArrays_BoundaryQs();

	//////////////////////////////////////////////////////////////////////////
	// initialize the grid
	output << "init lattice...";
	initLattice(para);
	output << "done.\n";

	//////////////////////////////////////////////////////////////////////////
	// print initialized grid
	output << "Print files Init...";
	dataWriter->writeInit(para, cudaManager);
	output << "done.\n";

	//////////////////////////////////////////////////////////////////////////
	// print the amount of used device memory
	output << "used device memory: " << cudaManager->getMemsizeGPU() / 1000000.0 << " MB\n";
}

void Simulation::run()
{
	uint t;
	//////////////////////////////////////////////////////////////////////////
	// Timer
	auto timer = Timer::makeStart();
	real timeComplete = 0.0;
	// MNUPS(million node updates per second)
	output << "Processing time (s) \t NUPS * 10^6\n";

	////////////////////////////////////////////////////////////////////////////////
	// Time loop
	for (t = para->getTimestepStart(); t <= para->getTimestepEnd(); t++)
	{
		////////////////////////////////////////////////////////////////////////////////
		// LBM Kernel
		CumulantK17LBMDeviceKernel(
			para->getParD()->numberofthreads,
			para->getParD()->omega,
			para->getParD()->typeOfGridNode,
			para->getParD()->neighborX,
			para->getParD()->neighborY,
			para->getParD()->neighborZ,
			para->getParD()->distributions.f[0],
			para->getParD()->numberOfNodes,
			para->getParD()->forcing,
			para->getParD()->isEvenTimestep);

		////////////////////////////////////////////////////////////////////////////////
		//velocity boundary condition
		QVelDevicePlainBB27(
			para->getParD()->numberofthreads,
			para->getParD()->inflowBC.Vx,
			para->getParD()->inflowBC.Vy,
			para->getParD()->inflowBC.Vz,
			para->getParD()->distributions.f[0],
			para->getParD()->inflowBC.k,
			para->getParD()->inflowBC.q27[0],
			para->getParD()->numberOfInflowBCnodes,
			para->getParD()->inflowBC.kArray,
			para->getParD()->neighborX,
			para->getParD()->neighborY,
			para->getParD()->neighborZ,
			para->getParD()->numberOfNodes,
			para->getParD()->isEvenTimestep);

		////////////////////////////////////////////////////////////////////////////////
		if (para->getParD()->isEvenTimestep)	para->getParD()->isEvenTimestep = false;
		else									para->getParD()->isEvenTimestep = true;
		////////////////////////////////////////////////////////////////////////////////
		// File IO and calculation of MNUPS(million node updates per second)
		if (para->getTimestepOut() > 0 && t%para->getTimestepOut() == 0 && t > para->getTimestepStartOut())
		{
			checkCudaErrors(cudaDeviceSynchronize());
			//////////////////////////////////////////////////////////////////////////
			// Timer
			timer->end();
			real timeInterval = timer->getTimeInSeconds();
			timeComplete += timeInterval;
			real fnups = ((real)(t - para->getTimestepStart()) * para->getParH()->numberOfNodes) / (timeComplete*1.0E6);
			output << timeInterval << " / " << timeComplete << " \t " << fnups << "\n";
			//////////////////////////////////////////////////////////////////////////
			//IO
			if (para->getPrintFiles())
			{
				output << "File IO for t=" << t << "...";
				////////////////////////////////////////////////////////////////////////////////
				CalcMacCompSP27(
					para->getParD()->velocityX,
					para->getParD()->velocityY,
					para->getParD()->velocityZ,
					para->getParD()->rho,
					para->getParD()->pressure,
					para->getParD()->typeOfGridNode,
					para->getParD()->neighborX,
					para->getParD()->neighborY,
					para->getParD()->neighborZ,
					para->getParD()->numberOfNodes,
					para->getParD()->numberofthreads,
					para->getParD()->distributions.f[0],
					para->getParD()->isEvenTimestep);
				////////////////////////////////////////////////////////////////////////
				cudaManager->cudaCopyDataToHost();
				////////////////////////////////////////////////////////////////////////
				dataWriter->writeTimestep(para, t);
				////////////////////////////////////////////////////////////////////////
				output << "done.\n";
			}
			timer->start();
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Timer
	timer->end();
	real timeInterval = timer->getTimeInSeconds();
	timeComplete += timeInterval;
	real fnups = ((real)(t - para->getTimestepStart()) * para->getParH()->numberOfNodes) / (timeComplete*1.0E6);
	output << "Processing time: " << timeComplete << "(ms)\n";
	output << "NUPS: " << fnups << " * 10^6\n";
}

void Simulation::free()
{
	//CudaFreeHostMemory
	cudaManager->cudaFreeCoord();
	cudaManager->cudaFreeSP();
	cudaManager->cudaFreeNeighborWSB();
	cudaManager->cudaFreeVeloBC();
	cudaManager->cudaFreeForcing();

	para->~Parameter();
	gridProvider->~GridProvider();
	dataWriter->~DataWriter();
}
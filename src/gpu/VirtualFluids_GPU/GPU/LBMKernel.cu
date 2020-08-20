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
//! \file LBMKernel.cu
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
#include <helper_cuda.h>

#include "LBM/LB.h"
#include "GPU/GPU_Kernels.cuh"
//////////////////////////////////////////////////////////////////////////
extern "C" void CumulantK17LBMDeviceKernel(
	uint numberOfThreads,
	real omega,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	int size_Mat,
	real* forces,
	bool isEvenTimestep)
{
	int Grid = (size_Mat / numberOfThreads) + 1;
	dim3 grid(Grid, 1, 1);
	dim3 threads(numberOfThreads, 1, 1);

	Cumulant_K17_LBM_Device_Kernel <<< grid, threads >>> (
		omega,
		typeOfGridNode,
		neighborX,
		neighborY,
		neighborZ,
		distributions,
		size_Mat,
		forces,
		isEvenTimestep);
	getLastCudaError("Cumulant_K17_AA2016_chim_Comp_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
extern "C" void LB_Init(
	uint numberOfThreads,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint* typeOfGridNode,
	real* rho,
	real* velocityX,
	real* velocityY,
	real* velocityZ,
	uint size_Mat,
	real* distributions,
	bool isEvenTimestep)
{
	int Grid = (size_Mat / numberOfThreads) + 1;
	dim3 grid(Grid, 1, 1);
	dim3 threads(numberOfThreads, 1, 1);

	LBInit <<< grid, threads >>> (
		neighborX,
		neighborY,
		neighborZ,
		typeOfGridNode,
		rho,
		velocityX,
		velocityY,
		velocityZ,
		size_Mat,
		distributions,
		isEvenTimestep);
	getLastCudaError("LBInit execution failed");
}
//////////////////////////////////////////////////////////////////////////
extern "C" void CalcMacCompSP27(
	real* velocityX,
	real* velocityY,
	real* velocityZ,
	real* rho,
	real* pressure,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint size_Mat,
	uint numberOfThreads,
	real* distributions,
	bool isEvenTimestep)
{
	int Grid = (size_Mat / numberOfThreads) + 1;
	dim3 grid(Grid, 1, 1);
	dim3 threads(numberOfThreads, 1, 1);

	LBCalcMacCompSP27 <<< grid, threads >>> (
		velocityX,
		velocityY,
		velocityZ,
		rho,
		pressure,
		typeOfGridNode,
		neighborX,
		neighborY,
		neighborZ,
		size_Mat,
		distributions,
		isEvenTimestep);
	getLastCudaError("LBCalcMacSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
extern "C" void QVelDevicePlainBB27(
	unsigned int numberOfThreads,
	real* velocityX,
	real* velocityY,
	real* velocityZ,
	real* distributions,
	int* k_Q,
	real* QQ,
	uint sizeQ,
	int kQ,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint size_Mat,
	bool isEvenTimestep)
{
	int Grid = (kQ / numberOfThreads) + 1;
	dim3 gridQ(Grid, 1, 1);
	dim3 threads(numberOfThreads, 1, 1);

	QVelDevPlainBB27 <<< gridQ, threads >>> (
		velocityX,
		velocityY,
		velocityZ,
		distributions,
		k_Q,
		QQ,
		sizeQ,
		kQ,
		neighborX,
		neighborY,
		neighborZ,
		size_Mat,
		isEvenTimestep);
	getLastCudaError("QVelDevicePlainBB27 execution failed");
}














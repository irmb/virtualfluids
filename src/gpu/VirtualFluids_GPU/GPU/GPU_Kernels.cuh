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
//! \file GPU_Kernels.cuh
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
#ifndef D3Q27_KERNELS_H
#define D3Q27_KERNELS_H

#include "LBM/LB.h"

//////////////////////////////////////////////////////////////////////////
//! \brief \ref Cumulant_K17_LBM_Device_Kernel : Cumulant K17 lattice Boltzmann device kernel function 
extern "C" __global__ void Cumulant_K17_LBM_Device_Kernel(
	real omega,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	int size_Mat,
	real* forces,
	bool isEvenTimestep);

//////////////////////////////////////////////////////////////////////////
//! \brief  \ref LBInit : device function to initialize the distributions
extern "C" __global__ void LBInit(
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
	bool isEvenTimestep);


//////////////////////////////////////////////////////////////////////////
//! \brief \ref LBCalcMacCompSP27 : device function to calculate macroscopic values from distributions
extern "C" __global__ void LBCalcMacCompSP27(
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
	real* distributions,
	bool isEvenTimestep);

//////////////////////////////////////////////////////////////////////////
//! \brief \ref QVelDevPlainBB27 : device function for the velocity boundary condition (plain bounce-back + velocity)
extern "C" __global__ void QVelDevPlainBB27(
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
	bool isEvenTimestep);



#endif
							 
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
//! \file GPU_Interface.h
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
#include "LBM/LB.h"

//////////////////////////////////////////////////////////////////////////
//! \brief Cumulant LBM kernel
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
	bool isEvenTimestep);

//////////////////////////////////////////////////////////////////////////
//! \brief initialize the lattice (distribution functions)
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
	bool isEvenTimestep);
	
//////////////////////////////////////////////////////////////////////////
//! \brief calculates the macroscopic values from distribution functions
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
	bool isEvenTimestep);

//////////////////////////////////////////////////////////////////////////
//! \brief defines the behavior of a velocity boundary condition (plain bounce-back plus velocity)
extern "C" void QVelDevicePlainBB27(
	uint numberOfThreads,
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
	


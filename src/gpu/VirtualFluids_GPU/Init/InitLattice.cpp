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
//! \file InitLattice.h
//! \ingroup Init
//! \author Martin Schoenherr
//=======================================================================================
#include "Init/InitLattice.h"
#include "Parameter/Parameter.h"
#include "GPU/GPU_Interface.h"

////////////////////////////////////////////////////////////////////////////////
void initLattice(SPtr<Parameter> para)
{
	unsigned lev = 0;
	//////////////////////////////////////////////////////////////////////////
	para->getParD()->isEvenTimestep = true;
	//////////////////////////////////////////////////////////////////////////
	LB_Init(
		para->getParD()->numberofthreads, 
		para->getParD()->neighborX,
		para->getParD()->neighborY,
		para->getParD()->neighborZ,
		para->getParD()->typeOfGridNode,
		para->getParD()->rho,
		para->getParD()->velocityX,
		para->getParD()->velocityY,
		para->getParD()->velocityZ,
		para->getParD()->numberOfNodes,
		para->getParD()->distributions.f[0],
		para->getParD()->isEvenTimestep);
	//////////////////////////////////////////////////////////////////////////
	para->getParD()->isEvenTimestep = false;
	//////////////////////////////////////////////////////////////////////////
	LB_Init(
		para->getParD()->numberofthreads,
		para->getParD()->neighborX,
		para->getParD()->neighborY,
		para->getParD()->neighborZ,
		para->getParD()->typeOfGridNode,
		para->getParD()->rho,
		para->getParD()->velocityX,
		para->getParD()->velocityY,
		para->getParD()->velocityZ,
		para->getParD()->numberOfNodes,
		para->getParD()->distributions.f[0],
		para->getParD()->isEvenTimestep);

	//////////////////////////////////////////////////////////////////////////
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
}

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
//! \file CudaKernelManager.cpp
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "CudaKernelManager.h"
#include "GPU_Interface.h"
#include <Parameter/Parameter.h>


void CudaKernelManager::runLBMKernel(SPtr<Parameter> para)
{
    if (para->getIsADcalculationOn()) {
		CumulantK17LBMDeviceKernelAD(
			para->getParD()->numberofthreads,
			para->getParD()->omega,
			para->getParD()->typeOfGridNode,
			para->getParD()->neighborX,
			para->getParD()->neighborY,
			para->getParD()->neighborZ,
			para->getParD()->distributions.f[0],
			para->getParD()->distributionsAD.f[0],
			para->getParD()->numberOfNodes,
			para->getParD()->forcing,
			para->getParD()->isEvenTimestep);
    } else {
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
	}
}

void CudaKernelManager::runVelocityBCKernel(SPtr<Parameter> para)
{
	if (para->getParD()->numberOfVeloBCnodes > 0)
	{
		QVelDevicePlainBB27(
			para->getParD()->numberofthreads,
			para->getParD()->veloBC.Vx,
			para->getParD()->veloBC.Vy,
			para->getParD()->veloBC.Vz,
			para->getParD()->distributions.f[0],
			para->getParD()->veloBC.k,
			para->getParD()->veloBC.q27[0],
			para->getParD()->numberOfVeloBCnodes,
			para->getParD()->veloBC.kArray,
			para->getParD()->neighborX,
			para->getParD()->neighborY,
			para->getParD()->neighborZ,
			para->getParD()->numberOfNodes,
			para->getParD()->isEvenTimestep);
	}
}

void CudaKernelManager::runGeoBCKernel(SPtr<Parameter> para)
{
    if (para->getParD()->numberOfGeoBCnodes > 0)
    {
        // ...
    }
}


void runSlipBCKernel(SPtr<Parameter> para){
    if (para->getParD()->numberOfSlipBCnodes > 0)
    {
        // ...
    }
}

void runNoSlipBCKernel(SPtr<Parameter> para){
    if (para->getParD()->numberOfNoSlipBCnodes > 0)
    {
        // ...
    }
}

void runPressureBCKernel(SPtr<Parameter> para){
    if (para->getParD()->numberOfPressureBCnodes > 0)
    {
        // ...
    }
}

void CudaKernelManager::calculateMacroscopicValues(SPtr<Parameter> para)
{
    if (para->getIsADcalculationOn()) {
		CalcMacADCompSP27(
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
			para->getParD()->distributionsAD.f[0],
            para->getParD()->forcing,
			para->getParD()->isEvenTimestep);
    } else {
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
}











SPtr<CudaKernelManager> CudaKernelManager::make(SPtr<Parameter> parameter)
{
    return SPtr<CudaKernelManager>(new CudaKernelManager(parameter));
}

CudaKernelManager::CudaKernelManager(SPtr<Parameter> parameter)
{
    this->parameter = parameter;
}

CudaKernelManager::CudaKernelManager(const CudaKernelManager&)
{

}

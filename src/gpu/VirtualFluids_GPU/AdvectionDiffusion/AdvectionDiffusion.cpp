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
//! \file AdvectionDiffusion.h
//! \ingroup AdvectionDiffusion
//! \author Martin Schoenherr
//=======================================================================================
#include "AdvectionDiffusion/AdvectionDiffusion.h"
#include "GPU/CudaMemoryManager.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

////////////////////////////////////////////////////////////////////////////////
void AdvectionDiffusion::initAD(int level)
{
    //////////////////////////////////////////////////////////////////////////
    // calculation of omega for diffusivity
    para->getParD(level)->omegaDiffusivity = (real)2.0 / ((real)6.0 * para->getParD(level)->diffusivity + (real)1.0);
    //////////////////////////////////////////////////////////////////////////
    para->getParD(level)->isEvenTimestep = true;
    //////////////////////////////////////////////////////////////////////////
    InitADDev27(
        para->getParD(level)->numberofthreads, 
        para->getParD(level)->neighborX, 
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ, 
        para->getParD(level)->typeOfGridNode, 
        para->getParD(level)->concentration,
        para->getParD(level)->velocityX, 
        para->getParD(level)->velocityY, 
        para->getParD(level)->velocityZ,
        para->getParD(level)->numberOfNodes, 
        para->getParD(level)->distributionsAD27.f[0],
        para->getParD(level)->isEvenTimestep);
    //////////////////////////////////////////////////////////////////////////
    para->getParD(level)->isEvenTimestep = false;
    //////////////////////////////////////////////////////////////////////////
    InitADDev27(
        para->getParD(level)->numberofthreads, 
        para->getParD(level)->neighborX, 
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ, 
        para->getParD(level)->typeOfGridNode, 
        para->getParD(level)->concentration,
        para->getParD(level)->velocityX, 
        para->getParD(level)->velocityY, 
        para->getParD(level)->velocityZ,
        para->getParD(level)->numberOfNodes, 
        para->getParD(level)->distributionsAD27.f[0],
        para->getParD(level)->isEvenTimestep);
    //////////////////////////////////////////////////////////////////////////
    CalcConcentration27(
        para->getParD(level)->numberofthreads,
        para->getParD(level)->concentration,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->distributionsAD27.f[0],
        para->getParD(level)->isEvenTimestep);
}

////////////////////////////////////////////////////////////////////////////////
void AdvectionDiffusion::setInitialNodeValuesAD(int level, SPtr<CudaMemoryManager> cudaMemoryManager)
{
    for (uint j = 1; j <= para->getParH(level)->numberOfNodes; j++) {
        const real coordX = para->getParH(level)->coordinateX[j];
        const real coordY = para->getParH(level)->coordinateY[j];
        const real coordZ = para->getParH(level)->coordinateZ[j];

        real concentration;

        // call functor object with initial condition
        if (para->getInitialConditionAD()) {
            para->getInitialConditionAD()(coordX, coordY, coordZ, concentration);
        } else {
            concentration = real(0.0);
        }

        para->getParH(level)->concentration[j] = concentration;
    }

    cudaMemoryManager->cudaCopyConcentrationHostToDevice();
}

////////////////////////////////////////////////////////////////////////////////
void AdvectionDiffusion::runADcollisionKernel(int level)
{
    FactorizedCentralMomentsAdvectionDiffusionDeviceKernel(
        para->getParD(level)->numberofthreads,
        para->getParD(level)->omegaDiffusivity,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->distributionsAD27.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->forcing,
        para->getParD(level)->isEvenTimestep);
}

void AdvectionDiffusion::runADslipBCKernel(int level){
    if (para->getParD(level)->numberOfSlipBCnodes > 1) {
        ADSlipVelDevComp(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->slipBC.normalX,
            para->getParD(level)->slipBC.normalY,
            para->getParD(level)->slipBC.normalZ,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->distributionsAD27.f[0],
            para->getParD(level)->slipBC.k,
            para->getParD(level)->slipBC.q27[0],
            para->getParD(level)->numberOfSlipBCnodes,
            para->getParD(level)->omegaDiffusivity,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);
    }
}

////////////////////////////////////////////////////////////////////////////////
void AdvectionDiffusion::printAD(int level, SPtr<CudaMemoryManager> cudaMemoryManager)
{
    CalcConcentration27(
        para->getParD(level)->numberofthreads,
        para->getParD(level)->concentration,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->distributionsAD27.f[0],
        para->getParD(level)->isEvenTimestep);

    cudaMemoryManager->cudaCopyConcentrationDeviceToHost();
}

SPtr<AdvectionDiffusion> AdvectionDiffusion::make(SPtr<Parameter> parameter){
    return SPtr<AdvectionDiffusion>(new AdvectionDiffusion(parameter));
}

AdvectionDiffusion::AdvectionDiffusion(SPtr<Parameter> parameter): para(parameter){}

AdvectionDiffusion::AdvectionDiffusion(const AdvectionDiffusion&)
{

}



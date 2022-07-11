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
//! \file ADKernelManager.h
//! \ingroup KernelManager
//! \author Martin Schoenherr
//=======================================================================================
#include "KernelManager/ADKernelManager.h"
#include "GPU/CudaMemoryManager.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

ADKernelManager::ADKernelManager(SPtr<Parameter> parameter): para(parameter){}

void ADKernelManager::initAD(const int level) const
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
void ADKernelManager::setInitialNodeValuesAD(const int level, SPtr<CudaMemoryManager> cudaMemoryManager) const
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

    cudaMemoryManager->cudaCopyConcentrationHostToDevice(level);
}

////////////////////////////////////////////////////////////////////////////////
void ADKernelManager::runADcollisionKernel(const int level)const
{
    if (para->getDiffMod() == 7)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // incompressible
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // KernelADincomp7(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->diffusivity,
        //     para->getParD(level)->typeOfGridNode,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->distributionsAD7.f[0],
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // compressible
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // KernelThS7(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->diffusivity,
        //     para->getParD(level)->typeOfGridNode,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->distributionsAD7.f[0],
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);
    }
    else if (para->getDiffMod() == 27)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // incompressible
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // KernelADincomp27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->diffusivity,
        //     para->getParD(level)->typeOfGridNode,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->distributionsAD27.f[0],
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // compressible
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
}

void ADKernelManager::runADslipBCKernel(const int level) const{
    if (para->getParD(level)->slipBC.numberOfBCnodes > 1) {
        ADSlipVelDevComp(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->slipBC.normalX,
            para->getParD(level)->slipBC.normalY,
            para->getParD(level)->slipBC.normalZ,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->distributionsAD27.f[0],
            para->getParD(level)->slipBC.k,
            para->getParD(level)->slipBC.q27[0],
            para->getParD(level)->slipBC.numberOfBCnodes,
            para->getParD(level)->omegaDiffusivity,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);
    }
}

void ADKernelManager::runADpressureBCKernel(const int level) const{
    if (para->getParD(level)->TempPress.kTemp > 0){
        if (para->getDiffMod() == 7) {
            // QADPressIncompDev7( 
            //     para->getParD(level)->numberofthreads,
            //     para->getParD(level)->distributions.f[0],
            //     para->getParD(level)->distributionsAD7.f[0],
            //     para->getParD(level)->TempPress.temp,
            //     para->getParD(level)->TempPress.velo,
            //     para->getParD(level)->diffusivity,
            //     para->getParD(level)->TempPress.k,
            //     para->getParD(level)->pressureBC.q27[0],
            //     para->getParD(level)->TempPress.kTemp,
            //     para->getParD(level)->omega,
            //     para->getParD(level)->neighborX,
            //     para->getParD(level)->neighborY,
            //     para->getParD(level)->neighborZ,
            //     para->getParD(level)->numberOfNodes,
            //     para->getParD(level)->isEvenTimestep);

             //////////////////////////////////////////////////////////////////////////
             // C O M P R E S S I B L E
             //////////////////////////////////////////////////////////////////////////
            QADPressDev7(
                para->getParD(level)->numberofthreads,
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->distributionsAD7.f[0],
                para->getParD(level)->TempPress.temp,
                para->getParD(level)->TempPress.velo,
                para->getParD(level)->diffusivity,
                para->getParD(level)->TempPress.k,
                para->getParD(level)->pressureBC.q27[0],
                para->getParD(level)->TempPress.kTemp,
                para->getParD(level)->omega,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep);

        } else if (para->getDiffMod() == 27) {
            // QADPressIncompDev27(
            //     para->getParD(level)->numberofthreads,
            //     para->getParD(level)->distributions.f[0],
            //     para->getParD(level)->distributionsAD27.f[0],
            //     para->getParD(level)->TempPress.temp,
            //     para->getParD(level)->TempPress.velo,
            //     para->getParD(level)->diffusivity,
            //     para->getParD(level)->TempPress.k,
            //     para->getParD(level)->pressureBC.q27[0],
            //     para->getParD(level)->TempPress.kTemp,
            //     para->getParD(level)->omega,
            //     para->getParD(level)->neighborX,
            //     para->getParD(level)->neighborY,
            //     para->getParD(level)->neighborZ,
            //     para->getParD(level)->numberOfNodes,
            //     para->getParD(level)->isEvenTimestep);

            //////////////////////////////////////////////////////////////////////////
            // C O M P R E S S I B L E
            //////////////////////////////////////////////////////////////////////////
            QADPressDev27(
                para->getParD(level)->numberofthreads,
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->distributionsAD27.f[0],
                para->getParD(level)->TempPress.temp,
                para->getParD(level)->TempPress.velo,
                para->getParD(level)->diffusivity,
                para->getParD(level)->TempPress.k,
                para->getParD(level)->pressureBC.q27[0],
                para->getParD(level)->TempPress.kTemp,
                para->getParD(level)->omega,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep);
        }
    }
}

void ADKernelManager::runADgeometryBCKernel(const int level) const
{
    if (para->getParD(level)->geometryBC.numberOfBCnodes > 0) {
        if (para->getDiffMod() == 7) {
            // QNoSlipADincompDev7(
            //     para->getParD(level)->numberofthreads,
            //     para->getParD(level)->distributions.f[0],
            //     para->getParD(level)->distributionsAD7.f[0],
            //     para->getParD(level)->Temp.temp,
            //     para->getParD(level)->diffusivity,
            //     para->getParD(level)->Temp.k,
            //     para->getParD(level)->geometryBC.q27[0],
            //     para->getParD(level)->Temp.kTemp,
            //     para->getParD(level)->omega,
            //     para->getParD(level)->neighborX,
            //     para->getParD(level)->neighborY,
            //     para->getParD(level)->neighborZ,
            //     para->getParD(level)->numberOfNodes,
            //     para->getParD(level)->isEvenTimestep);

            //////////////////////////////////////////////////////////////////////////
            // C O M P R E S S I B L E
            //////////////////////////////////////////////////////////////////////////

            QADDev7(
                para->getParD(level)->numberofthreads,
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->distributionsAD7.f[0],
                para->getParD(level)->Temp.temp,
                para->getParD(level)->diffusivity,
                para->getParD(level)->Temp.k,
                para->getParD(level)->geometryBC.q27[0],
                para->getParD(level)->Temp.kTemp,
                para->getParD(level)->omega,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep);

        } else if (para->getDiffMod() == 27) {
            // QNoSlipADincompDev27(
            //     para->getParD(level)->numberofthreads,
            //     para->getParD(level)->distributions.f[0],
            //     para->getParD(level)->distributionsAD27.f[0],
            //     para->getParD(level)->Temp.temp,
            //     para->getParD(level)->diffusivity,
            //     para->getParD(level)->Temp.k,
            //     para->getParD(level)->geometryBC.q27[0],
            //     para->getParD(level)->Temp.kTemp,
            //     para->getParD(level)->omega,
            //     para->getParD(level)->neighborX,
            //     para->getParD(level)->neighborY,
            //     para->getParD(level)->neighborZ,
            //     para->getParD(level)->numberOfNodes,
            //     para->getParD(level)->isEvenTimestep);

            //////////////////////////////////////////////////////////////////////////
            // C O M P R E S S I B L E
            //////////////////////////////////////////////////////////////////////////

            QADBBDev27(
                para->getParD(level)->numberofthreads,
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->distributionsAD27.f[0],
                para->getParD(level)->Temp.temp,
                para->getParD(level)->diffusivity,
                para->getParD(level)->Temp.k,
                para->getParD(level)->geometryBC.q27[0],
                para->getParD(level)->Temp.kTemp,
                para->getParD(level)->omega,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep);
        }
    }
}

void ADKernelManager::runADveloBCKernel(const int level) const{
    if (para->getParD(level)->TempVel.kTemp > 0){
        if (para->getDiffMod() == 7)
        {
            // QADVeloIncompDev7(
            //     para->getParD(level)->numberofthreads,
            //     para->getParD(level)->distributions.f[0],
            //     para->getParD(level)->distributionsAD7.f[0],
            //     para->getParD(level)->TempVel.tempPulse,
            //     para->getParD(level)->TempVel.velo,
            //     para->getParD(level)->diffusivity,
            //     para->getParD(level)->TempVel.k,
            //     para->getParD(level)->velocityBC.q27[0],
            //     para->getParD(level)->TempVel.kTemp,
            //     para->getParD(level)->omega,
            //     para->getParD(level)->neighborX,
            //     para->getParD(level)->neighborY,
            //     para->getParD(level)->neighborZ,
            //     para->getParD(level)->numberOfNodes,
            //     para->getParD(level)->isEvenTimestep);

            //////////////////////////////////////////////////////////////////////////
            // C O M P R E S S I B L E
            //////////////////////////////////////////////////////////////////////////

            QADVelDev7(
                para->getParD(level)->numberofthreads,
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->distributionsAD7.f[0],
                para->getParD(level)->TempVel.temp,
                para->getParD(level)->TempVel.velo,
                para->getParD(level)->diffusivity,
                para->getParD(level)->TempVel.k,
                para->getParD(level)->velocityBC.q27[0],
                para->getParD(level)->TempVel.kTemp,
                para->getParD(level)->omega,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep);

        } else if (para->getDiffMod() == 27) {
            // QADVeloIncompDev27(
            //     para->getParD(level)->numberofthreads,
            //     para->getParD(level)->distributions.f[0],
            //     para->getParD(level)->distributionsAD27.f[0],
            //     para->getParD(level)->TempVel.temp,
            //     para->getParD(level)->TempVel.velo,
            //     para->getParD(level)->diffusivity,
            //     para->getParD(level)->TempVel.k,
            //     para->getParD(level)->velocityBC.q27[0],
            //     para->getParD(level)->TempVel.kTemp,
            //     para->getParD(level)->omega,
            //     para->getParD(level)->neighborX,
            //     para->getParD(level)->neighborY,
            //     para->getParD(level)->neighborZ,
            //     para->getParD(level)->numberOfNodes,
            //     para->getParD(level)->isEvenTimestep);

            //////////////////////////////////////////////////////////////////////////
            // C O M P R E S S I B L E
            //////////////////////////////////////////////////////////////////////////
            QADVelDev27(
                para->getParD(level)->numberofthreads,
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->distributionsAD27.f[0],
                para->getParD(level)->TempVel.tempPulse,
                para->getParD(level)->TempVel.velo,
                para->getParD(level)->diffusivity,
                para->getParD(level)->velocityBC.k,
                para->getParD(level)->velocityBC.q27[0],
                para->getParD(level)->velocityBC.numberOfBCnodes,
                para->getParD(level)->omega,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep);

            //////////////////////////////////////////////////////////////////////////
            // W T G _ R U B
            //////////////////////////////////////////////////////////////////////////
            // if (t<1000)//(t>100000 && t<103895)//(t>1600000 && t<1662317)//(t>500000 && t<515580)//(t<1000)//(t<15580)//(t>400000 && t<415580)//
            // {
            //   QADVelDev27(
            //     para->getParD(level)->numberofthreads,
            //     para->getParD(level)->distributions.f[0],
            //     para->getParD(level)->distributionsAD27.f[0],
            //     para->getParD(level)->TempVel.tempPulse,
            //     para->getParD(level)->TempVel.velo,
            //     para->getParD(level)->diffusivity,
            //     para->getParD(level)->velocityBC.k,
            //     para->getParD(level)->velocityBC.q27[0],
            //     para->getParD(level)->velocityBC.numberOfBCnodes,
            //     para->getParD(level)->omega,
            //     para->getParD(level)->neighborX,
            //     para->getParD(level)->neighborY,
            //     para->getParD(level)->neighborZ,
            //     para->getParD(level)->numberOfNodes,
            //     para->getParD(level)->isEvenTimestep);
            // }
            // else
            // {
            //   QADVelDev27(
            //     para->getParD(level)->numberofthreads,
            //     para->getParD(level)->distributions.f[0],
            //     para->getParD(level)->distributionsAD27.f[0],
            //     para->getParD(level)->TempVel.temp,
            //     para->getParD(level)->TempVel.velo,
            //     para->getParD(level)->diffusivity,
            //     para->getParD(level)->velocityBC.k,
            //     para->getParD(level)->velocityBC.q27[0],
            //     para->getParD(level)->velocityBC.numberOfBCnodes,
            //     para->getParD(level)->omega,
            //     para->getParD(level)->neighborX,
            //     para->getParD(level)->neighborY,
            //     para->getParD(level)->neighborZ,
            //     para->getParD(level)->numberOfNodes,
            //     para->getParD(level)->isEvenTimestep);
            // }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void ADKernelManager::printAD(const int level, SPtr<CudaMemoryManager> cudaMemoryManager) const
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

    cudaMemoryManager->cudaCopyConcentrationDeviceToHost(level);
}

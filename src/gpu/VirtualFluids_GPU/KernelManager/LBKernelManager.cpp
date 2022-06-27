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
//! \file LBKernelManager.cpp
//! \ingroup KernelManager
//! \author Martin Schoenherr
//=======================================================================================
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "LBKernelManager.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "Calculation/DragLift.h"
#include "Calculation/Cp.h"
#include "PreCollisionInteractor/VelocitySetter.h"



void LBKernelManager::runLBMKernel(int level)
{
    // if (para->getIsADcalculationOn()) {
    //       CumulantK17LBMDeviceKernelAD(
    //            para->getParD()->numberofthreads,
    //            para->getParD()->omega,
    //            para->getParD()->typeOfGridNode,
    //            para->getParD()->neighborX,
    //            para->getParD()->neighborY,
    //            para->getParD()->neighborZ,
    //            para->getParD()->distributions.f[0],
    //            para->getParD()->distributionsAD.f[0],
    //            para->getParD()->numberOfNodes,
    //            para->getParD()->forcing,
    //            para->getParD()->isEvenTimestep);
    // } else {
    //       CumulantK17LBMDeviceKernel(
    //            para->getParD()->numberofthreads,
    //            para->getParD()->omega,
    //            para->getParD()->typeOfGridNode,
    //            para->getParD()->neighborX,
    //            para->getParD()->neighborY,
    //            para->getParD()->neighborZ,
    //            para->getParD()->distributions.f[0],
    //            para->getParD()->numberOfNodes,
    //            para->getParD()->forcing,
    //            para->getParD()->isEvenTimestep);
    //  }
}

void LBKernelManager::runVelocityBCKernelPre(int level)
{
    if (para->getParD(level)->numberOfVeloBCnodes > 0)
    {
        // TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
        // if ( myid == 0)
        // {
        //    VelSchlaffer27(para->getParD(level)->numberofthreads, t,
        //                   para->getParD(level)->distributions.f[0],       para->getParD(level)->velocityBC.Vz,
        //                   para->getParD(level)->velocityBC.deltaVz, para->getParD(level)->velocityBC.k,
        //                   para->getParD(level)->velocityBC.kN,      para->getParD(level)->numberOfVeloBCnodes,
        //                   para->getParD(level)->omega,           para->getParD(level)->neighborX,
        //                   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
        //                   para->getParD(level)->numberOfNodes,     para->getParD(level)->isEvenTimestep);
        //    getLastCudaError("VelSchlaffer27 execution failed");
        // }
        ////////////////////////////////////////////////////////////////////////////
        // high viscosity incompressible
        // QVelDevIncompHighNu27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->velocityBC.Vx,
        //     para->getParD(level)->velocityBC.Vy,
        //     para->getParD(level)->velocityBC.Vz,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->velocityBC.k,
        //     para->getParD(level)->velocityBC.q27[0],
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

        ////////////////////////////////////////////////////////////////////////////
        // high viscosity compressible
        // QVelDevCompHighNu27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->velocityBC.Vx,
        //     para->getParD(level)->velocityBC.Vy,
        //     para->getParD(level)->velocityBC.Vz,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->velocityBC.k,
        //     para->getParD(level)->velocityBC.q27[0],
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);
    }
}

void LBKernelManager::runVelocityBCKernelPost(int level)
{
     if (para->getParD(level)->numberOfVeloBCnodes > 0)
     {
        //   QVelDevicePlainBB27(
        //     para->getParD()->numberofthreads,
        //     para->getParD()->velocityBC.Vx,
        //     para->getParD()->velocityBC.Vy,
        //     para->getParD()->velocityBC.Vz,
        //     para->getParD()->distributions.f[0],
        //     para->getParD()->velocityBC.k,
        //     para->getParD()->velocityBC.q27[0],
        //     para->getParD()->numberOfVeloBCnodes,
        //     para->getParD()->velocityBC.kArray,
        //     para->getParD()->neighborX,
        //     para->getParD()->neighborY,
        //     para->getParD()->neighborZ,
        //     para->getParD()->numberOfNodes,
        //     para->getParD()->isEvenTimestep);

        // QVelDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->velocityBC.Vx,
        //     para->getParD(level)->velocityBC.Vy,
        //     para->getParD(level)->velocityBC.Vz,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->velocityBC.k,
        //     para->getParD(level)->velocityBC.q27[0],
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        // QVelDevComp27(
        //     para->getParD(level)->numberofthreads, para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->velocityBC.Vx,
        //     para->getParD(level)->velocityBC.Vy,
        //     para->getParD(level)->velocityBC.Vz,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->velocityBC.k,
        //     para->getParD(level)->velocityBC.q27[0],
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        QVelDevCompZeroPress27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->nx,
            para->getParD(level)->ny,
            para->getParD(level)->velocityBC.Vx,
            para->getParD(level)->velocityBC.Vy,
            para->getParD(level)->velocityBC.Vz,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->velocityBC.k,
            para->getParD(level)->velocityBC.q27[0],
            para->getParD(level)->numberOfVeloBCnodes,
            para->getParD(level)->velocityBC.kArray,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);

        //////////////////////////////////////////////////////////////////////////
        // D E P R E C A T E D
        //////////////////////////////////////////////////////////////////////////

        //QVelDevice1h27( para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
        //                para->getParD(level)->velocityBC.Vx,      para->getParD(level)->velocityBC.Vy,   para->getParD(level)->velocityBC.Vz,
        //                para->getParD(level)->distributions.f[0],       para->getParD(level)->velocityBC.k,    para->getParD(level)->velocityBC.q27[0],
        //                para->getParD(level)->numberOfVeloBCnodes,        para->getParD(level)->numberOfVeloBCnodes,     para->getParD(level)->omega,
        //                para->getPhi(),                        para->getAngularVelocity(),
        //                para->getParD(level)->neighborX,    para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
        //                para->getParD(level)->coordinateX,       para->getParD(level)->coordinateY,    para->getParD(level)->coordinateZ,
        //                para->getParD(level)->size_Mat_SP,     para->getParD(level)->isEvenTimestep);
        //getLastCudaError("QVelDev27 execution failed");
     }
}

void LBKernelManager::runGeoBCKernelPre(int level, unsigned int t, CudaMemoryManager* cudaMemoryManager){
    if (para->getParD(level)->geometryBC.numberOfBCnodes > 0){
        if (para->getCalcDragLift())
        {
            //Drag and Lift Part II
            DragLiftPreD27(
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->geometryBC.k,
                para->getParD(level)->geometryBC.q27[0],
                para->getParD(level)->geometryBC.numberOfBCnodes,
                para->getParD(level)->DragPreX,
                para->getParD(level)->DragPreY,
                para->getParD(level)->DragPreZ,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep,
                para->getParD(level)->numberofthreads);
            ////////////////////////////////////////////////////////////////////////////////
            //Calculation of Drag and Lift
            ////////////////////////////////////////////////////////////////////////////////
            calcDragLift(para.get(), cudaMemoryManager, level);
            ////////////////////////////////////////////////////////////////////////////////
        }

        if (para->getCalcCp())
        {
            ////////////////////////////////////////////////////////////////////////////////
            //Calculation of cp
            ////////////////////////////////////////////////////////////////////////////////

            if(t > para->getTStartOut())
            {
                ////////////////////////////////////////////////////////////////////////////////
                CalcCPtop27(
                    para->getParD(level)->distributions.f[0],
                    para->getParD(level)->cpTopIndex,
                    para->getParD(level)->numberOfPointsCpTop,
                    para->getParD(level)->cpPressTop,
                    para->getParD(level)->neighborX,
                    para->getParD(level)->neighborY,
                    para->getParD(level)->neighborZ,
                    para->getParD(level)->numberOfNodes,
                    para->getParD(level)->isEvenTimestep,
                    para->getParD(level)->numberofthreads);
                //////////////////////////////////////////////////////////////////////////////////
                CalcCPbottom27(
                    para->getParD(level)->distributions.f[0],
                    para->getParD(level)->cpBottomIndex,
                    para->getParD(level)->numberOfPointsCpBottom,
                    para->getParD(level)->cpPressBottom,
                    para->getParD(level)->neighborX,
                    para->getParD(level)->neighborY,
                    para->getParD(level)->neighborZ,
                    para->getParD(level)->numberOfNodes,
                    para->getParD(level)->isEvenTimestep,
                    para->getParD(level)->numberofthreads);
                //////////////////////////////////////////////////////////////////////////////////
                CalcCPbottom27(
                    para->getParD(level)->distributions.f[0],
                    para->getParD(level)->cpBottom2Index,
                    para->getParD(level)->numberOfPointsCpBottom2,
                    para->getParD(level)->cpPressBottom2,
                    para->getParD(level)->neighborX,
                    para->getParD(level)->neighborY,
                    para->getParD(level)->neighborZ,
                    para->getParD(level)->numberOfNodes,
                    para->getParD(level)->isEvenTimestep,
                    para->getParD(level)->numberofthreads);
                //////////////////////////////////////////////////////////////////////////////////
                calcCp(para.get(), cudaMemoryManager, level);
            }            
        }

        ////////////////////////////////////////////////////////////////////////////////
        // high viscosity incompressible
        // QDevIncompHighNu27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->geometryBC.k,
        //     para->getParD(level)->geometryBC.q27[0],
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

        //////////////////////////////////////////////////////////////////////////////////
        // high viscosity compressible
        // QDevCompHighNu27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->geometryBC.k,
        //     para->getParD(level)->geometryBC.q27[0],
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

    }
}

void LBKernelManager::runGeoBCKernelPost(int level)
{
    if (para->getParD(level)->geometryBC.numberOfBCnodes > 0)
    {
        if (para->getCalcDragLift())
        {
            //Drag and Lift Part I
            DragLiftPostD27(para->getParD(level)->distributions.f[0],
                            para->getParD(level)->geometryBC.k,
                            para->getParD(level)->geometryBC.q27[0],
                            para->getParD(level)->geometryBC.numberOfBCnodes,
                            para->getParD(level)->DragPostX,
                            para->getParD(level)->DragPostY,
                            para->getParD(level)->DragPostZ,
                            para->getParD(level)->neighborX,
                            para->getParD(level)->neighborY,
                            para->getParD(level)->neighborZ,
                            para->getParD(level)->numberOfNodes,
                            para->getParD(level)->isEvenTimestep,
                            para->getParD(level)->numberofthreads);
            getLastCudaError("DragLift27 execution failed");
        }

        // BBDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->geometryBC.k,
        //     para->getParD(level)->geometryBC.q27[0],
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        // QDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->geometryBC.k,
        //     para->getParD(level)->geometryBC.q27[0],
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        // QVelDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->geometryBC.Vx,
        //     para->getParD(level)->geometryBC.Vy,
        //     para->getParD(level)->geometryBC.Vz,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->geometryBC.k,
        //     para->getParD(level)->geometryBC.q27[0],
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        // QDevComp27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->geometryBC.k,
        //     para->getParD(level)->geometryBC.q27[0],
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        QVelDevComp27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->nx,
            para->getParD(level)->ny,
            para->getParD(level)->geometryBC.Vx,
            para->getParD(level)->geometryBC.Vy,
            para->getParD(level)->geometryBC.Vz,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->geometryBC.k,
            para->getParD(level)->geometryBC.q27[0],
            para->getParD(level)->geometryBC.numberOfBCnodes,
            para->getParD(level)->geometryBC.numberOfBCnodes,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);

    //     QVelDevCompZeroPress27(
    //         para->getParD(0)->numberofthreads, para->getParD(0)->nx,
    //         para->getParD(0)->ny,
    //         para->getParD(0)->geometryBC.Vx,
    //         para->getParD(0)->geometryBC.Vy,
    //         para->getParD(0)->geometryBC.Vz,
    //         para->getParD(0)->distributions.f[0],
    //         para->getParD(0)->geometryBC.k,
    //         para->getParD(0)->geometryBC.q27[0],
    //         para->getParD(0)->geometryBC.numberOfBCnodes,
    //         para->getParD(0)->geometryBC.numberOfBCnodes,
    //         para->getParD(0)->omega,
    //         para->getParD(0)->neighborX,
    //         para->getParD(0)->neighborY,
    //         para->getParD(0)->neighborZ,
    //         para->getParD(0)->size_Mat_SP,
    //         para->getParD(0)->isEvenTimestep);

    //     QDev3rdMomentsComp27(
    //         para->getParD(level)->numberofthreads,
    //         para->getParD(level)->nx,
    //         para->getParD(level)->ny,
    //         para->getParD(level)->distributions.f[0],
    //         para->getParD(level)->geometryBC.k,
    //         para->getParD(level)->geometryBC.q27[0],
    //         para->getParD(level)->geometryBC.numberOfBCnodes,
    //         para->getParD(level)->geometryBC.numberOfBCnodes,
    //         para->getParD(level)->omega,
    //         para->getParD(level)->neighborX,
    //         para->getParD(level)->neighborY,
    //         para->getParD(level)->neighborZ,
    //         para->getParD(level)->size_Mat_SP,
    //         para->getParD(level)->isEvenTimestep);

    //     QSlipDev27(
    //         para->getParD(level)->numberofthreads,
    //         para->getParD(level)->distributions.f[0],
    //         para->getParD(level)->geometryBC.k,
    //         para->getParD(level)->geometryBC.q27[0],
    //         para->getParD(level)->geometryBC.numberOfBCnodes,
    //         para->getParD(level)->omega,
    //         para->getParD(level)->neighborX,
    //         para->getParD(level)->neighborY,
    //         para->getParD(level)->neighborZ,
    //         para->getParD(level)->size_Mat_SP,
    //         para->getParD(level)->isEvenTimestep);

    //////////////////////////////////////////////////////////////////////////
    // D E P R E C A T E D
    //////////////////////////////////////////////////////////////////////////
    // the GridGenerator does currently not provide normals!

    //     QSlipGeomDevComp27(
    //         para->getParD(level)->numberofthreads,
    //         para->getParD(level)->distributions.f[0],
    //         para->getParD(level)->geometryBC.k,
    //         para->getParD(level)->geometryBC.q27[0],
    //         para->getParD(level)->geometryBC.numberOfBCnodes,
    //         para->getParD(level)->omega,
    //         para->getParD(level)->geometryBCnormalX.q27[0],
    //         para->getParD(level)->geometryBCnormalY.q27[0],
    //         para->getParD(level)->geometryBCnormalZ.q27[0],
    //         para->getParD(level)->neighborX,
    //         para->getParD(level)->neighborY,
    //         para->getParD(level)->neighborZ,
    //         para->getParD(level)->size_Mat_SP,
    //         para->getParD(level)->isEvenTimestep);

    //     QSlipNormDevComp27(
    //         para->getParD(level)->numberofthreads,
    //         para->getParD(level)->distributions.f[0],
    //         para->getParD(level)->geometryBC.k,
    //         para->getParD(level)->geometryBC.q27[0],
    //         para->getParD(level)->geometryBC.numberOfBCnodes,
    //         para->getParD(level)->omega,
    //         para->getParD(level)->geometryBCnormalX.q27[0],
    //         para->getParD(level)->geometryBCnormalY.q27[0],
    //         para->getParD(level)->geometryBCnormalZ.q27[0],
    //         para->getParD(level)->neighborX,
    //         para->getParD(level)->neighborY,
    //         para->getParD(level)->neighborZ,
    //         para->getParD(level)->size_Mat_SP,
    //         para->getParD(level)->isEvenTimestep);
    }
}

void LBKernelManager::runOutflowBCKernelPre(int level){
    if (para->getParD(level)->numberOfOutflowBCnodes > 0)
    {
        QPressNoRhoDev27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->outflowBC.RhoBC,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->outflowBC.k,
            para->getParD(level)->outflowBC.kN,
            para->getParD(level)->outflowBC.numberOfBCnodes,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);

        // TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
        // if (  myid == numprocs - 1)
        // PressSchlaffer27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->outflowBC.RhoBC,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->outflowBC.Vx,
        //     para->getParD(level)->outflowBC.Vy,
        //     para->getParD(level)->outflowBC.Vz,
        //     para->getParD(level)->outflowBC.deltaVz,
        //     para->getParD(level)->outflowBC.k,
        //     para->getParD(level)->outflowBC.kN,
        //     para->getParD(level)->numberOfOutflowBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);
    }
}

void LBKernelManager::runPressureBCKernelPre(int level){
    if (para->getParD(level)->pressureBC.numberOfBCnodes > 0)
    {
        QPressNoRhoDev27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->pressureBC.RhoBC,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->pressureBC.k,
            para->getParD(level)->pressureBC.kN,
            para->getParD(level)->pressureBC.numberOfBCnodes,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);

        // QPressDevEQZ27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->pressureBC.RhoBC,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->pressureBC.k,
        //     para->getParD(level)->pressureBC.kN,
        //     para->getParD(level)->kDistTestRE.f[0],
        //     para->getParD(level)->pressureBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(0)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

        // QInflowScaleByPressDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->pressureBC.RhoBC,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->pressureBC.k,
        //     para->getParD(level)->pressureBC.kN,
        //     para->getParD(level)->pressureBC.numberOfBCnodes,
        //     para->getParD(0)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(0)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

        ////////////////////////////////////////////////////////////////////////////////
        //press NEQ incompressible
        // QPressDevIncompNEQ27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->pressureBC.RhoBC,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->pressureBC.k,
        //     para->getParD(level)->pressureBC.kN,
        //     para->getParD(level)->pressureBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

        //////////////////////////////////////////////////////////////////////////////////
        //press NEQ compressible
        // QPressDevNEQ27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->pressureBC.RhoBC,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->pressureBC.k,
        //     para->getParD(level)->pressureBC.kN,
        //     para->getParD(level)->pressureBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);
    }
}

void LBKernelManager::runPressureBCKernelPost(int level){
    if (para->getParD(level)->numberOfPressureBCnodes > 0)
    {
        QPressDev27_IntBB(
            para->getParD(level)->numberofthreads, 
            para->getParD(level)->pressureBC.RhoBC,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->pressureBC.k,
            para->getParD(level)->pressureBC.q27[0],
            para->getParD(level)->pressureBC.numberOfBCnodes,
            para->getParD(level)->pressureBC.numberOfBCnodes,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);
    }
}

void LBKernelManager::runStressWallModelKernel(int level){
    if (para->getParD(level)->numberOfStressBCnodes > 0)
    {
        // QStressDevComp27(para->getParD(level)->numberofthreads, para->getParD(level)->distributions.f[0],
        //                 para->getParD(level)->stressBC.k,       para->getParD(level)->stressBC.kN,
        //                 para->getParD(level)->stressBC.q27[0],  para->getParD(level)->numberOfStressBCnodes,
        //                 para->getParD(level)->omega,           para->getParD(level)->turbViscosity,
        //                 para->getParD(level)->velocityX,           para->getParD(level)->velocityY,             para->getParD(level)->velocityY,
        //                 para->getParD(level)->stressBC.normalX, para->getParD(level)->stressBC.normalY,   para->getParD(level)->stressBC.normalZ,
        //                 para->getParD(level)->stressBC.Vx,      para->getParD(level)->stressBC.Vy,        para->getParD(level)->stressBC.Vz,
        //                 para->getParD(level)->stressBC.Vx1,     para->getParD(level)->stressBC.Vy1,       para->getParD(level)->stressBC.Vz1,
        //                 para->getParD(level)->wallModel.samplingOffset, para->getParD(level)->wallModel.z0,
        //                 para->getHasWallModelMonitor(),        para->getParD(level)->wallModel.u_star,
        //                 para->getParD(level)->wallModel.Fx,    para->getParD(level)->wallModel.Fy,      para->getParD(level)->wallModel.Fz,
        //                 para->getParD(level)->neighborX,    para->getParD(level)->neighborY,      para->getParD(level)->neighborZ,
        //                 para->getParD(level)->size_Mat_SP,     para->getParD(level)->isEvenTimestep);

        BBStressDev27( para->getParD(level)->numberofthreads, para->getParD(level)->distributions.f[0],
                        para->getParD(level)->stressBC.k,       para->getParD(level)->stressBC.kN,
                        para->getParD(level)->stressBC.q27[0],  para->getParD(level)->numberOfStressBCnodes,
                        para->getParD(level)->velocityX,           para->getParD(level)->velocityY,             para->getParD(level)->velocityY,
                        para->getParD(level)->stressBC.normalX, para->getParD(level)->stressBC.normalY,   para->getParD(level)->stressBC.normalZ,
                        para->getParD(level)->stressBC.Vx,      para->getParD(level)->stressBC.Vy,        para->getParD(level)->stressBC.Vz,
                        para->getParD(level)->stressBC.Vx1,     para->getParD(level)->stressBC.Vy1,       para->getParD(level)->stressBC.Vz1,
                        para->getParD(level)->wallModel.samplingOffset, para->getParD(level)->wallModel.z0,
                        para->getHasWallModelMonitor(),        para->getParD(level)->wallModel.u_star,
                        para->getParD(level)->wallModel.Fx,    para->getParD(level)->wallModel.Fy,      para->getParD(level)->wallModel.Fz,
                        para->getParD(level)->neighborX,    para->getParD(level)->neighborY,      para->getParD(level)->neighborZ,
                        para->getParD(level)->numberOfNodes,     para->getParD(level)->isEvenTimestep);
    }
}

void LBKernelManager::runPrecursorBCKernelPost(int level, uint t, CudaMemoryManager* cudaMemoryManager)
{
    if(para->getParH(level)->numberOfPrecursorBCnodes > 0)
    {
        auto precursorStruct = para->getParD(level)->precursorBC;

        int tLastRead = precursorStruct.nReads*precursorStruct.nTRead;

        if(t>tLastRead+precursorStruct.nTRead)
        {
            real nextTime = (t+1)*pow(2,-level)*para->getTimeRatio();
            // TODO add sync barrier here
            //cycle pointers
            real* tmp = precursorStruct.vxLast;
            precursorStruct.vxLast = precursorStruct.vxCurrent;
            precursorStruct.vxCurrent = precursorStruct.vxNext;
            precursorStruct.vxNext = tmp;


            tmp = precursorStruct.vyLast;
            precursorStruct.vyLast = precursorStruct.vyCurrent;
            precursorStruct.vyCurrent = precursorStruct.vyNext;
            precursorStruct.vyNext = tmp;
            
            tmp = precursorStruct.vzLast;
            precursorStruct.vzLast = precursorStruct.vzCurrent;
            precursorStruct.vzCurrent = precursorStruct.vzNext;
            precursorStruct.vzNext = tmp;

            for(auto reader : para->getParH(level)->velocityReader)
            {   
                reader->getNextVelocities(precursorStruct.vxNext, precursorStruct.vyNext, precursorStruct.vzNext, nextTime);
            }

            cudaMemoryManager->cudaCopyPrecursorVelocities(level);
            precursorStruct.nReads++;
            tLastRead += precursorStruct.nTRead;  
        }

        real tRatio = real(t-tLastRead)/precursorStruct.nTRead;
        QPrecursorDevCompZeroPress( para->getParD(level)->numberofthreads, tRatio, para->getParD(level)->distributions.f[0], para->getParD(level)->precursorBC.q27[0],
                                    para->getParD(level)->precursorBC.k, para->getParD(level)->precursorBC.kQ, para->getParD(level)->precursorBC.kArray, 
                                    para->getParD(level)->precursorBC.nVelocityPoints, para->getParH(level)->omega, para->getVelocityRatio(),
                                    para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                    para->getParD(level)->precursorBC.planeNeighborNT, para->getParD(level)->precursorBC.planeNeighborNB, para->getParD(level)->precursorBC.planeNeighborST, para->getParD(level)->precursorBC.planeNeighborSB, 
                                    para->getParD(level)->precursorBC.weightsNT, para->getParD(level)->precursorBC.weightsNB, para->getParD(level)->precursorBC.weightsST, para->getParD(level)->precursorBC.weightsSB, 
                                    para->getParD(level)->precursorBC.vxLast, para->getParD(level)->precursorBC.vyLast, para->getParD(level)->precursorBC.vzLast, 
                                    para->getParD(level)->precursorBC.vxCurrent, para->getParD(level)->precursorBC.vyCurrent, para->getParD(level)->precursorBC.vzCurrent, 
                                    para->getParD(level)->numberOfNodes,     para->getParD(level)->isEvenTimestep);
        getLastCudaError("QPrecursorDevCompZeroPress execution failed");
    }
}

void LBKernelManager::runSlipBCKernel(int level){
    if (para->getParD(level)->numberOfSlipBCnodes > 0)
    {
        // QSlipDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->slipBC.k,
        //     para->getParD(level)->slipBC.q27[0],
        //     para->getParD(level)->numberOfSlipBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        QSlipDevComp27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->slipBC.k,
            para->getParD(level)->slipBC.q27[0],
            para->getParD(level)->numberOfSlipBCnodes,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->turbViscosity,
            para->getUseTurbulentViscosity(),
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);
    }
}

void LBKernelManager::runNoSlipBCKernel(int level){
    if (para->getParD(level)->numberOfNoSlipBCnodes > 0)
    {
        // QDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->noSlipBC.k,
        //     para->getParD(level)->noSlipBC.q27[0],
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        // BBDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->noSlipBC.k,
        //     para->getParD(level)->noSlipBC.q27[0],
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        // QDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->noSlipBC.k,
        //     para->getParD(level)->noSlipBC.q27[0],
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->isEvenTimestep);

        QDevComp27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->nx,
            para->getParD(level)->ny,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->noSlipBC.k,
            para->getParD(level)->noSlipBC.q27[0],
            para->getParD(level)->numberOfNoSlipBCnodes,
            para->getParD(level)->numberOfNoSlipBCnodes,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);
    }
}

// void LBKernelManager::runPressureBCKernelPre(int level){
//     if (para->getParD()->numberOfPressureBCnodes > 0)
//     {
//         // ...
//     }
// }

// void LBKernelManager::calculateMacroscopicValues(int level)
// {
//     if (para->getIsADcalculationOn()) {
//           CalcMacADCompSP27(
//                para->getParD()->velocityX,
//                para->getParD()->velocityY,
//                para->getParD()->velocityZ,
//                para->getParD()->rho,
//                para->getParD()->pressure,
//                para->getParD()->typeOfGridNode,
//                para->getParD()->neighborX,
//                para->getParD()->neighborY,
//                para->getParD()->neighborZ,
//                para->getParD()->numberOfNodes,
//                para->getParD()->numberofthreads,
//                para->getParD()->distributions.f[0],
//                para->getParD()->distributionsAD.f[0],
//             para->getParD()->forcing,
//                para->getParD()->isEvenTimestep);
//     } else {
//           CalcMacCompSP27(
//                para->getParD()->velocityX,
//                para->getParD()->velocityY,
//                para->getParD()->velocityZ,
//                para->getParD()->rho,
//                para->getParD()->pressure,
//                para->getParD()->typeOfGridNode,
//                para->getParD()->neighborX,
//                para->getParD()->neighborY,
//                para->getParD()->neighborZ,
//                para->getParD()->numberOfNodes,
//                para->getParD()->numberofthreads,
//                para->getParD()->distributions.f[0],
//                para->getParD()->isEvenTimestep);
//      }
// }









SPtr<LBKernelManager> LBKernelManager::make(SPtr<Parameter> parameter)
{
    return SPtr<LBKernelManager>(new LBKernelManager(parameter));
}

LBKernelManager::LBKernelManager(SPtr<Parameter> parameter)
{
    this->para = parameter;
}

LBKernelManager::LBKernelManager(const LBKernelManager&)
{

}

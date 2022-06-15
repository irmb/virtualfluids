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


void CudaKernelManager::runLBMKernel(int level)
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

void CudaKernelManager::runVelocityBCKernel(int level)
{
     if (para->getParD(level)->numberOfVeloBCnodes > 0)
     {
        //   QVelDevicePlainBB27(
        //     para->getParD()->numberofthreads,
        //     para->getParD()->veloBC.Vx,
        //     para->getParD()->veloBC.Vy,
        //     para->getParD()->veloBC.Vz,
        //     para->getParD()->distributions.f[0],
        //     para->getParD()->veloBC.k,
        //     para->getParD()->veloBC.q27[0],
        //     para->getParD()->numberOfVeloBCnodes,
        //     para->getParD()->veloBC.kArray,
        //     para->getParD()->neighborX,
        //     para->getParD()->neighborY,
        //     para->getParD()->neighborZ,
        //     para->getParD()->numberOfNodes,
        //     para->getParD()->isEvenTimestep);

        // QVelDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->Qinflow.Vx,
        //     para->getParD(level)->Qinflow.Vy,
        //     para->getParD(level)->Qinflow.Vz,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->Qinflow.k,
        //     para->getParD(level)->Qinflow.q27[0],
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        // QVelDevComp27(
        //     para->getParD(level)->numberofthreads, para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->Qinflow.Vx,
        //     para->getParD(level)->Qinflow.Vy,
        //     para->getParD(level)->Qinflow.Vz,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->Qinflow.k,
        //     para->getParD(level)->Qinflow.q27[0],
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->numberOfVeloBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        QVelDevCompZeroPress27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->nx,
            para->getParD(level)->ny,
            para->getParD(level)->Qinflow.Vx,
            para->getParD(level)->Qinflow.Vy,
            para->getParD(level)->Qinflow.Vz,
            para->getParD(level)->d0SP.f[0],
            para->getParD(level)->Qinflow.k,
            para->getParD(level)->Qinflow.q27[0],
            para->getParD(level)->numberOfVeloBCnodes,
            para->getParD(level)->Qinflow.kArray,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd);

        //////////////////////////////////////////////////////////////////////////
        // D E P R E C A T E D
        //////////////////////////////////////////////////////////////////////////

        //QVelDevice1h27( para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
        //                para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
        //                para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0],
        //                para->getParD(level)->numberOfVeloBCnodes,        para->getParD(level)->numberOfVeloBCnodes,     para->getParD(level)->omega,
        //                para->getPhi(),                        para->getAngularVelocity(),
        //                para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                para->getParD(level)->coordX_SP,       para->getParD(level)->coordY_SP,    para->getParD(level)->coordZ_SP,
        //                para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QVelDev27 execution failed");
     }
}

void CudaKernelManager::runGeoBCKernelPost(int level)
{
    if (para->getParD(level)->QGeom.numberOfBCnodes > 0)
    {
        if (para->getCalcDragLift())
        {
            //Drag and Lift Part I
            DragLiftPostD27(para->getParD(level)->d0SP.f[0],
                            para->getParD(level)->QGeom.k,
                            para->getParD(level)->QGeom.q27[0],
                            para->getParD(level)->QGeom.numberOfBCnodes,
                            para->getParD(level)->DragPostX,
                            para->getParD(level)->DragPostY,
                            para->getParD(level)->DragPostZ,
                            para->getParD(level)->neighborX_SP,
                            para->getParD(level)->neighborY_SP,
                            para->getParD(level)->neighborZ_SP,
                            para->getParD(level)->size_Mat_SP,
                            para->getParD(level)->evenOrOdd,
                            para->getParD(level)->numberofthreads);
            getLastCudaError("DragLift27 execution failed");
        }

        // BBDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->QGeom.k,
        //     para->getParD(level)->QGeom.q27[0],
        //     para->getParD(level)->QGeom.numberOfBCnodes,
        //     para->getParD(level)->QGeom.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        // QDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->QGeom.k,
        //     para->getParD(level)->QGeom.q27[0],
        //     para->getParD(level)->QGeom.numberOfBCnodes,
        //     para->getParD(level)->QGeom.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        // QVelDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->QGeom.Vx,
        //     para->getParD(level)->QGeom.Vy,
        //     para->getParD(level)->QGeom.Vz,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->QGeom.k,
        //     para->getParD(level)->QGeom.q27[0],
        //     para->getParD(level)->QGeom.numberOfBCnodes,
        //     para->getParD(level)->QGeom.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        // QDevComp27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->QGeom.k,
        //     para->getParD(level)->QGeom.q27[0],
        //     para->getParD(level)->QGeom.numberOfBCnodes,
        //     para->getParD(level)->QGeom.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        QVelDevComp27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->nx,
            para->getParD(level)->ny,
            para->getParD(level)->QGeom.Vx,
            para->getParD(level)->QGeom.Vy,
            para->getParD(level)->QGeom.Vz,
            para->getParD(level)->d0SP.f[0],
            para->getParD(level)->QGeom.k,
            para->getParD(level)->QGeom.q27[0],
            para->getParD(level)->QGeom.numberOfBCnodes,
            para->getParD(level)->QGeom.numberOfBCnodes,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd);

    //     QVelDevCompZeroPress27(
    //         para->getParD(0)->numberofthreads, para->getParD(0)->nx,
    //         para->getParD(0)->ny,
    //         para->getParD(0)->QGeom.Vx,
    //         para->getParD(0)->QGeom.Vy,
    //         para->getParD(0)->QGeom.Vz,
    //         para->getParD(0)->d0SP.f[0],
    //         para->getParD(0)->QGeom.k,
    //         para->getParD(0)->QGeom.q27[0],
    //         para->getParD(0)->QGeom.numberOfBCnodes,
    //         para->getParD(0)->QGeom.numberOfBCnodes,
    //         para->getParD(0)->omega,
    //         para->getParD(0)->neighborX_SP,
    //         para->getParD(0)->neighborY_SP,
    //         para->getParD(0)->neighborZ_SP,
    //         para->getParD(0)->size_Mat_SP,
    //         para->getParD(0)->evenOrOdd);

    //     QDev3rdMomentsComp27(
    //         para->getParD(level)->numberofthreads,
    //         para->getParD(level)->nx,
    //         para->getParD(level)->ny,
    //         para->getParD(level)->d0SP.f[0],
    //         para->getParD(level)->QGeom.k,
    //         para->getParD(level)->QGeom.q27[0],
    //         para->getParD(level)->QGeom.numberOfBCnodes,
    //         para->getParD(level)->QGeom.numberOfBCnodes,
    //         para->getParD(level)->omega,
    //         para->getParD(level)->neighborX_SP,
    //         para->getParD(level)->neighborY_SP,
    //         para->getParD(level)->neighborZ_SP,
    //         para->getParD(level)->size_Mat_SP,
    //         para->getParD(level)->evenOrOdd);

    //     QSlipDev27(
    //         para->getParD(level)->numberofthreads,
    //         para->getParD(level)->d0SP.f[0],
    //         para->getParD(level)->QGeom.k,
    //         para->getParD(level)->QGeom.q27[0],
    //         para->getParD(level)->QGeom.numberOfBCnodes,
    //         para->getParD(level)->omega,
    //         para->getParD(level)->neighborX_SP,
    //         para->getParD(level)->neighborY_SP,
    //         para->getParD(level)->neighborZ_SP,
    //         para->getParD(level)->size_Mat_SP,
    //         para->getParD(level)->evenOrOdd);

    //////////////////////////////////////////////////////////////////////////
    // D E P R E C A T E D
    //////////////////////////////////////////////////////////////////////////
    // the GridGenerator does currently not provide normals!

    //     QSlipGeomDevComp27(
    //         para->getParD(level)->numberofthreads,
    //         para->getParD(level)->d0SP.f[0],
    //         para->getParD(level)->QGeom.k,
    //         para->getParD(level)->QGeom.q27[0],
    //         para->getParD(level)->QGeom.numberOfBCnodes,
    //         para->getParD(level)->omega,
    //         para->getParD(level)->QGeomNormalX.q27[0],
    //         para->getParD(level)->QGeomNormalY.q27[0],
    //         para->getParD(level)->QGeomNormalZ.q27[0],
    //         para->getParD(level)->neighborX_SP,
    //         para->getParD(level)->neighborY_SP,
    //         para->getParD(level)->neighborZ_SP,
    //         para->getParD(level)->size_Mat_SP,
    //         para->getParD(level)->evenOrOdd);

    //     QSlipNormDevComp27(
    //         para->getParD(level)->numberofthreads,
    //         para->getParD(level)->d0SP.f[0],
    //         para->getParD(level)->QGeom.k,
    //         para->getParD(level)->QGeom.q27[0],
    //         para->getParD(level)->QGeom.numberOfBCnodes,
    //         para->getParD(level)->omega,
    //         para->getParD(level)->QGeomNormalX.q27[0],
    //         para->getParD(level)->QGeomNormalY.q27[0],
    //         para->getParD(level)->QGeomNormalZ.q27[0],
    //         para->getParD(level)->neighborX_SP,
    //         para->getParD(level)->neighborY_SP,
    //         para->getParD(level)->neighborZ_SP,
    //         para->getParD(level)->size_Mat_SP,
    //         para->getParD(level)->evenOrOdd);
    }
}

void CudaKernelManager::runStressWallModelKernel(int level){
    if (para->getParD(level)->numberOfStressBCnodes > 0)
    {
        // QStressDevComp27(para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],
        //                 para->getParD(level)->QStress.k,       para->getParD(level)->QStress.kN,
        //                 para->getParD(level)->QStress.q27[0],  para->getParD(level)->numberOfStressBCnodes,
        //                 para->getParD(level)->omega,           para->getParD(level)->turbViscosity,
        //                 para->getParD(level)->vx_SP,           para->getParD(level)->vy_SP,             para->getParD(level)->vy_SP,
        //                 para->getParD(level)->QStress.normalX, para->getParD(level)->QStress.normalY,   para->getParD(level)->QStress.normalZ,
        //                 para->getParD(level)->QStress.Vx,      para->getParD(level)->QStress.Vy,        para->getParD(level)->QStress.Vz,
        //                 para->getParD(level)->QStress.Vx1,     para->getParD(level)->QStress.Vy1,       para->getParD(level)->QStress.Vz1,
        //                 para->getParD(level)->wallModel.samplingOffset, para->getParD(level)->wallModel.z0,
        //                 para->getHasWallModelMonitor(),        para->getParD(level)->wallModel.u_star,
        //                 para->getParD(level)->wallModel.Fx,    para->getParD(level)->wallModel.Fy,      para->getParD(level)->wallModel.Fz,
        //                 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,      para->getParD(level)->neighborZ_SP,
        //                 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);

        BBStressDev27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],
                        para->getParD(level)->QStress.k,       para->getParD(level)->QStress.kN,
                        para->getParD(level)->QStress.q27[0],  para->getParD(level)->numberOfStressBCnodes,
                        para->getParD(level)->vx_SP,           para->getParD(level)->vy_SP,             para->getParD(level)->vy_SP,
                        para->getParD(level)->QStress.normalX, para->getParD(level)->QStress.normalY,   para->getParD(level)->QStress.normalZ,
                        para->getParD(level)->QStress.Vx,      para->getParD(level)->QStress.Vy,        para->getParD(level)->QStress.Vz,
                        para->getParD(level)->QStress.Vx1,     para->getParD(level)->QStress.Vy1,       para->getParD(level)->QStress.Vz1,
                        para->getParD(level)->wallModel.samplingOffset, para->getParD(level)->wallModel.z0,
                        para->getHasWallModelMonitor(),        para->getParD(level)->wallModel.u_star,
                        para->getParD(level)->wallModel.Fx,    para->getParD(level)->wallModel.Fy,      para->getParD(level)->wallModel.Fz,
                        para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,      para->getParD(level)->neighborZ_SP,
                        para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
    }
}


void CudaKernelManager::runSlipBCKernel(int level){
    if (para->getParD(level)->numberOfSlipBCnodes > 0)
    {
        // QSlipDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->QSlip.k,
        //     para->getParD(level)->QSlip.q27[0],
        //     para->getParD(level)->numberOfSlipBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        QSlipDevComp27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->d0SP.f[0],
            para->getParD(level)->QSlip.k,
            para->getParD(level)->QSlip.q27[0],
            para->getParD(level)->numberOfSlipBCnodes,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->turbViscosity,
            para->getUseTurbulentViscosity(),
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd);
    }
}

void CudaKernelManager::runNoSlipBCKernel(int level){
    if (para->getParD(level)->numberOfNoSlipBCnodes > 0)
    {
        // QDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->QWall.k,
        //     para->getParD(level)->QWall.q27[0],
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        // BBDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->QWall.k,
        //     para->getParD(level)->QWall.q27[0],
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        // QDev27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->nx,
        //     para->getParD(level)->ny,
        //     para->getParD(level)->d0SP.f[0],
        //     para->getParD(level)->QWall.k,
        //     para->getParD(level)->QWall.q27[0],
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->numberOfNoSlipBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX_SP,
        //     para->getParD(level)->neighborY_SP,
        //     para->getParD(level)->neighborZ_SP,
        //     para->getParD(level)->size_Mat_SP,
        //     para->getParD(level)->evenOrOdd);

        QDevComp27(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->nx,
            para->getParD(level)->ny,
            para->getParD(level)->d0SP.f[0],
            para->getParD(level)->QWall.k,
            para->getParD(level)->QWall.q27[0],
            para->getParD(level)->numberOfNoSlipBCnodes,
            para->getParD(level)->numberOfNoSlipBCnodes,
            para->getParD(level)->omega,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd);
    }
}

// void CudaKernelManager::runPressureBCKernel(int level){
//     if (para->getParD()->numberOfPressureBCnodes > 0)
//     {
//         // ...
//     }
// }

// void CudaKernelManager::calculateMacroscopicValues(int level)
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









SPtr<CudaKernelManager> CudaKernelManager::make(SPtr<Parameter> parameter)
{
    return SPtr<CudaKernelManager>(new CudaKernelManager(parameter));
}

CudaKernelManager::CudaKernelManager(SPtr<Parameter> parameter)
{
    this->para = parameter;
}

CudaKernelManager::CudaKernelManager(const CudaKernelManager&)
{

}

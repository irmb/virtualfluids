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
//! \author Martin Schoenherr
//=======================================================================================
#include "CudaMemoryManager.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <basics/constants/NumericConstants.h>

#include "Parameter/Parameter.h"
#include "CudaStreamManager.h"
#include "PreCollisionInteractor/Actuator/ActuatorFarm.h"
#include "PreCollisionInteractor/Probes/Probe.h"
#include "PreCollisionInteractor/PrecursorWriter.h"

void CudaMemoryManager::cudaCopyPrint(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityX   , parameter->getParD(lev)->velocityX   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityY   , parameter->getParD(lev)->velocityY   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityZ   , parameter->getParD(lev)->velocityZ   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->rho         , parameter->getParD(lev)->rho         , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->pressure    , parameter->getParD(lev)->pressure    , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));

    if(parameter->getIsBodyForce())
    {
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->forceX_SP   , parameter->getParD(lev)->forceX_SP   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->forceY_SP   , parameter->getParD(lev)->forceY_SP   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->forceZ_SP   , parameter->getParD(lev)->forceZ_SP   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    }

    if(parameter->getUseTurbulentViscosity())
    {
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->turbViscosity   , parameter->getParD(lev)->turbViscosity   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    }
}
void CudaMemoryManager::cudaCopyMeanPrint(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vx_SP_Med   , parameter->getParD(lev)->vx_SP_Med   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vy_SP_Med   , parameter->getParD(lev)->vy_SP_Med   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vz_SP_Med   , parameter->getParD(lev)->vz_SP_Med   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->rho_SP_Med  , parameter->getParD(lev)->rho_SP_Med  , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->press_SP_Med, parameter->getParD(lev)->press_SP_Med, parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaAllocCoord(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordinateX      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordinateY      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordinateZ      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    //Device (spinning ship + uppsala)
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordinateX      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordinateY      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordinateZ      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)parameter->getParH(lev)->memSizeRealLBnodes;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyCoord(int lev)
{
    //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordinateX,  parameter->getParH(lev)->coordinateX,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordinateY,  parameter->getParH(lev)->coordinateY,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordinateZ,  parameter->getParH(lev)->coordinateZ,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeCoord(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->coordinateX   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->coordinateY   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->coordinateZ   ));
}
void CudaMemoryManager::cudaAllocBodyForce(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->forceX_SP      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->forceY_SP      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->forceZ_SP      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->forceX_SP      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->forceY_SP      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->forceZ_SP      ), parameter->getParH(lev)->memSizeRealLBnodes  ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)parameter->getParH(lev)->memSizeRealLBnodes;
    setMemsizeGPU(tmp, false);

}
void CudaMemoryManager::cudaCopyBodyForce(int lev)
{
       //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->forceX_SP,  parameter->getParH(lev)->forceX_SP,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->forceY_SP,  parameter->getParH(lev)->forceY_SP,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->forceZ_SP,  parameter->getParH(lev)->forceZ_SP,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));

}
void CudaMemoryManager::cudaFreeBodyForce(int lev)
{
       checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->forceX_SP   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->forceY_SP   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->forceZ_SP   ));

}
//print
void CudaMemoryManager::cudaCopyDataToHost(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityX   , parameter->getParD(lev)->velocityX   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityY   , parameter->getParD(lev)->velocityY   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityZ   , parameter->getParD(lev)->velocityZ   , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->rho         , parameter->getParD(lev)->rho         , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->pressure    , parameter->getParD(lev)->pressure    , parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
}
//sparse
void CudaMemoryManager::cudaAllocSP(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->typeOfGridNode), parameter->getParH(lev)->memSizeLonglongLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborX     ), parameter->getParH(lev)->memSizeLonglongLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborY     ), parameter->getParH(lev)->memSizeLonglongLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborZ     ), parameter->getParH(lev)->memSizeLonglongLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->rho           ), parameter->getParH(lev)->memSizeRealLBnodes    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityX     ), parameter->getParH(lev)->memSizeRealLBnodes    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityY     ), parameter->getParH(lev)->memSizeRealLBnodes    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityZ     ), parameter->getParH(lev)->memSizeRealLBnodes    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->pressure      ), parameter->getParH(lev)->memSizeRealLBnodes    ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->typeOfGridNode    ), parameter->getParD(lev)->memSizeLonglongLBnodes));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborX         ), parameter->getParD(lev)->memSizeLonglongLBnodes));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborY         ), parameter->getParD(lev)->memSizeLonglongLBnodes));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborZ         ), parameter->getParD(lev)->memSizeLonglongLBnodes));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->rho               ), parameter->getParD(lev)->memSizeRealLBnodes    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityX         ), parameter->getParD(lev)->memSizeRealLBnodes    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityY         ), parameter->getParD(lev)->memSizeRealLBnodes    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityZ         ), parameter->getParD(lev)->memSizeRealLBnodes    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressure          ), parameter->getParD(lev)->memSizeRealLBnodes    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->distributions.f[0]), (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParD(lev)->memSizeRealLBnodes));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 4. * (double)parameter->getParH(lev)->memSizeLonglongLBnodes + 5. * (double)parameter->getParH(lev)->memSizeRealLBnodes + (double)parameter->getD3Qxx() * (double)parameter->getParH(lev)->memSizeRealLBnodes;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopySP(int lev)
{
    //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->typeOfGridNode, parameter->getParH(lev)->typeOfGridNode,  parameter->getParH(lev)->memSizeLonglongLBnodes , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborX     , parameter->getParH(lev)->neighborX     ,  parameter->getParH(lev)->memSizeLonglongLBnodes , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborY     , parameter->getParH(lev)->neighborY     ,  parameter->getParH(lev)->memSizeLonglongLBnodes , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborZ     , parameter->getParH(lev)->neighborZ     ,  parameter->getParH(lev)->memSizeLonglongLBnodes , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->rho           , parameter->getParH(lev)->rho           ,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityX     , parameter->getParH(lev)->velocityX     ,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityY     , parameter->getParH(lev)->velocityY     ,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityZ     , parameter->getParH(lev)->velocityZ     ,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->pressure      , parameter->getParH(lev)->pressure      ,  parameter->getParH(lev)->memSizeRealLBnodes     , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeSP(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->typeOfGridNode ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityX      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityY      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityZ      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->rho            ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->pressure       ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborX      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborY      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborZ      ));
}





















//Velo
void CudaMemoryManager::cudaAllocVeloBC(int lev)
{
    unsigned int mem_size_inflow_Q_k = sizeof(int)*parameter->getParH(lev)->velocityBC.numberOfBCnodes;
    unsigned int mem_size_inflow_Q_q = sizeof(real)*parameter->getParH(lev)->velocityBC.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.q27[0]),  parameter->getD3Qxx()*mem_size_inflow_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.k),                             mem_size_inflow_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.Vx),                            mem_size_inflow_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.Vy),                            mem_size_inflow_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.Vz),                            mem_size_inflow_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.deltaVz),                       mem_size_inflow_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.RhoBC),                         mem_size_inflow_Q_q ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.q27[0]),      parameter->getD3Qxx()*mem_size_inflow_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.k),                                 mem_size_inflow_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.Vx),                                mem_size_inflow_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.Vy),                                mem_size_inflow_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.Vz),                                mem_size_inflow_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.deltaVz),                           mem_size_inflow_Q_q ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_inflow_Q_k + 4. * (double)mem_size_inflow_Q_q + (double)parameter->getD3Qxx() * (double)mem_size_inflow_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyVeloBC(int lev)
{
    unsigned int mem_size_inflow_Q_k = sizeof(int)*parameter->getParH(lev)->velocityBC.numberOfBCnodes;
    unsigned int mem_size_inflow_Q_q = sizeof(real)*parameter->getParH(lev)->velocityBC.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.q27[0],  parameter->getParH(lev)->velocityBC.q27[0], parameter->getD3Qxx()* mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.k,       parameter->getParH(lev)->velocityBC.k,                             mem_size_inflow_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.Vx,      parameter->getParH(lev)->velocityBC.Vx,                            mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.Vy,      parameter->getParH(lev)->velocityBC.Vy,                            mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.Vz,      parameter->getParH(lev)->velocityBC.Vz,                            mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.deltaVz, parameter->getParH(lev)->velocityBC.deltaVz,                       mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));

}

void CudaMemoryManager::cudaFreeVeloBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityBC.q27[0] ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityBC.k      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityBC.Vx     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityBC.Vy     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityBC.Vz     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityBC.deltaVz));
}
//Press
void CudaMemoryManager::cudaAllocOutflowBC(int lev)
{
    unsigned int mem_size_outflow_Q_k = sizeof(int)*parameter->getParH(lev)->outflowBC.numberOfBCnodes;
    unsigned int mem_size_outflow_Q_q = sizeof(real)*parameter->getParH(lev)->outflowBC.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBC.q27[0]), parameter->getD3Qxx()*mem_size_outflow_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBC.k),                            mem_size_outflow_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBC.kN),                           mem_size_outflow_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBC.RhoBC),                        mem_size_outflow_Q_q ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBC.q27[0]),     parameter->getD3Qxx()* mem_size_outflow_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBC.k),                                 mem_size_outflow_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBC.kN),                                mem_size_outflow_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBC.RhoBC),                             mem_size_outflow_Q_q ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_outflow_Q_q + 2. * (double)mem_size_outflow_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_outflow_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyOutflowBC(int lev)
{
    unsigned int mem_size_outflow_Q_k = sizeof(int)*parameter->getParH(lev)->outflowBC.numberOfBCnodes;
    unsigned int mem_size_outflow_Q_q = sizeof(real)*parameter->getParH(lev)->outflowBC.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBC.q27[0],  parameter->getParH(lev)->outflowBC.q27[0], parameter->getD3Qxx()* mem_size_outflow_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBC.k,       parameter->getParH(lev)->outflowBC.k,                             mem_size_outflow_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBC.kN,      parameter->getParH(lev)->outflowBC.kN,                            mem_size_outflow_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBC.RhoBC,   parameter->getParH(lev)->outflowBC.RhoBC,                         mem_size_outflow_Q_q,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeOutflowBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBC.q27[0] ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBC.k      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBC.kN     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBC.RhoBC  ));
}
//No-Slip
void CudaMemoryManager::cudaAllocNoSlipBC(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->noSlipBC.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->noSlipBC.numberOfBCnodes;
    unsigned int mem_size_Q_value  = sizeof(long long)*parameter->getParH(lev)->noSlipBC.numberOfBCnodes; //Geller
    unsigned int mem_size_Q_q_read = sizeof(real)*parameter->getParH(lev)->numberOfNoSlipBCnodesRead;     //Geller

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->noSlipBC.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->noSlipBC.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->noSlipBC.qread),                        mem_size_Q_q_read ));//Geller
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->noSlipBC.valueQ),                       mem_size_Q_value  ));//Geller

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->noSlipBC.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->noSlipBC.k),                                 mem_size_Q_k     ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyNoSlipBC(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->noSlipBC.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->noSlipBC.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->noSlipBC.q27[0], parameter->getParH(lev)->noSlipBC.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->noSlipBC.k,      parameter->getParH(lev)->noSlipBC.k,                             mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeNoSlipBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->noSlipBC.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->noSlipBC.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->noSlipBC.valueQ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->noSlipBC.qread));
}
//Geometrie
void CudaMemoryManager::cudaAllocGeomBC(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->geometryBC.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->geometryBC.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBC.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBC.k),                            mem_size_Q_k      ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBC.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBC.k),                                 mem_size_Q_k     ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyGeomBC(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->geometryBC.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->geometryBC.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBC.q27[0], parameter->getParH(lev)->geometryBC.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBC.k,      parameter->getParH(lev)->geometryBC.k,                             mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeGeomBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBC.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBC.k));
}
//Press
void CudaMemoryManager::cudaAllocPress(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->pressureBC.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->pressureBC.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->pressureBC.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->pressureBC.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->pressureBC.kN),                           mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->pressureBC.RhoBC),                        mem_size_Q_q      ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressureBC.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressureBC.k),                                 mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressureBC.kN),                                mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressureBC.RhoBC),                             mem_size_Q_q     ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 2. * (double)mem_size_Q_k + (double)mem_size_Q_q + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyPress(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->pressureBC.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->pressureBC.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->pressureBC.q27[0], parameter->getParH(lev)->pressureBC.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->pressureBC.k,      parameter->getParH(lev)->pressureBC.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->pressureBC.kN,     parameter->getParH(lev)->pressureBC.kN,                 mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->pressureBC.RhoBC,  parameter->getParH(lev)->pressureBC.RhoBC,              mem_size_Q_q,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreePress(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->pressureBC.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->pressureBC.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->pressureBC.kN));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->pressureBC.RhoBC));
}
//Forcing
void CudaMemoryManager::cudaAllocForcing()
{
    unsigned int mem_size = sizeof(real) * 3;
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->forcingH), mem_size));
    parameter->forcingH[0] = parameter->getForcesDouble()[0];
    parameter->forcingH[1] = parameter->getForcesDouble()[1];
    parameter->forcingH[2] = parameter->getForcesDouble()[2];
    //Device
    checkCudaErrors( cudaMalloc((void**) &parameter->forcingD, mem_size));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyForcingToDevice()
{
    unsigned int mem_size = sizeof(real) * 3;
    checkCudaErrors( cudaMemcpy(parameter->forcingD, parameter->forcingH, mem_size, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyForcingToHost()
{
    unsigned int mem_size = sizeof(real) * 3;
    checkCudaErrors( cudaMemcpy(parameter->forcingH, parameter->forcingD, mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeForcing()
{
    checkCudaErrors( cudaFreeHost(parameter->getForcesHost()));
}

void CudaMemoryManager::cudaAllocLevelForcing(int level)
{
    real fx_t{ 1. }, fy_t{ 1. }, fz_t{ 1. };
    for (int i = 0; i < level; i++) {
        fx_t *= vf::basics::constant::c2o1;
        fy_t *= vf::basics::constant::c2o1;
        fz_t *= vf::basics::constant::c2o1;
    }

    const unsigned int mem_size = sizeof(real) * 3;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(level)->forcing), mem_size));
    parameter->getParH(level)->forcing[0] = parameter->forcingH[0] / fx_t;
    parameter->getParH(level)->forcing[1] = parameter->forcingH[1] / fy_t;
    parameter->getParH(level)->forcing[2] = parameter->forcingH[2] / fz_t;

    //Device
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(level)->forcing, mem_size));
    //////////////////////////////////////////////////////////////////////////
    const double tmp = (double)mem_size;
    setMemsizeGPU(tmp, false);
}

void CudaMemoryManager::cudaCopyLevelForcingToDevice(int level)
{
    unsigned int mem_size = sizeof(real) * 3;
    checkCudaErrors( cudaMemcpy(parameter->getParD(level)->forcing, parameter->getParH(level)->forcing, mem_size, cudaMemcpyHostToDevice));
}

void CudaMemoryManager::cudaFreeLevelForcing(int level)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(level)->forcing));
    checkCudaErrors( cudaFree(parameter->getParD(level)->forcing));
}


//quadric Limiters
void CudaMemoryManager::cudaAllocQuadricLimiters()
{
    unsigned int mem_size = sizeof(real) * 3;
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->quadricLimitersH), mem_size));
    parameter->quadricLimitersH[0] = parameter->getQuadricLimitersDouble()[0];
    parameter->quadricLimitersH[1] = parameter->getQuadricLimitersDouble()[1];
    parameter->quadricLimitersH[2] = parameter->getQuadricLimitersDouble()[2];
    //Device
    checkCudaErrors( cudaMalloc((void**) &parameter->quadricLimitersD, mem_size));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyQuadricLimitersToDevice()
{
    unsigned int mem_size = sizeof(real) * 3;
    checkCudaErrors( cudaMemcpy(parameter->quadricLimitersD, parameter->quadricLimitersH, mem_size, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeQuadricLimiters()
{
    checkCudaErrors( cudaFreeHost(parameter->getQuadricLimitersHost()));
}

//////////////////////////////////////////////////////////////////////////
//Process Neighbors
//  3D domain decomposition
//  X  /////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocProcessNeighborX(int lev, unsigned int processNeighbor)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].index ),                             parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].index ),                             parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].index ),                                 parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].index ),                                 parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].memsizeIndex + (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].memsizeFs +
                 (double)parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].memsizeIndex + (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].memsizeFs;
    setMemsizeGPU(tmp, false);
    //printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void CudaMemoryManager::cudaCopyProcessNeighborXIndex(int lev, unsigned int processNeighbor)
{
    //copy send Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
    //copy recv Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborXFsHD(int lev, unsigned int processNeighbor,
                                                     const unsigned int &memsizeFsRecv)
{
    if (!parameter->getStreamManager()->streamIsRegistered(CudaStreamIndex::SubDomainBorder))
        checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].f[0],
                         parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].f[0],
                         parameter->getD3Qxx() * memsizeFsRecv,
                         cudaMemcpyHostToDevice));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].f[0],
                                         parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].f[0],
                                         parameter->getD3Qxx() * memsizeFsRecv,
                                         cudaMemcpyHostToDevice,
                                         parameter->getStreamManager()->getStream(CudaStreamIndex::SubDomainBorder)));
}
void CudaMemoryManager::cudaCopyProcessNeighborXFsDH(int lev, unsigned int processNeighbor,
                                                     const unsigned int &memsizeFsSend)
{  
    if (!parameter->getStreamManager()->streamIsRegistered(CudaStreamIndex::SubDomainBorder))
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].f[0],
                                    parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].f[0],
                                    parameter->getD3Qxx() * memsizeFsSend,
                                    cudaMemcpyDeviceToHost));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].f[0],
                                         parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].f[0],
                                         parameter->getD3Qxx() * memsizeFsSend,
                                         cudaMemcpyDeviceToHost,
                                         parameter->getStreamManager()->getStream(CudaStreamIndex::SubDomainBorder)));
}
void CudaMemoryManager::cudaFreeProcessNeighborX(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].index ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].f[0]     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].index  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].f[0]     ));
}
//  Y  /////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocProcessNeighborY(int lev, unsigned int processNeighbor)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].index ),                             parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].index ),                             parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].index ),                                 parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].index ),                                 parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].memsizeIndex + (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].memsizeFs +
                 (double)parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].memsizeIndex + (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].memsizeFs;
    setMemsizeGPU(tmp, false);
    //printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void CudaMemoryManager::cudaCopyProcessNeighborYIndex(int lev, unsigned int processNeighbor)
{
    //copy send Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
    //copy recv Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborYFsHD(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsRecv)
{
    if (!parameter->getStreamManager()->streamIsRegistered(CudaStreamIndex::SubDomainBorder))
        checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].f[0],
                                    parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].f[0],
                                    parameter->getD3Qxx() * memsizeFsRecv,
                                    cudaMemcpyHostToDevice));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].f[0],
                                         parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].f[0],
                                         parameter->getD3Qxx() * memsizeFsRecv,
                                         cudaMemcpyHostToDevice,
                                         parameter->getStreamManager()->getStream(CudaStreamIndex::SubDomainBorder)));
}
void CudaMemoryManager::cudaCopyProcessNeighborYFsDH(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsSend)
{
    if (!parameter->getStreamManager()->streamIsRegistered(CudaStreamIndex::SubDomainBorder))
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].f[0],
                                    parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].f[0],
                                    parameter->getD3Qxx() * memsizeFsSend,
                                    cudaMemcpyDeviceToHost));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].f[0],
                                         parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].f[0],
                                         parameter->getD3Qxx() * memsizeFsSend,
                                         cudaMemcpyDeviceToHost, 
                                         parameter->getStreamManager()->getStream(CudaStreamIndex::SubDomainBorder)));
}
void CudaMemoryManager::cudaFreeProcessNeighborY(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].index ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].f[0]     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].index  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].f[0]     ));
}
//  Z  /////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocProcessNeighborZ(int lev, unsigned int processNeighbor)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].index ),                             parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].index ),                             parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].index ),                                 parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].index ),                                 parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].memsizeIndex + (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].memsizeFs +
                 (double)parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].memsizeIndex + (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].memsizeFs;
    setMemsizeGPU(tmp, false);
    //printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void CudaMemoryManager::cudaCopyProcessNeighborZIndex(int lev, unsigned int processNeighbor)
{
    //copy send Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
    //copy recv Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborZFsHD(int lev, unsigned int processNeighbor,
                                                     const unsigned int &memsizeFsRecv)
{
    if (!parameter->getStreamManager()->streamIsRegistered(CudaStreamIndex::SubDomainBorder))
        checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].f[0],
                                    parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].f[0],
                                    parameter->getD3Qxx() * memsizeFsRecv,
                                    cudaMemcpyHostToDevice));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].f[0],
                                         parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].f[0],
                                         parameter->getD3Qxx() * memsizeFsRecv,
                                         cudaMemcpyHostToDevice,
                                         parameter->getStreamManager()->getStream(CudaStreamIndex::SubDomainBorder)));
}
void CudaMemoryManager::cudaCopyProcessNeighborZFsDH(int lev, unsigned int processNeighbor,
                                                     const unsigned int &memsizeFsSend)
{
    if (!parameter->getStreamManager()->streamIsRegistered(CudaStreamIndex::SubDomainBorder))
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].f[0],
                                    parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].f[0],
                                    parameter->getD3Qxx() * memsizeFsSend,
                                    cudaMemcpyDeviceToHost));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].f[0],
                                         parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].f[0],
                                         parameter->getD3Qxx() * memsizeFsSend,
                                         cudaMemcpyDeviceToHost,
                                         parameter->getStreamManager()->getStream(CudaStreamIndex::SubDomainBorder)));
}
void CudaMemoryManager::cudaFreeProcessNeighborZ(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].index ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].f[0]     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].index  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].f[0]     ));
}



void CudaMemoryManager::cudaAllocNeighborWSB(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborInverse    ), parameter->getParH(lev)->memSizeLonglongLBnodes    ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborInverse        ), parameter->getParD(lev)->memSizeLonglongLBnodes    ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->memSizeLonglongLBnodes;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyNeighborWSB(int lev)
{
    //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborInverse,  parameter->getParH(lev)->neighborInverse,  parameter->getParH(lev)->memSizeLonglongLBnodes     , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeNeighborWSB(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborInverse));
}
//turbulent viscosity
void CudaMemoryManager::cudaAllocTurbulentViscosity(int lev)
{
    //Host
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->turbViscosity), parameter->getParH(lev)->memSizeRealLBnodes));

    //Device
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->turbViscosity), parameter->getParD(lev)->memSizeRealLBnodes));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->memSizeRealLBnodes;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyTurbulentViscosityHD(int lev)
{
    //copy host to device
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->turbViscosity, parameter->getParH(lev)->turbViscosity, parameter->getParH(lev)->memSizeRealLBnodes, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyTurbulentViscosityDH(int lev)
{
    //copy device to host
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->turbViscosity, parameter->getParD(lev)->turbViscosity, parameter->getParH(lev)->memSizeRealLBnodes, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeTurbulentViscosity(int lev)
{
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->turbViscosity));
}
//turbulence intensity
void CudaMemoryManager::cudaAllocTurbulenceIntensity(int lev, uint size)
{
    uint mem_size = sizeof(real) * size;
    // Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vxx        ), mem_size));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vyy        ), mem_size));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vzz        ), mem_size));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vxy        ), mem_size));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vxz        ), mem_size));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vyz        ), mem_size));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vx_mean    ), mem_size));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vy_mean    ), mem_size));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vz_mean    ), mem_size));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vxx            ), mem_size));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vyy            ), mem_size));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vzz            ), mem_size));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vxy            ), mem_size));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vxz            ), mem_size));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vyz            ), mem_size));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vx_mean        ), mem_size));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vy_mean        ), mem_size));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vz_mean        ), mem_size));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 9. * (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyTurbulenceIntensityHD(int lev, uint size)
{
    uint mem_size = sizeof(real) * size;
    //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vxx    ,  parameter->getParH(lev)->vxx    ,  mem_size , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vyy    ,  parameter->getParH(lev)->vyy    ,  mem_size , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vzz    ,  parameter->getParH(lev)->vzz    ,  mem_size , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vxy    ,  parameter->getParH(lev)->vxy    ,  mem_size , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vxz    ,  parameter->getParH(lev)->vxz    ,  mem_size , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vyz    ,  parameter->getParH(lev)->vyz    ,  mem_size , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vx_mean,  parameter->getParH(lev)->vx_mean,  mem_size , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vy_mean,  parameter->getParH(lev)->vy_mean,  mem_size , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vz_mean,  parameter->getParH(lev)->vz_mean,  mem_size , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyTurbulenceIntensityDH(int lev, uint size)
{
    uint mem_size = sizeof(real) * size;
    //copy device to host
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vxx    ,  parameter->getParD(lev)->vxx    ,  mem_size , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vyy    ,  parameter->getParD(lev)->vyy    ,  mem_size , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vzz    ,  parameter->getParD(lev)->vzz    ,  mem_size , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vxy    ,  parameter->getParD(lev)->vxy    ,  mem_size , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vxz    ,  parameter->getParD(lev)->vxz    ,  mem_size , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vyz    ,  parameter->getParD(lev)->vyz    ,  mem_size , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vx_mean,  parameter->getParD(lev)->vx_mean,  mem_size , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vy_mean,  parameter->getParD(lev)->vy_mean,  mem_size , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vz_mean,  parameter->getParD(lev)->vz_mean,  mem_size , cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeTurbulenceIntensity(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vxx     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vyy     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vzz     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vxy     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vxz     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vyz     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vx_mean ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vy_mean ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vz_mean ));
}
//mean
void CudaMemoryManager::cudaAllocMeanSP(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->rho_SP_Med      ), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vx_SP_Med       ), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vy_SP_Med       ), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vz_SP_Med       ), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->press_SP_Med    ), parameter->getParH(lev)->memSizeRealLBnodes));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->rho_SP_Med          ), parameter->getParD(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vx_SP_Med           ), parameter->getParD(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vy_SP_Med           ), parameter->getParD(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vz_SP_Med           ), parameter->getParD(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->press_SP_Med        ), parameter->getParD(lev)->memSizeRealLBnodes));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 5. * (double)parameter->getParH(lev)->memSizeRealLBnodes;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyMeanSP(int lev)
{
    //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->rho_SP_Med  ,  parameter->getParH(lev)->rho_SP_Med  ,  parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vx_SP_Med   ,  parameter->getParH(lev)->vx_SP_Med   ,  parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vy_SP_Med   ,  parameter->getParH(lev)->vy_SP_Med   ,  parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vz_SP_Med   ,  parameter->getParH(lev)->vz_SP_Med   ,  parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->press_SP_Med,  parameter->getParH(lev)->press_SP_Med,  parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeMeanSP(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vx_SP_Med   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vy_SP_Med   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vz_SP_Med   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->rho_SP_Med  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->press_SP_Med));
}
void CudaMemoryManager::cudaAllocMeanOut(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->rho_SP_Med_Out      ), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vx_SP_Med_Out       ), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vy_SP_Med_Out       ), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vz_SP_Med_Out       ), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->press_SP_Med_Out    ), parameter->getParH(lev)->memSizeRealLBnodes));
}
void CudaMemoryManager::cudaFreeMeanOut(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vx_SP_Med_Out   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vy_SP_Med_Out   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vz_SP_Med_Out   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->rho_SP_Med_Out  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->press_SP_Med_Out));
}
//Interface CF
void CudaMemoryManager::cudaAllocInterfaceCF(int lev)
{
    uint mem_size_kCF = sizeof(uint) * parameter->getParH(lev)->coarseToFine.numberOfCells;
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coarseToFine.coarseCellIndices), mem_size_kCF  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coarseToFine.fineCellIndices), mem_size_kCF  ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coarseToFine.coarseCellIndices), mem_size_kCF  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coarseToFine.fineCellIndices), mem_size_kCF  ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 2. * (double)mem_size_kCF;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInterfaceCF(int lev)
{
    uint mem_size_kCF = sizeof(uint) * parameter->getParH(lev)->coarseToFine.numberOfCells;

    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->coarseToFine.coarseCellIndices, parameter->getParH(lev)->coarseToFine.coarseCellIndices, mem_size_kCF, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->coarseToFine.fineCellIndices, parameter->getParH(lev)->coarseToFine.fineCellIndices, mem_size_kCF, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeInterfaceCF(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->coarseToFine.coarseCellIndices));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->coarseToFine.fineCellIndices));
}
//Interface FC
void CudaMemoryManager::cudaAllocInterfaceFC(int lev)
{
    uint mem_size_kFC = sizeof(uint) * parameter->getParH(lev)->fineToCoarse.numberOfCells;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->fineToCoarse.fineCellIndices), mem_size_kFC  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->fineToCoarse.coarseCellIndices), mem_size_kFC  ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->fineToCoarse.fineCellIndices), mem_size_kFC  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->fineToCoarse.coarseCellIndices), mem_size_kFC  ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 2. * (double)mem_size_kFC;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInterfaceFC(int lev)
{
    uint mem_size_kFC = sizeof(uint) * parameter->getParH(lev)->fineToCoarse.numberOfCells;

    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->fineToCoarse.fineCellIndices, parameter->getParH(lev)->fineToCoarse.fineCellIndices, mem_size_kFC, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->fineToCoarse.coarseCellIndices, parameter->getParH(lev)->fineToCoarse.coarseCellIndices, mem_size_kFC, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCheckInterfaceFCBulk(int lev)
{
    // only use for testing!
    size_t memsize = sizeof(uint) * parameter->getParH(lev)->fineToCoarseBulk.numberOfCells;
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->fineToCoarseBulk.coarseCellIndices, parameter->getParH(lev)->fineToCoarseBulk.coarseCellIndices, memsize, cudaMemcpyDeviceToDevice));
    for (uint i = 0; i < parameter->getParH(lev)->fineToCoarseBulk.numberOfCells; i++)
        printf("%d %d\n", i, parameter->getParH(lev)->fineToCoarseBulk.coarseCellIndices[i]);
}
void CudaMemoryManager::cudaFreeInterfaceFC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->fineToCoarse.fineCellIndices));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->fineToCoarse.coarseCellIndices));
}
//Interface Offset CF
void CudaMemoryManager::cudaAllocInterfaceOffCF(int lev)
{
    uint mem_size_kCF_off = sizeof(real) * parameter->getParH(lev)->coarseToFine.numberOfCells;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborCoarseToFine.x), mem_size_kCF_off  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborCoarseToFine.y), mem_size_kCF_off  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborCoarseToFine.z), mem_size_kCF_off  ));
    getLastCudaError("Allocate host memory");
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborCoarseToFine.x), mem_size_kCF_off  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborCoarseToFine.y), mem_size_kCF_off  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborCoarseToFine.z), mem_size_kCF_off  ));
    getLastCudaError("Allocate device memory");
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)mem_size_kCF_off;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInterfaceOffCF(int lev)
{
    uint mem_size_kCF_off = sizeof(real) * parameter->getParH(lev)->coarseToFine.numberOfCells;

    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->neighborCoarseToFine.x, parameter->getParH(lev)->neighborCoarseToFine.x, mem_size_kCF_off, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->neighborCoarseToFine.y, parameter->getParH(lev)->neighborCoarseToFine.y, mem_size_kCF_off, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->neighborCoarseToFine.z, parameter->getParH(lev)->neighborCoarseToFine.z, mem_size_kCF_off, cudaMemcpyHostToDevice));
    getLastCudaError("Copy host memory to device");
}
void CudaMemoryManager::cudaFreeInterfaceOffCF(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborCoarseToFine.x));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborCoarseToFine.y));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborCoarseToFine.z));
}
//Interface Offset FC
void CudaMemoryManager::cudaAllocInterfaceOffFC(int lev)
{
    uint mem_size_kFC_off = sizeof(real) * parameter->getParH(lev)->fineToCoarse.numberOfCells;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborFineToCoarse.x), mem_size_kFC_off  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborFineToCoarse.y), mem_size_kFC_off  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborFineToCoarse.z), mem_size_kFC_off  ));
    getLastCudaError("Allocate host memory");
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborFineToCoarse.x), mem_size_kFC_off  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborFineToCoarse.y), mem_size_kFC_off  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborFineToCoarse.z), mem_size_kFC_off  ));
    getLastCudaError("Allocate device memory");
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)mem_size_kFC_off;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInterfaceOffFC(int lev)
{
    uint mem_size_kFC_off = sizeof(real) * parameter->getParH(lev)->fineToCoarse.numberOfCells;

    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->neighborFineToCoarse.x, parameter->getParH(lev)->neighborFineToCoarse.x, mem_size_kFC_off, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->neighborFineToCoarse.y, parameter->getParH(lev)->neighborFineToCoarse.y, mem_size_kFC_off, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->neighborFineToCoarse.z, parameter->getParH(lev)->neighborFineToCoarse.z, mem_size_kFC_off, cudaMemcpyHostToDevice));
    getLastCudaError("Copy host memory to device");
}
void CudaMemoryManager::cudaFreeInterfaceOffFC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborFineToCoarse.x));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborFineToCoarse.y));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborFineToCoarse.z));
}

//Inlet
void CudaMemoryManager::cudaAllocInlet(int lev)
{
    unsigned int mem_size_inlet_Q_k = sizeof(int)*parameter->getParH(lev)->QInlet.numberOfBCnodes;
    unsigned int mem_size_inlet_Q_q = sizeof(real)*parameter->getParH(lev)->QInlet.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInlet.q27[0]), parameter->getD3Qxx()*mem_size_inlet_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInlet.k),                 mem_size_inlet_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInlet.kN),                mem_size_inlet_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInlet.RhoBC),             mem_size_inlet_Q_q ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInlet.q27[0]), parameter->getD3Qxx()* mem_size_inlet_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInlet.k),                      mem_size_inlet_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInlet.kN),                     mem_size_inlet_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInlet.RhoBC),                  mem_size_inlet_Q_q ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 2. * (double)mem_size_inlet_Q_k + (double)mem_size_inlet_Q_q + (double)parameter->getD3Qxx()*(double)mem_size_inlet_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInlet(int lev)
{
    unsigned int mem_size_inlet_Q_k = sizeof(int)*parameter->getParH(lev)->QInlet.numberOfBCnodes;
    unsigned int mem_size_inlet_Q_q = sizeof(real)*parameter->getParH(lev)->QInlet.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInlet.q27[0],  parameter->getParH(lev)->QInlet.q27[0], parameter->getD3Qxx()* mem_size_inlet_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInlet.k,       parameter->getParH(lev)->QInlet.k,                  mem_size_inlet_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInlet.kN,      parameter->getParH(lev)->QInlet.kN,                 mem_size_inlet_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInlet.RhoBC,   parameter->getParH(lev)->QInlet.RhoBC,              mem_size_inlet_Q_q,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeInlet(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInlet.q27[0] ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInlet.k      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInlet.kN     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInlet.RhoBC  ));
}
//Outlet
void CudaMemoryManager::cudaAllocOutlet(int lev)
{
    unsigned int mem_size_outlet_Q_k = sizeof(int)*parameter->getParH(lev)->QOutlet.numberOfBCnodes;
    unsigned int mem_size_outlet_Q_q = sizeof(real)*parameter->getParH(lev)->QOutlet.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutlet.q27[0]), parameter->getD3Qxx()*mem_size_outlet_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutlet.k),                 mem_size_outlet_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutlet.kN),                mem_size_outlet_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutlet.RhoBC),             mem_size_outlet_Q_q ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutlet.q27[0]),     parameter->getD3Qxx()* mem_size_outlet_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutlet.k),                      mem_size_outlet_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutlet.kN),                     mem_size_outlet_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutlet.RhoBC),                  mem_size_outlet_Q_q ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 2. * (double)mem_size_outlet_Q_k + (double)mem_size_outlet_Q_q + (double)parameter->getD3Qxx()*(double)mem_size_outlet_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyOutlet(int lev)
{
    unsigned int mem_size_outlet_Q_k = sizeof(int)*parameter->getParH(lev)->QOutlet.numberOfBCnodes;
    unsigned int mem_size_outlet_Q_q = sizeof(real)*parameter->getParH(lev)->QOutlet.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutlet.q27[0],  parameter->getParH(lev)->QOutlet.q27[0], parameter->getD3Qxx()* mem_size_outlet_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutlet.k,       parameter->getParH(lev)->QOutlet.k,                  mem_size_outlet_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutlet.kN,      parameter->getParH(lev)->QOutlet.kN,                 mem_size_outlet_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutlet.RhoBC,   parameter->getParH(lev)->QOutlet.RhoBC,              mem_size_outlet_Q_q,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeOutlet(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutlet.q27[0] ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutlet.k      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutlet.kN     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutlet.RhoBC  ));
}



//Geometrie inkl. Values
void CudaMemoryManager::cudaAllocGeomValuesBC(int lev)
{
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->geometryBC.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBC.Vx),  mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBC.Vy),  mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBC.Vz),  mem_size_Q_q ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBC.Vx),      mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBC.Vy),      mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBC.Vz),      mem_size_Q_q ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyGeomValuesBC(int lev)
{
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->geometryBC.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBC.Vx, parameter->getParH(lev)->geometryBC.Vx,  mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBC.Vy, parameter->getParH(lev)->geometryBC.Vy,  mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBC.Vz, parameter->getParH(lev)->geometryBC.Vz,  mem_size_Q_q,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeGeomValuesBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBC.Vx));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBC.Vy));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBC.Vz));
}
//Geometrie inkl. Normale f¸r Slip
void CudaMemoryManager::cudaAllocGeomNormals(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->geometryBCnormalX.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->geometryBCnormalX.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBCnormalX.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBCnormalX.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBCnormalY.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBCnormalY.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBCnormalZ.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBCnormalZ.k),                            mem_size_Q_k      ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBCnormalX.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBCnormalX.k),                                 mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBCnormalY.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBCnormalY.k),                                 mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBCnormalZ.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBCnormalZ.k),                                 mem_size_Q_k     ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 3.0 * (double)mem_size_Q_k + 3.0 * (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyGeomNormals(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->geometryBCnormalX.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->geometryBCnormalX.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBCnormalX.q27[0], parameter->getParH(lev)->geometryBCnormalX.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBCnormalX.k,      parameter->getParH(lev)->geometryBCnormalX.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBCnormalY.q27[0], parameter->getParH(lev)->geometryBCnormalY.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBCnormalY.k,      parameter->getParH(lev)->geometryBCnormalY.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBCnormalZ.q27[0], parameter->getParH(lev)->geometryBCnormalZ.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBCnormalZ.k,      parameter->getParH(lev)->geometryBCnormalZ.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeGeomNormals(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBCnormalX.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBCnormalX.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBCnormalY.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBCnormalY.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBCnormalZ.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geometryBCnormalZ.k));
}
//Geometrie inkl. Normale f¸r Inflow
void CudaMemoryManager::cudaAllocInflowNormals(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->inflowBCnormalX.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->inflowBCnormalX.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->inflowBCnormalX.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->inflowBCnormalX.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->inflowBCnormalY.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->inflowBCnormalY.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->inflowBCnormalZ.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->inflowBCnormalZ.k),                            mem_size_Q_k      ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->inflowBCnormalX.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->inflowBCnormalX.k),                                 mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->inflowBCnormalY.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->inflowBCnormalY.k),                                 mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->inflowBCnormalZ.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->inflowBCnormalZ.k),                                 mem_size_Q_k     ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 3.0 * (double)mem_size_Q_k + 3.0 * (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInflowNormals(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->inflowBCnormalX.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->inflowBCnormalX.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->inflowBCnormalX.q27[0], parameter->getParH(lev)->inflowBCnormalX.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->inflowBCnormalX.k,      parameter->getParH(lev)->inflowBCnormalX.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->inflowBCnormalY.q27[0], parameter->getParH(lev)->inflowBCnormalY.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->inflowBCnormalY.k,      parameter->getParH(lev)->inflowBCnormalY.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->inflowBCnormalZ.q27[0], parameter->getParH(lev)->inflowBCnormalZ.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->inflowBCnormalZ.k,      parameter->getParH(lev)->inflowBCnormalZ.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeInflowNormals(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->inflowBCnormalX.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->inflowBCnormalX.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->inflowBCnormalY.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->inflowBCnormalY.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->inflowBCnormalZ.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->inflowBCnormalZ.k));
}
//Geometrie inkl. Normale f¸r Outflow
void CudaMemoryManager::cudaAllocOutflowNormals(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->outflowBCnormalX.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->outflowBCnormalX.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBCnormalX.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBCnormalX.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBCnormalY.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBCnormalY.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBCnormalZ.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBCnormalZ.k),                            mem_size_Q_k      ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBCnormalX.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBCnormalX.k),                                 mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBCnormalY.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBCnormalY.k),                                 mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBCnormalZ.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBCnormalZ.k),                                 mem_size_Q_k     ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 3.0 * (double)mem_size_Q_k + 3.0 * (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyOutflowNormals(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->outflowBCnormalX.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->outflowBCnormalX.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBCnormalX.q27[0], parameter->getParH(lev)->outflowBCnormalX.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBCnormalX.k,      parameter->getParH(lev)->outflowBCnormalX.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBCnormalY.q27[0], parameter->getParH(lev)->outflowBCnormalY.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBCnormalY.k,      parameter->getParH(lev)->outflowBCnormalY.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBCnormalZ.q27[0], parameter->getParH(lev)->outflowBCnormalZ.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBCnormalZ.k,      parameter->getParH(lev)->outflowBCnormalZ.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeOutflowNormals(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBCnormalX.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBCnormalX.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBCnormalY.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBCnormalY.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBCnormalZ.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->outflowBCnormalZ.k));
}
//Slip
void CudaMemoryManager::cudaAllocSlipBC(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->slipBC.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->slipBC.numberOfBCnodes;
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->slipBC.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->slipBC.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->slipBC.normalX),                      mem_size_Q_q    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->slipBC.normalY),                      mem_size_Q_q    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->slipBC.normalZ),                      mem_size_Q_q    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->slipBC.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->slipBC.k),                                 mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->slipBC.normalX),                           mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->slipBC.normalY),                           mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->slipBC.normalZ),                           mem_size_Q_q     ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q + 3.0*(double)mem_size_Q_q;;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopySlipBC(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->slipBC.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->slipBC.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->slipBC.q27[0], parameter->getParH(lev)->slipBC.q27[0],  parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->slipBC.k,      parameter->getParH(lev)->slipBC.k,                              mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->slipBC.normalX, parameter->getParH(lev)->slipBC.normalX,                       mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->slipBC.normalY, parameter->getParH(lev)->slipBC.normalY,                       mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->slipBC.normalZ, parameter->getParH(lev)->slipBC.normalZ,                       mem_size_Q_q,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeSlipBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->slipBC.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->slipBC.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->slipBC.normalX));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->slipBC.normalY));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->slipBC.normalZ));
}
//Stress
void CudaMemoryManager::cudaAllocStressBC(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->stressBC.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->stressBC.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.k),                            mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.kN),                           mem_size_Q_k     ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.normalX),                      mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.normalY),                      mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.normalZ),                      mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.Vx),                           mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.Vy),                           mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.Vz),                           mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.Vx1),                          mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.Vy1),                          mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->stressBC.Vz1),                          mem_size_Q_q      ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.k),                                 mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.kN),                                mem_size_Q_k    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.normalX),                           mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.normalY),                           mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.normalZ),                           mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.Vx),                                mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.Vy),                                mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.Vz),                                mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.Vx1),                               mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.Vy1),                               mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->stressBC.Vz1),                               mem_size_Q_q     ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 2*(double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q + 9.0*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyStressBC(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->stressBC.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->stressBC.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.q27[0],  parameter->getParH(lev)->stressBC.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.k,       parameter->getParH(lev)->stressBC.k,                             mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.kN,      parameter->getParH(lev)->stressBC.kN,                            mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.normalX, parameter->getParH(lev)->stressBC.normalX,                       mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.normalY, parameter->getParH(lev)->stressBC.normalY,                       mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.normalZ, parameter->getParH(lev)->stressBC.normalZ,                       mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.Vx,      parameter->getParH(lev)->stressBC.Vx,                            mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.Vy,      parameter->getParH(lev)->stressBC.Vy,                            mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.Vz,      parameter->getParH(lev)->stressBC.Vz,                            mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.Vx1,     parameter->getParH(lev)->stressBC.Vx1,                           mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.Vy1,     parameter->getParH(lev)->stressBC.Vy1,                           mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->stressBC.Vz1,     parameter->getParH(lev)->stressBC.Vz1,                           mem_size_Q_q,       cudaMemcpyHostToDevice));

}
void CudaMemoryManager::cudaFreeStressBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.kN));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.normalX));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.normalY));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.normalZ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.Vx));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.Vy));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.Vz));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.Vx1));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.Vy1));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->stressBC.Vz1));
}
// Wall model
void CudaMemoryManager::cudaAllocWallModel(int lev, bool hasWallModelMonitor)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->stressBC.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->stressBC.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->wallModel.samplingOffset),  mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->wallModel.z0),              mem_size_Q_q      ));
    if(hasWallModelMonitor)
    {
        checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->wallModel.u_star),      mem_size_Q_q      ));
        checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->wallModel.Fx),          mem_size_Q_q      ));
        checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->wallModel.Fy),          mem_size_Q_q      ));
        checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->wallModel.Fz),          mem_size_Q_q      ));
    }

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->wallModel.samplingOffset),  mem_size_Q_k));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->wallModel.z0),  mem_size_Q_q));
    if(hasWallModelMonitor)
    {
        checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->wallModel.u_star),      mem_size_Q_q      ));
        checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->wallModel.Fx),          mem_size_Q_q      ));
        checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->wallModel.Fy),          mem_size_Q_q      ));
        checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->wallModel.Fz),          mem_size_Q_q      ));
    }

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_Q_k + (double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyWallModel(int lev, bool hasWallModelMonitor)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->stressBC.numberOfBCnodes;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->stressBC.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->wallModel.samplingOffset,  parameter->getParH(lev)->wallModel.samplingOffset,  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->wallModel.z0,              parameter->getParH(lev)->wallModel.z0,              mem_size_Q_q,       cudaMemcpyHostToDevice));
    if(hasWallModelMonitor)
    {
        checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->wallModel.u_star,          parameter->getParH(lev)->wallModel.u_star,          mem_size_Q_k,       cudaMemcpyHostToDevice));
        checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->wallModel.Fx,              parameter->getParH(lev)->wallModel.Fx,              mem_size_Q_q,       cudaMemcpyHostToDevice));
        checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->wallModel.Fy,              parameter->getParH(lev)->wallModel.Fy,              mem_size_Q_q,       cudaMemcpyHostToDevice));
        checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->wallModel.Fz,              parameter->getParH(lev)->wallModel.Fz,              mem_size_Q_q,       cudaMemcpyHostToDevice));
    }
}
void CudaMemoryManager::cudaFreeWallModel(int lev, bool hasWallModelMonitor)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->wallModel.samplingOffset));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->wallModel.z0));
    if(hasWallModelMonitor)
    {
        checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->wallModel.u_star));
        checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->wallModel.Fx));
        checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->wallModel.Fy));
        checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->wallModel.Fz));
    }
}


//Precursor BC
void CudaMemoryManager::cudaAllocPrecursorBC(int lev)
{   
    uint memSizeQInt = parameter->getParH(lev)->precursorBC.numberOfBCnodes*sizeof(int);
    uint memSizeQUint = parameter->getParH(lev)->precursorBC.numberOfBCnodes*sizeof(uint);
    uint memSizeQReal = parameter->getParH(lev)->precursorBC.numberOfBCnodes*sizeof(real);

    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.k, memSizeQInt));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.q27[0], parameter->getD3Qxx()*memSizeQReal));
    

    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.planeNeighbor0PP, memSizeQUint));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.planeNeighbor0PM, memSizeQUint));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.planeNeighbor0MP, memSizeQUint));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.planeNeighbor0MM, memSizeQUint));

    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.weights0PP, memSizeQReal));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.weights0PM, memSizeQReal));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.weights0MP, memSizeQReal));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.weights0MM, memSizeQReal));

    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.k, memSizeQInt));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.q27[0], parameter->getD3Qxx()*memSizeQReal));

    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.planeNeighbor0PP, memSizeQUint));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.planeNeighbor0PM, memSizeQUint));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.planeNeighbor0MP, memSizeQUint));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.planeNeighbor0MM, memSizeQUint));

    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.weights0PP, memSizeQReal));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.weights0PM, memSizeQReal));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.weights0MP, memSizeQReal));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.weights0MM, memSizeQReal));

    real memSize = memSizeQInt+4*memSizeQUint+(4+parameter->getD3Qxx())*memSizeQReal;
    setMemsizeGPU(memSize, false);

}


void CudaMemoryManager::cudaAllocPrecursorData(int lev)
{
    size_t size = parameter->getParH(lev)->precursorBC.numberOfPrecursorNodes*sizeof(real)*parameter->getParH(lev)->precursorBC.numberOfQuantities;

    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.last, size));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.current, size));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->precursorBC.next, size));

    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.last, size));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.current, size));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->precursorBC.next, size));
    setMemsizeGPU(3*size, false);
}


void CudaMemoryManager::cudaCopyPrecursorBC(int lev)
{
    uint memSizeQInt = parameter->getParH(lev)->precursorBC.numberOfBCnodes*sizeof(int);
    uint memSizeQUint = parameter->getParH(lev)->precursorBC.numberOfBCnodes*sizeof(uint);
    uint memSizeQReal = parameter->getParH(lev)->precursorBC.numberOfBCnodes*sizeof(real);

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.k, parameter->getParH(lev)->precursorBC.k, memSizeQInt, cudaMemcpyHostToDevice));

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.q27[0], parameter->getParH(lev)->precursorBC.q27[0], memSizeQReal*parameter->getD3Qxx(), cudaMemcpyHostToDevice));

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.planeNeighbor0PP, parameter->getParH(lev)->precursorBC.planeNeighbor0PP, memSizeQUint, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.planeNeighbor0PM, parameter->getParH(lev)->precursorBC.planeNeighbor0PM, memSizeQUint, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.planeNeighbor0MP, parameter->getParH(lev)->precursorBC.planeNeighbor0MP, memSizeQUint, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.planeNeighbor0MM, parameter->getParH(lev)->precursorBC.planeNeighbor0MM, memSizeQUint, cudaMemcpyHostToDevice));

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.weights0PP, parameter->getParH(lev)->precursorBC.weights0PP, memSizeQReal, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.weights0PM, parameter->getParH(lev)->precursorBC.weights0PM, memSizeQReal, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.weights0MP, parameter->getParH(lev)->precursorBC.weights0MP, memSizeQReal, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->precursorBC.weights0MM, parameter->getParH(lev)->precursorBC.weights0MM, memSizeQReal, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyPrecursorData(int lev)
{
    auto precurser = &parameter->getParH(lev)->precursorBC;
    auto precurserStream = parameter->getStreamManager()->getStream(CudaStreamIndex::Precursor, precurser->streamIndex);
    size_t memSize = precurser->numberOfPrecursorNodes*sizeof(real)*precurser->numberOfQuantities;
    checkCudaErrors( cudaStreamSynchronize(precurserStream) );
    checkCudaErrors( cudaMemcpyAsync(parameter->getParD(lev)->precursorBC.next, precurser->next, memSize, cudaMemcpyHostToDevice, precurserStream) );
}


void CudaMemoryManager::cudaFreePrecursorBC(int lev)
{
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.k));

    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.q27[0]));

    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.planeNeighbor0PP));
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.planeNeighbor0PM));
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.planeNeighbor0MP));
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.planeNeighbor0MM));

    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.weights0PP));
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.weights0PM));
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.weights0MP));
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.weights0MM));

    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.k));

    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.q27[0]));

    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.planeNeighbor0PP));
    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.planeNeighbor0PM));
    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.planeNeighbor0MP));
    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.planeNeighbor0MM));

    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.weights0PP));
    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.weights0PM));
    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.weights0MP));
    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.weights0MM));
}

void CudaMemoryManager::cudaFreePrecursorData(int lev)
{
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.last));
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.current));
    checkCudaErrors( cudaFreeHost( parameter->getParH(lev)->precursorBC.next));

    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.last));
    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.current));
    checkCudaErrors( cudaFree( parameter->getParD(lev)->precursorBC.next));
}
//Test roundoff error
void CudaMemoryManager::cudaAllocTestRE(int lev, unsigned int size)
{
    unsigned int mem_size = sizeof(real)*size;
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->kDistTestRE.f[0]), (1+parameter->getD3Qxx())*mem_size));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->kDistTestRE.f[0]), (1+parameter->getD3Qxx())*mem_size));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)((1. + parameter->getD3Qxx()) * mem_size);
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyTestREtoDevice(int lev, unsigned int size)
{
    unsigned int mem_size = sizeof(real)*size;
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->kDistTestRE.f[0], parameter->getParH(lev)->kDistTestRE.f[0], (1+parameter->getD3Qxx())*mem_size, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyTestREtoHost(int lev, unsigned int size)
{
    unsigned int mem_size = sizeof(real)*size;
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->kDistTestRE.f[0], parameter->getParD(lev)->kDistTestRE.f[0], (1+parameter->getD3Qxx())*mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeTestRE(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->kDistTestRE.f[0]));
}
//PressX0 = X-inflow
void CudaMemoryManager::cudaAllocPressX0(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QpressX0.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QpressX0.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX0.q27[0]), parameter->getD3Qxx()*mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX0.k),                 mem_size_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX0.kN),                mem_size_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX0.Vx),                mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX0.Vy),                mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX0.Vz),                mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX0.RhoBC),             mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX0.deltaVz),           mem_size_Q_q ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX0.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX0.k),                      mem_size_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX0.kN),                     mem_size_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX0.Vx),                     mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX0.Vy),                     mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX0.Vz),                     mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX0.RhoBC),                  mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX0.deltaVz),                mem_size_Q_q ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 2. * (double)mem_size_Q_k + 5. * (double)mem_size_Q_q + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyPressX0(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QpressX0.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QpressX0.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX0.q27[0],  parameter->getParH(lev)->QpressX0.q27[0], parameter->getD3Qxx()* mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX0.k,       parameter->getParH(lev)->QpressX0.k,                  mem_size_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX0.kN,      parameter->getParH(lev)->QpressX0.kN,                 mem_size_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX0.Vx,      parameter->getParH(lev)->QpressX0.Vx,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX0.Vy,      parameter->getParH(lev)->QpressX0.Vy,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX0.Vz,      parameter->getParH(lev)->QpressX0.Vz,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX0.RhoBC,   parameter->getParH(lev)->QpressX0.RhoBC,              mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX0.deltaVz, parameter->getParH(lev)->QpressX0.deltaVz,            mem_size_Q_q,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreePressX0(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX0.q27[0] ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX0.k      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX0.kN     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX0.Vx     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX0.Vy     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX0.Vz     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX0.RhoBC  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX0.deltaVz));
}
//PressX1 = X-outflow
void CudaMemoryManager::cudaAllocPressX1(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QpressX1.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QpressX1.numberOfBCnodes;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX1.q27[0]), parameter->getD3Qxx()*mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX1.k),                 mem_size_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX1.kN),                mem_size_Q_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX1.Vx),                mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX1.Vy),                mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX1.Vz),                mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX1.RhoBC),             mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QpressX1.deltaVz),           mem_size_Q_q ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX1.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX1.k),                      mem_size_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX1.kN),                     mem_size_Q_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX1.Vx),                     mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX1.Vy),                     mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX1.Vz),                     mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX1.RhoBC),                  mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QpressX1.deltaVz),                mem_size_Q_q ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 2. * (double)mem_size_Q_k + 5. * (double)mem_size_Q_q + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyPressX1(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QpressX1.numberOfBCnodes;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QpressX1.numberOfBCnodes;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX1.q27[0],  parameter->getParH(lev)->QpressX1.q27[0], parameter->getD3Qxx()* mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX1.k,       parameter->getParH(lev)->QpressX1.k,                  mem_size_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX1.kN,      parameter->getParH(lev)->QpressX1.kN,                 mem_size_Q_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX1.Vx,      parameter->getParH(lev)->QpressX1.Vx,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX1.Vy,      parameter->getParH(lev)->QpressX1.Vy,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX1.Vz,      parameter->getParH(lev)->QpressX1.Vz,                 mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX1.RhoBC,   parameter->getParH(lev)->QpressX1.RhoBC,              mem_size_Q_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QpressX1.deltaVz, parameter->getParH(lev)->QpressX1.deltaVz,            mem_size_Q_q,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreePressX1(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX1.q27[0] ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX1.k      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX1.kN     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX1.Vx     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX1.Vy     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX1.Vz     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX1.RhoBC  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QpressX1.deltaVz));
}
void CudaMemoryManager::cudaAllocMeasurePointsIndex(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->kMP),                     parameter->getParH(lev)->memSizeIntkMP  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->VxMP),                    parameter->getParH(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->VyMP),                    parameter->getParH(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->VzMP),                    parameter->getParH(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->RhoMP),                   parameter->getParH(lev)->memSizerealkMP ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->kMP),                         parameter->getParD(lev)->memSizeIntkMP  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->VxMP),                        parameter->getParD(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->VyMP),                        parameter->getParD(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->VzMP),                        parameter->getParD(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->RhoMP),                       parameter->getParD(lev)->memSizerealkMP ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->memSizeIntkMP + 4. * (double)parameter->getParH(lev)->memSizerealkMP;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyMeasurePointsIndex(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->kMP,           parameter->getParH(lev)->kMP,                parameter->getParH(lev)->memSizeIntkMP,      cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->VxMP,          parameter->getParH(lev)->VxMP,               parameter->getParH(lev)->memSizerealkMP,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->VyMP,          parameter->getParH(lev)->VyMP,               parameter->getParH(lev)->memSizerealkMP,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->VzMP,          parameter->getParH(lev)->VzMP,               parameter->getParH(lev)->memSizerealkMP,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->RhoMP,         parameter->getParH(lev)->RhoMP,              parameter->getParH(lev)->memSizerealkMP,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyMeasurePointsToHost(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->kMP,           parameter->getParD(lev)->kMP,                parameter->getParH(lev)->memSizeIntkMP,      cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->VxMP,          parameter->getParD(lev)->VxMP,               parameter->getParH(lev)->memSizerealkMP,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->VyMP,          parameter->getParD(lev)->VyMP,               parameter->getParH(lev)->memSizerealkMP,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->VzMP,          parameter->getParD(lev)->VzMP,               parameter->getParH(lev)->memSizerealkMP,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->RhoMP,         parameter->getParD(lev)->RhoMP,              parameter->getParH(lev)->memSizerealkMP,  cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeMeasurePointsIndex(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->kMP));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->VxMP));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->VyMP));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->VzMP));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->RhoMP));
}
void CudaMemoryManager::cudaAllocFsForCheckPointAndRestart(int lev) const
{
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->distributions.f[0] ),           (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParH(lev)->memSizeRealLBnodes));
}
void CudaMemoryManager::cudaAllocFsForAllLevelsOnHost() const
{
    for (int level = 0; level <= parameter->getMaxLevel(); level++) {
        cudaAllocFsForCheckPointAndRestart(level);
    }
}
void CudaMemoryManager::cudaCopyFsForRestart(int lev) const
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->distributions.f[0],  parameter->getParH(lev)->distributions.f[0],     (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyFsForCheckPoint(int lev) const
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->distributions.f[0],  parameter->getParD(lev)->distributions.f[0],     (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaCopyFsForAllLevelsToHost() const
{
    for (int level = 0; level <= parameter->getMaxLevel(); level++)
        cudaCopyFsForCheckPoint(level);
}
void CudaMemoryManager::cudaFreeFsForCheckPointAndRestart(int lev) const
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->distributions.f[0]));
}
//DragLift
void CudaMemoryManager::cudaAllocDragLift(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(double)*numofelem;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPreX),  mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPreY),  mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPreZ),  mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPostX), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPostY), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPostZ), mem_size  ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->DragPreX),  mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->DragPreY),  mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->DragPreZ),  mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->DragPostX), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->DragPostY), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->DragPostZ), mem_size  ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 6. * (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyDragLift(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(double)*numofelem;

    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->DragPreX, parameter->getParD(lev)->DragPreX, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->DragPreY, parameter->getParD(lev)->DragPreY, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->DragPreZ, parameter->getParD(lev)->DragPreZ, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->DragPostX, parameter->getParD(lev)->DragPostX, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->DragPostY, parameter->getParD(lev)->DragPostY, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->DragPostZ, parameter->getParD(lev)->DragPostZ, mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeDragLift(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->DragPreX));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->DragPreY));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->DragPreZ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->DragPostX));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->DragPostY));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->DragPostZ));
}
//2ndMoments
void CudaMemoryManager::cudaAlloc2ndMoments(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->kxyFromfcNEQ   ), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->kyzFromfcNEQ   ), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->kxzFromfcNEQ   ), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->kxxMyyFromfcNEQ), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->kxxMzzFromfcNEQ), mem_size  ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->kxyFromfcNEQ   ), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->kyzFromfcNEQ   ), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->kxzFromfcNEQ   ), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->kxxMyyFromfcNEQ), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->kxxMzzFromfcNEQ), mem_size  ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 5. * (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopy2ndMoments(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;

    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->kxyFromfcNEQ   , parameter->getParD(lev)->kxyFromfcNEQ   , mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->kyzFromfcNEQ   , parameter->getParD(lev)->kyzFromfcNEQ   , mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->kxzFromfcNEQ   , parameter->getParD(lev)->kxzFromfcNEQ   , mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->kxxMyyFromfcNEQ, parameter->getParD(lev)->kxxMyyFromfcNEQ, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->kxxMzzFromfcNEQ, parameter->getParD(lev)->kxxMzzFromfcNEQ, mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFree2ndMoments(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->kxyFromfcNEQ   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->kyzFromfcNEQ   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->kxzFromfcNEQ   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->kxxMyyFromfcNEQ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->kxxMzzFromfcNEQ));
}
//3rdMoments
void CudaMemoryManager::cudaAlloc3rdMoments(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMbbb ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMabc ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMbac ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMbca ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMcba ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMacb ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMcab ), mem_size ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMbbb ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMabc ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMbac ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMbca ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMcba ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMacb ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMcab ), mem_size ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 7. * (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopy3rdMoments(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;

    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMbbb, parameter->getParD(lev)->CUMbbb, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMabc, parameter->getParD(lev)->CUMabc, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMbac, parameter->getParD(lev)->CUMbac, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMbca, parameter->getParD(lev)->CUMbca, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMcba, parameter->getParD(lev)->CUMcba, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMacb, parameter->getParD(lev)->CUMacb, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMcab, parameter->getParD(lev)->CUMcab, mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFree3rdMoments(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMbbb ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMabc ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMbac ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMbca ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMcba ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMacb ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMcab ));
}
//higher order moments
void CudaMemoryManager::cudaAllocHigherMoments(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMcbb ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMbcb ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMbbc ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMcca ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMcac ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMacc ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMbcc ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMcbc ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMccb ), mem_size ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->CUMccc ), mem_size ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMcbb ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMbcb ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMbbc ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMcca ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMcac ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMacc ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMbcc ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMcbc ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMccb ), mem_size ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->CUMccc ), mem_size ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 10. * (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyHigherMoments(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;

    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMcbb, parameter->getParD(lev)->CUMcbb, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMbcb, parameter->getParD(lev)->CUMbcb, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMbbc, parameter->getParD(lev)->CUMbbc, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMcca, parameter->getParD(lev)->CUMcca, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMcac, parameter->getParD(lev)->CUMcac, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMacc, parameter->getParD(lev)->CUMacc, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMbcc, parameter->getParD(lev)->CUMbcc, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMcbc, parameter->getParD(lev)->CUMcbc, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMccb, parameter->getParD(lev)->CUMccb, mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->CUMccc, parameter->getParD(lev)->CUMccc, mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeHigherMoments(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMcbb ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMbcb ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMbbc ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMcca ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMcac ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMacc ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMbcc ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMcbc ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMccb ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->CUMccc ));
}
//Velcities to fit the Forcing
void CudaMemoryManager::cudaAllocForceVelo(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->VxForce   ), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->VyForce   ), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->VzForce   ), mem_size  ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->VxForce   ), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->VyForce   ), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->VzForce   ), mem_size  ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyForceVelo(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;

    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->VxForce   , parameter->getParD(lev)->VxForce   , mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->VyForce   , parameter->getParD(lev)->VyForce   , mem_size, cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->VzForce   , parameter->getParD(lev)->VzForce   , mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeForceVelo(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->VxForce   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->VyForce   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->VzForce   ));
}
//cp Top
void CudaMemoryManager::cudaAllocCpTop(int lev)
{
    unsigned int mem_size_double = sizeof(double)       * parameter->getParH(lev)->numberOfPointsCpTop;
    unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsCpTop;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->cpPressTop), mem_size_double  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->cpTopIndex), mem_size_int     ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->cpPressTop), mem_size_double      ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->cpTopIndex), mem_size_int         ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_double + (double)mem_size_int;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyCpTopInit(int lev)
{
    unsigned int mem_size_double = sizeof(double)       * parameter->getParH(lev)->numberOfPointsCpTop;
    unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsCpTop;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->cpPressTop, parameter->getParH(lev)->cpPressTop, mem_size_double, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->cpTopIndex, parameter->getParH(lev)->cpTopIndex, mem_size_int,    cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyCpTop(int lev)
{
    unsigned int mem_size_double = sizeof(double)       * parameter->getParH(lev)->numberOfPointsCpTop;
    //unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsCpTop;

    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->cpPressTop, parameter->getParD(lev)->cpPressTop, mem_size_double, cudaMemcpyDeviceToHost));
    //checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->cpTopIndex, parameter->getParD(lev)->cpTopIndex, mem_size_int,    cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeCpTop(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->cpPressTop));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->cpTopIndex));
}
//cp Bottom
void CudaMemoryManager::cudaAllocCpBottom(int lev)
{
    unsigned int mem_size_double = sizeof(double)       * parameter->getParH(lev)->numberOfPointsCpBottom;
    unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsCpBottom;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->cpPressBottom), mem_size_double  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->cpBottomIndex), mem_size_int     ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->cpPressBottom), mem_size_double      ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->cpBottomIndex), mem_size_int         ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_double + (double)mem_size_int;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyCpBottomInit(int lev)
{
    unsigned int mem_size_double = sizeof(double)       * parameter->getParH(lev)->numberOfPointsCpBottom;
    unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsCpBottom;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->cpPressBottom, parameter->getParH(lev)->cpPressBottom, mem_size_double, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->cpBottomIndex, parameter->getParH(lev)->cpBottomIndex, mem_size_int,    cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyCpBottom(int lev)
{
    unsigned int mem_size_double = sizeof(double)       * parameter->getParH(lev)->numberOfPointsCpBottom;
    //unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsCpBottom;

    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->cpPressBottom, parameter->getParD(lev)->cpPressBottom, mem_size_double, cudaMemcpyDeviceToHost));
    //checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->cpBottomIndex, parameter->getParD(lev)->cpBottomIndex, mem_size_int,    cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeCpBottom(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->cpPressBottom));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->cpBottomIndex));
}
//cp Bottom 2
void CudaMemoryManager::cudaAllocCpBottom2(int lev)
{
    unsigned int mem_size_double = sizeof(double)       * parameter->getParH(lev)->numberOfPointsCpBottom2;
    unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsCpBottom2;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->cpPressBottom2), mem_size_double  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->cpBottom2Index), mem_size_int     ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->cpPressBottom2), mem_size_double      ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->cpBottom2Index), mem_size_int         ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_double + (double)mem_size_int;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyCpBottom2Init(int lev)
{
    unsigned int mem_size_double = sizeof(double)       * parameter->getParH(lev)->numberOfPointsCpBottom2;
    unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsCpBottom2;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->cpPressBottom2, parameter->getParH(lev)->cpPressBottom2, mem_size_double, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->cpBottom2Index, parameter->getParH(lev)->cpBottom2Index, mem_size_int,    cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyCpBottom2(int lev)
{
    unsigned int mem_size_double = sizeof(double)       * parameter->getParH(lev)->numberOfPointsCpBottom2;

    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->cpPressBottom2, parameter->getParD(lev)->cpPressBottom2, mem_size_double, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeCpBottom2(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->cpPressBottom2));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->cpBottom2Index));
}
//////////////////////////////////////////////////////////////////////////
//advection diffusion
void CudaMemoryManager::cudaAllocConcentration(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->concentration), parameter->getParH(lev)->memSizeRealLBnodes));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->concentration), parameter->getParD(lev)->memSizeRealLBnodes));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->memSizeRealLBnodes;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyConcentrationDeviceToHost(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->concentration, parameter->getParD(lev)->concentration,  parameter->getParH(lev)->memSizeRealLBnodes , cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaCopyConcentrationHostToDevice(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->concentration, parameter->getParH(lev)->concentration, parameter->getParH(lev)->memSizeRealLBnodes, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeConcentration(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->concentration));
}
//////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocTempFs(int lev)
{
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->distributionsAD.f[0]), 27*parameter->getParH(lev)->memSizeRealLBnodes));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)(27 * parameter->getParH(lev)->memSizeRealLBnodes);
    setMemsizeGPU(tmp, false);
}
//////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocTempPressBC(int lev)
{
    unsigned int mem_size_TempPress_k = sizeof(int)*parameter->getParH(lev)->TempPress.kTemp;
    unsigned int mem_size_TempPress_q = sizeof(real)*parameter->getParH(lev)->TempPress.kTemp;

    // Host Memory
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->TempPress.temp, mem_size_TempPress_q ));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->TempPress.velo, mem_size_TempPress_q ));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->TempPress.k,    mem_size_TempPress_k ));

    // Device Memory
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->TempPress.temp, mem_size_TempPress_q));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->TempPress.velo, mem_size_TempPress_q));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->TempPress.k,    mem_size_TempPress_k));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)(2.0 * mem_size_TempPress_q) + (double)(mem_size_TempPress_k);
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyTempPressBCHD(int lev)
{
    unsigned int mem_size_TempPress_k = sizeof(int)*parameter->getParH(lev)->TempPress.kTemp;
    unsigned int mem_size_TempPress_q = sizeof(real)*parameter->getParH(lev)->TempPress.kTemp;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->TempPress.temp, parameter->getParH(lev)->TempPress.temp, mem_size_TempPress_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->TempPress.velo, parameter->getParH(lev)->TempPress.velo, mem_size_TempPress_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->TempPress.k,    parameter->getParH(lev)->TempPress.k,    mem_size_TempPress_k,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeTempPressBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->TempPress.temp));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->TempPress.velo));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->TempPress.k   ));
}
//////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocTempVeloBC(int lev)
{
    unsigned int mem_size_TempVel_k = sizeof(int)*parameter->getParH(lev)->TempVel.kTemp;
    unsigned int mem_size_TempVel_q = sizeof(real)*parameter->getParH(lev)->TempVel.kTemp;

    printf("mem_size_TempVel_k = %d,  mem_size_TempVel_q = %d \n", mem_size_TempVel_k, mem_size_TempVel_q);
    // Host Memory
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->TempVel.temp,      mem_size_TempVel_q ));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->TempVel.tempPulse, mem_size_TempVel_q ));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->TempVel.velo,      mem_size_TempVel_q ));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->TempVel.k,         mem_size_TempVel_k ));

    // Device Memory
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->TempVel.temp,      mem_size_TempVel_q));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->TempVel.tempPulse, mem_size_TempVel_q));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->TempVel.velo,      mem_size_TempVel_q));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->TempVel.k,         mem_size_TempVel_k));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)(3.0 * mem_size_TempVel_q) + (double)(mem_size_TempVel_k);
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyTempVeloBCHD(int lev)
{
    unsigned int mem_size_TempVel_k = sizeof(int)*parameter->getParH(lev)->TempVel.kTemp;
    unsigned int mem_size_TempVel_q = sizeof(real)*parameter->getParH(lev)->TempVel.kTemp;

    printf("mem_size_TempVel_k = %d,  mem_size_TempVel_q = %d \n", mem_size_TempVel_k, mem_size_TempVel_q);
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->TempVel.temp,      parameter->getParH(lev)->TempVel.temp,      mem_size_TempVel_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->TempVel.tempPulse, parameter->getParH(lev)->TempVel.tempPulse, mem_size_TempVel_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->TempVel.velo,      parameter->getParH(lev)->TempVel.velo,      mem_size_TempVel_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->TempVel.k,         parameter->getParH(lev)->TempVel.k,         mem_size_TempVel_k,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeTempVeloBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->TempVel.temp     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->TempVel.tempPulse));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->TempVel.velo     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->TempVel.k        ));
}
//////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocTempNoSlipBC(int lev)
{
    unsigned int mem_size_Temp_k = sizeof(int)*parameter->getParH(lev)->Temp.kTemp;
    unsigned int mem_size_Temp_q = sizeof(real)*parameter->getParH(lev)->Temp.kTemp;

    // Host Memory
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->Temp.temp, mem_size_Temp_q ));
    checkCudaErrors( cudaMallocHost((void**) &parameter->getParH(lev)->Temp.k,    mem_size_Temp_k ));

    // Device Memory
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->Temp.temp, mem_size_Temp_q));
    checkCudaErrors( cudaMalloc((void**) &parameter->getParD(lev)->Temp.k,    mem_size_Temp_k));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)(mem_size_Temp_q) + (double)(mem_size_Temp_k);
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyTempNoSlipBCHD(int lev)
{
    unsigned int mem_size_Temp_k = sizeof(int)*parameter->getParH(lev)->Temp.kTemp;
    unsigned int mem_size_Temp_q = sizeof(real)*parameter->getParH(lev)->Temp.kTemp;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Temp.temp, parameter->getParH(lev)->Temp.temp, mem_size_Temp_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Temp.k,    parameter->getParH(lev)->Temp.k,    mem_size_Temp_k,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeTempNoSlipBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Temp.temp));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Temp.k   ));
}
//PlaneConc
void CudaMemoryManager::cudaAllocPlaneConcIn(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->ConcPlaneIn), mem_size  ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->ConcPlaneIn), mem_size  ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyPlaneConcIn(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->ConcPlaneIn,   parameter->getParD(lev)->ConcPlaneIn,   mem_size, cudaMemcpyDeviceToHost));
}
//////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocPlaneConcOut1(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->ConcPlaneOut1), mem_size  ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->ConcPlaneOut1), mem_size  ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyPlaneConcOut1(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->ConcPlaneOut1, parameter->getParD(lev)->ConcPlaneOut1, mem_size, cudaMemcpyDeviceToHost));
}
//////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocPlaneConcOut2(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->ConcPlaneOut2), mem_size  ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->ConcPlaneOut2), mem_size  ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyPlaneConcOut2(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(real)*numofelem;
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->ConcPlaneOut2, parameter->getParD(lev)->ConcPlaneOut2, mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreePlaneConc(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->ConcPlaneIn));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->ConcPlaneOut1));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->ConcPlaneOut2));
}
//////////////////////////////////////////////////////////////////////////
//concentration file
void CudaMemoryManager::cudaAllocConcFile(int lev)
{
    unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsConc;

    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->concIndex), mem_size_int     ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->concIndex), mem_size_int         ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_int;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyConcFile(int lev)
{
    unsigned int mem_size_int    = sizeof(unsigned int) * parameter->getParH(lev)->numberOfPointsConc;

    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->concIndex, parameter->getParH(lev)->concIndex, mem_size_int,    cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeConcFile(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->concIndex));
}
//////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocMeanOutAD(int lev)
{
    //Host
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->rho_SP_Med_Out),   parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->vx_SP_Med_Out),    parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->vy_SP_Med_Out),    parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->vz_SP_Med_Out),    parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->press_SP_Med_Out), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->Conc_Med_Out),     parameter->getParH(lev)->memSizeRealLBnodes));
}
void CudaMemoryManager::cudaFreeMeanOutAD(int lev)
{
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->vx_SP_Med_Out));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->vy_SP_Med_Out));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->vz_SP_Med_Out));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->rho_SP_Med_Out));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->press_SP_Med_Out));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->Conc_Med_Out));
}


//////////////////////////////////////////////////////////////////////////
//Process Neighbors
//1D domain decomposition
void CudaMemoryManager::cudaAllocProcessNeighbor(int lev, unsigned int processNeighbor)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].index ),                             parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].index ),                             parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].index ),                                 parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].index ),                                 parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp =
        (double)parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].memsizeIndex +
        (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].memsizeFs +
        (double)parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].memsizeIndex +
        (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].memsizeFs;
    setMemsizeGPU(tmp, false);
    //printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void CudaMemoryManager::cudaCopyProcessNeighborIndex(int lev, unsigned int processNeighbor)
{
    //copy send Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
    //copy recv Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborFsHD(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].f[0],
                                parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].f[0],
                                parameter->getD3Qxx() * parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].memsizeFs,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborFsDH(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].f[0],
                                parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].f[0],
                                parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].memsizeFs,
                                cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeProcessNeighbor(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].index ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].f[0]     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].index  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].f[0]     ));
}

////////////////////////////////////////////////////////////////////////////////////
//  3D domain decomposition convection diffusion
//  X  /////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocProcessNeighborADX(int lev, unsigned int processNeighbor)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].index ),      parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].f[0]  ), 27 * parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].index ),      parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].f[0]  ), 27 * parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].index ),      parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].f[0]  ), 27 * parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].index ), parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].f[0]  ), 27 * parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp =
        (double)parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeIndex +
        (double)27.0*(double)parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeFs +
        (double)parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeIndex +
        (double)27.0*(double)parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeFs;
    setMemsizeGPU(tmp, false);
    //printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void CudaMemoryManager::cudaCopyProcessNeighborADXIndex(int lev, unsigned int processNeighbor)
{
    //copy send Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
    //copy recv Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborADXFsHD(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].f[0],
                                parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].f[0],
                                27 * parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].memsizeFs,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborADXFsDH(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].f[0],
                                parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].f[0],
                                27 * parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].memsizeFs,
                                cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeProcessNeighborADX(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].index ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].f[0]     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].index  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].f[0]     ));
}
//  Y  /////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocProcessNeighborADY(int lev, unsigned int processNeighbor)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].index ),      parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].f[0]  ), 27 * parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].index ),      parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].f[0]  ), 27 * parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].index ),      parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].f[0]  ), 27 * parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].index ),      parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].f[0]  ), 27 * parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp =
        (double)parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex +
        (double)27.0*(double)parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs +
        (double)parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex +
        (double)27.0*(double)parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs;
    setMemsizeGPU(tmp, false);
    //printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void CudaMemoryManager::cudaCopyProcessNeighborADYIndex(int lev, unsigned int processNeighbor)
{
    //copy send Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
    //copy recv Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborADYFsHD(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].f[0],
                                parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].f[0],
                                27 * parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborADYFsDH(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].f[0],
                                parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].f[0],
                                27 * parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs,
                                cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeProcessNeighborADY(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].index ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].f[0]     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].index  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].f[0]     ));
}
//  Z  /////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocProcessNeighborADZ(int lev, unsigned int processNeighbor)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].index ),      parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].f[0]  ), 27 * parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].index ),      parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].f[0]  ), 27 * parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].index ),      parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].f[0]  ), 27 * parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].index ),      parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].f[0]  ), 27 * parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp =
        (double)parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex +
        (double)27.0*(double)parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs +
        (double)parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex +
        (double)27.0*(double)parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs;
    setMemsizeGPU(tmp, false);
    //printf("memsize GPU for neighbors %f \n",tmp/1000000.0);
}
void CudaMemoryManager::cudaCopyProcessNeighborADZIndex(int lev, unsigned int processNeighbor)
{
    //copy send Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].index,
                                parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
    //copy recv Index
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].index,
                                parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborADZFsHD(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].f[0],
                                parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].f[0],
                                27 * parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborADZFsDH(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].f[0],
                                parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].f[0],
                                27 * parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs,
                                cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeProcessNeighborADZ(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].index ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].f[0]     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].index  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].f[0]     ));
}

void CudaMemoryManager::cudaAlloc2ndOrderDerivitivesIsoTest(int lev)
{
    //Host
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->dxxUx), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->dyyUy), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->dzzUz), parameter->getParH(lev)->memSizeRealLBnodes));
    //Device (spinning ship)
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->dxxUx), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->dyyUy), parameter->getParH(lev)->memSizeRealLBnodes));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->dzzUz), parameter->getParH(lev)->memSizeRealLBnodes));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)parameter->getParH(lev)->memSizeRealLBnodes;
    setMemsizeGPU(tmp, false);
    //printf("Coord = %f MB",tmp/1000000.);
}
void CudaMemoryManager::cudaCopy2ndOrderDerivitivesIsoTestDH(int lev)
{
    //copy device to host
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->dxxUx, parameter->getParD(lev)->dxxUx, parameter->getParH(lev)->memSizeRealLBnodes, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->dyyUy, parameter->getParD(lev)->dyyUy, parameter->getParH(lev)->memSizeRealLBnodes, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->dzzUz, parameter->getParD(lev)->dzzUz, parameter->getParH(lev)->memSizeRealLBnodes, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaCopy2ndOrderDerivitivesIsoTestHD(int lev)
{
    //copy host to device
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->dxxUx, parameter->getParH(lev)->dxxUx, parameter->getParH(lev)->memSizeRealLBnodes, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->dyyUy, parameter->getParH(lev)->dyyUy, parameter->getParH(lev)->memSizeRealLBnodes, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->dzzUz, parameter->getParH(lev)->dzzUz, parameter->getParH(lev)->memSizeRealLBnodes, cudaMemcpyHostToDevice));

}
void CudaMemoryManager::cudaFree2ndOrderDerivitivesIsoTest(int lev)
{
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->dxxUx));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->dyyUy));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->dzzUz));

}

void CudaMemoryManager::cudaAllocTaggedFluidNodeIndices(CollisionTemplate tag, int lev) {
    uint mem_size_tagged_fluid_nodes = sizeof(uint) * parameter->getParH(lev)->numberOfTaggedFluidNodes[tag];
    // Host
    checkCudaErrors(cudaMallocHost((void **)&(parameter->getParH(lev)->taggedFluidNodeIndices[tag]), mem_size_tagged_fluid_nodes));
    // Device
    checkCudaErrors(cudaMalloc((void **)&(parameter->getParD(lev)->taggedFluidNodeIndices[tag]), mem_size_tagged_fluid_nodes));
    //////////////////////////////////////////////////////////////////////////
    setMemsizeGPU((double)mem_size_tagged_fluid_nodes, false);
}

void CudaMemoryManager::cudaCopyTaggedFluidNodeIndices(CollisionTemplate tag, int lev) {
    uint mem_size_tagged_fluid_nodes = sizeof(uint) * parameter->getParH(lev)->numberOfTaggedFluidNodes[tag];
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->taggedFluidNodeIndices[tag],
                               parameter->getParH(lev)->taggedFluidNodeIndices[tag],
                               mem_size_tagged_fluid_nodes, cudaMemcpyHostToDevice));
}

void CudaMemoryManager::cudaFreeTaggedFluidNodeIndices(CollisionTemplate tag, int lev) {
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->taggedFluidNodeIndices[tag]));
}

////////////////////////////////////////////////////////////////////////////////////
//  ActuatorFarm
///////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocBladeGeometries(ActuatorFarm* actuatorFarm)
{
    const uint sizeRealTurbine = sizeof(real)*actuatorFarm->getNumberOfTurbines();
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->turbinePosXH, sizeRealTurbine) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->turbinePosYH, sizeRealTurbine) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->turbinePosZH, sizeRealTurbine) );

    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->turbinePosXD, sizeRealTurbine) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->turbinePosYD, sizeRealTurbine) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->turbinePosZD, sizeRealTurbine) );
    setMemsizeGPU(sizeof(real)*4.f*actuatorFarm->getNumberOfTurbines(), false);

}
void CudaMemoryManager::cudaCopyBladeGeometriesHtoD(ActuatorFarm* actuatorFarm)
{
    const uint sizeRealTurbine = sizeof(real)*actuatorFarm->getNumberOfTurbines();
    checkCudaErrors( cudaMemcpy(actuatorFarm->turbinePosXD, actuatorFarm->turbinePosXH, sizeRealTurbine, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->turbinePosYD, actuatorFarm->turbinePosYH, sizeRealTurbine, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->turbinePosZD, actuatorFarm->turbinePosZH, sizeRealTurbine, cudaMemcpyHostToDevice) );

}
void CudaMemoryManager::cudaCopyBladeGeometriesDtoH(ActuatorFarm* actuatorFarm)
{
    uint sizeRealTurbine = sizeof(real)*actuatorFarm->getNumberOfTurbines();
    checkCudaErrors( cudaMemcpy(actuatorFarm->turbinePosXH, actuatorFarm->turbinePosXD, sizeRealTurbine, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->turbinePosYH, actuatorFarm->turbinePosYD, sizeRealTurbine, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->turbinePosZH, actuatorFarm->turbinePosZD, sizeRealTurbine, cudaMemcpyDeviceToHost) );

}
void CudaMemoryManager::cudaFreeBladeGeometries(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaFree(actuatorFarm->turbinePosXD) );
    checkCudaErrors( cudaFree(actuatorFarm->turbinePosYD) );
    checkCudaErrors( cudaFree(actuatorFarm->turbinePosZD) );    
    
    checkCudaErrors( cudaFreeHost(actuatorFarm->turbinePosXH) );
    checkCudaErrors( cudaFreeHost(actuatorFarm->turbinePosYH) );
    checkCudaErrors( cudaFreeHost(actuatorFarm->turbinePosZH) );
}

void CudaMemoryManager::cudaAllocBladeCoords(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeCoordsXH, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeCoordsYH, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeCoordsZH, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeCoordsXDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeCoordsYDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeCoordsZDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );    
    
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeCoordsXDPreviousTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeCoordsYDPreviousTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeCoordsZDPreviousTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );

    setMemsizeGPU(6.f*actuatorFarm->getNumberOfGridNodes(), false);
}

void CudaMemoryManager::cudaCopyBladeCoordsHtoD(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeCoordsXDCurrentTimestep, actuatorFarm->bladeCoordsXH, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeCoordsYDCurrentTimestep, actuatorFarm->bladeCoordsYH, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeCoordsZDCurrentTimestep, actuatorFarm->bladeCoordsZH, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyBladeCoordsDtoH(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeCoordsXH, actuatorFarm->bladeCoordsXDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeCoordsYH, actuatorFarm->bladeCoordsYDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeCoordsZH, actuatorFarm->bladeCoordsZDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeBladeCoords(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaFree(actuatorFarm->bladeCoordsXDCurrentTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeCoordsYDCurrentTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeCoordsZDCurrentTimestep) );

    checkCudaErrors( cudaFree(actuatorFarm->bladeCoordsXDPreviousTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeCoordsYDPreviousTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeCoordsZDPreviousTimestep) );

    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeCoordsXH) );
    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeCoordsYH) );
    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeCoordsZH) );
}

void CudaMemoryManager::cudaAllocBladeIndices(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeIndicesH, sizeof(uint)*actuatorFarm->getNumberOfGridNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeIndicesD, sizeof(uint)*actuatorFarm->getNumberOfGridNodes()) );

    setMemsizeGPU(sizeof(uint)*actuatorFarm->getNumberOfGridNodes(), false);
}

void CudaMemoryManager::cudaCopyBladeIndicesHtoD(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeIndicesD, actuatorFarm->bladeIndicesH, sizeof(uint)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaFreeBladeIndices(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaFree(actuatorFarm->bladeIndicesD) );

    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeIndicesH) );
}

void CudaMemoryManager::cudaAllocBladeVelocities(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeVelocitiesXH, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeVelocitiesYH, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeVelocitiesZH, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeVelocitiesXDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeVelocitiesYDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeVelocitiesZDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeVelocitiesXDPreviousTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeVelocitiesYDPreviousTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeVelocitiesZDPreviousTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );

    setMemsizeGPU(3.*sizeof(real)*actuatorFarm->getNumberOfGridNodes(), false);
}

void CudaMemoryManager::cudaCopyBladeVelocitiesHtoD(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeVelocitiesXDCurrentTimestep, actuatorFarm->bladeVelocitiesXH, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeVelocitiesYDCurrentTimestep, actuatorFarm->bladeVelocitiesYH, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeVelocitiesZDCurrentTimestep, actuatorFarm->bladeVelocitiesZH, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyBladeVelocitiesDtoH(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeVelocitiesXH, actuatorFarm->bladeVelocitiesXDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeVelocitiesYH, actuatorFarm->bladeVelocitiesYDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeVelocitiesZH, actuatorFarm->bladeVelocitiesZDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeBladeVelocities(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaFree(actuatorFarm->bladeVelocitiesXDCurrentTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeVelocitiesYDCurrentTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeVelocitiesZDCurrentTimestep) );    
    
    checkCudaErrors( cudaFree(actuatorFarm->bladeVelocitiesXDPreviousTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeVelocitiesYDPreviousTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeVelocitiesZDPreviousTimestep) );

    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeVelocitiesXH) );
    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeVelocitiesYH) );
    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeVelocitiesZH) );
}

void CudaMemoryManager::cudaAllocBladeForces(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeForcesXH, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeForcesYH, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorFarm->bladeForcesZH, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeForcesXDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeForcesYDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeForcesZDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeForcesXDPreviousTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeForcesYDPreviousTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorFarm->bladeForcesZDPreviousTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes()) );

    setMemsizeGPU(3.*sizeof(real)*actuatorFarm->getNumberOfGridNodes(), false);
}

void CudaMemoryManager::cudaCopyBladeForcesHtoD(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeForcesXDCurrentTimestep, actuatorFarm->bladeForcesXH, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeForcesYDCurrentTimestep, actuatorFarm->bladeForcesYH, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeForcesZDCurrentTimestep, actuatorFarm->bladeForcesZH, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyBladeForcesDtoH(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeForcesXH, actuatorFarm->bladeForcesXDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeForcesYH, actuatorFarm->bladeForcesYDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorFarm->bladeForcesZH, actuatorFarm->bladeForcesZDCurrentTimestep, sizeof(real)*actuatorFarm->getNumberOfGridNodes(), cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeBladeForces(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaFree(actuatorFarm->bladeForcesXDCurrentTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeForcesYDCurrentTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeForcesZDCurrentTimestep) );

    checkCudaErrors( cudaFree(actuatorFarm->bladeForcesXDPreviousTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeForcesYDPreviousTimestep) );
    checkCudaErrors( cudaFree(actuatorFarm->bladeForcesZDPreviousTimestep) );

    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeForcesXH) );
    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeForcesYH) );
    checkCudaErrors( cudaFreeHost(actuatorFarm->bladeForcesZH) );
}

void CudaMemoryManager::cudaAllocSphereIndices(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMallocHost((void**) &(actuatorFarm->boundingSphereIndicesH), sizeof(int)*actuatorFarm->getNumberOfIndices()));
    checkCudaErrors( cudaMalloc((void**) &(actuatorFarm->boundingSphereIndicesD), sizeof(int)*actuatorFarm->getNumberOfIndices()));
    setMemsizeGPU(sizeof(int)*actuatorFarm->getNumberOfIndices(), false);
}

void CudaMemoryManager::cudaCopySphereIndicesHtoD(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaMemcpy(actuatorFarm->boundingSphereIndicesD, actuatorFarm->boundingSphereIndicesH, sizeof(int)*actuatorFarm->getNumberOfIndices(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaFreeSphereIndices(ActuatorFarm* actuatorFarm)
{
    checkCudaErrors( cudaFreeHost(actuatorFarm->boundingSphereIndicesH) );
    checkCudaErrors( cudaFree(actuatorFarm->boundingSphereIndicesD) );
}

////////////////////////////////////////////////////////////////////////////////////
//  Probe
///////////////////////////////////////////////////////////////////////////////

void CudaMemoryManager::cudaAllocProbeDistances(Probe* probe, int level)
{
    size_t tmp = sizeof(real)*probe->getProbeStruct(level)->nPoints;

    checkCudaErrors( cudaMallocHost((void**) &probe->getProbeStruct(level)->distXH, tmp) );
    checkCudaErrors( cudaMallocHost((void**) &probe->getProbeStruct(level)->distYH, tmp) );
    checkCudaErrors( cudaMallocHost((void**) &probe->getProbeStruct(level)->distZH, tmp) );

    checkCudaErrors( cudaMalloc    ((void**) &probe->getProbeStruct(level)->distXD, tmp) );
    checkCudaErrors( cudaMalloc    ((void**) &probe->getProbeStruct(level)->distYD, tmp) );
    checkCudaErrors( cudaMalloc    ((void**) &probe->getProbeStruct(level)->distZD, tmp) );
    setMemsizeGPU(3.f*tmp, false);
}
void CudaMemoryManager::cudaCopyProbeDistancesHtoD(Probe* probe, int level)
{
    size_t tmp = sizeof(real)*probe->getProbeStruct(level)->nPoints;

    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->distXD, probe->getProbeStruct(level)->distXH, tmp, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->distYD, probe->getProbeStruct(level)->distYH, tmp, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->distZD, probe->getProbeStruct(level)->distZH, tmp, cudaMemcpyHostToDevice) );
}
void CudaMemoryManager::cudaCopyProbeDistancesDtoH(Probe* probe, int level)
{
    size_t tmp = sizeof(real)*probe->getProbeStruct(level)->nPoints;

    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->distXH, probe->getProbeStruct(level)->distXD, tmp, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->distYH, probe->getProbeStruct(level)->distXD, tmp, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->distZH, probe->getProbeStruct(level)->distXD, tmp, cudaMemcpyDeviceToHost) );
}
void CudaMemoryManager::cudaFreeProbeDistances(Probe* probe, int level)
{
    checkCudaErrors( cudaFreeHost(probe->getProbeStruct(level)->distXH) );
    checkCudaErrors( cudaFreeHost(probe->getProbeStruct(level)->distYH) );
    checkCudaErrors( cudaFreeHost(probe->getProbeStruct(level)->distZH) );

    checkCudaErrors( cudaFree    (probe->getProbeStruct(level)->distXD) );
    checkCudaErrors( cudaFree    (probe->getProbeStruct(level)->distYD) );
    checkCudaErrors( cudaFree    (probe->getProbeStruct(level)->distZD) );
}

void CudaMemoryManager::cudaAllocProbeIndices(Probe* probe, int level)
{
    size_t tmp = sizeof(int)*probe->getProbeStruct(level)->nIndices;
    checkCudaErrors( cudaMallocHost((void**) &probe->getProbeStruct(level)->pointIndicesH, tmp) );
    checkCudaErrors( cudaMalloc    ((void**) &probe->getProbeStruct(level)->pointIndicesD, tmp) );
    setMemsizeGPU(1.f*tmp, false);
}

void CudaMemoryManager::cudaCopyProbeIndicesHtoD(Probe* probe, int level)
{
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->pointIndicesD, probe->getProbeStruct(level)->pointIndicesH, sizeof(int)*probe->getProbeStruct(level)->nIndices, cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyProbeIndicesDtoH(Probe* probe, int level)
{
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->pointIndicesH, probe->getProbeStruct(level)->pointIndicesD, sizeof(int)*probe->getProbeStruct(level)->nIndices, cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeProbeIndices(Probe* probe, int level)
{
    checkCudaErrors( cudaFreeHost(probe->getProbeStruct(level)->pointIndicesH) );
    checkCudaErrors( cudaFree    (probe->getProbeStruct(level)->pointIndicesD) );
}

void CudaMemoryManager::cudaAllocProbeQuantityArray(Probe* probe, int level)
{
    auto probeStruct = probe->getProbeStruct(level);
    size_t tmp = sizeof(real)*probeStruct->nArrays*probeStruct->nPoints*probeStruct->nTimesteps;

    checkCudaErrors( cudaMallocHost((void**) &probeStruct->quantitiesArrayH, tmp) );
    if(probe->getHasDeviceQuantityArray())
    {
        checkCudaErrors( cudaMalloc    ((void**) &probeStruct->quantitiesArrayD, tmp) );
        setMemsizeGPU(1.f*tmp, false);
    }
}

void CudaMemoryManager::cudaCopyProbeQuantityArrayHtoD(Probe* probe, int level)
{
    auto probeStruct = probe->getProbeStruct(level);
    size_t tmp = sizeof(real)*probeStruct->nArrays*probeStruct->nPoints*probeStruct->nTimesteps;
    checkCudaErrors( cudaMemcpy(probeStruct->quantitiesArrayD, probeStruct->quantitiesArrayH, tmp, cudaMemcpyHostToDevice) );
}
void CudaMemoryManager::cudaCopyProbeQuantityArrayDtoH(Probe* probe, int level)
{
    auto probeStruct = probe->getProbeStruct(level);
    size_t tmp = sizeof(real)*probeStruct->nArrays*probeStruct->nPoints*probeStruct->nTimesteps;
    checkCudaErrors( cudaMemcpy(probeStruct->quantitiesArrayH, probeStruct->quantitiesArrayD, tmp, cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeProbeQuantityArray(Probe* probe, int level)
{
    checkCudaErrors( cudaFreeHost(probe->getProbeStruct(level)->quantitiesArrayH) );
    if(probe->getHasDeviceQuantityArray())
        checkCudaErrors( cudaFree    (probe->getProbeStruct(level)->quantitiesArrayD) );
}

void CudaMemoryManager::cudaAllocProbeQuantitiesAndOffsets(Probe* probe, int level)
{
    size_t tmpA = int(Statistic::LAST)*sizeof(int);
    size_t tmpQ = int(Statistic::LAST)*sizeof(bool);
    checkCudaErrors( cudaMallocHost((void**) &probe->getProbeStruct(level)->quantitiesH, tmpQ) );
    checkCudaErrors( cudaMalloc    ((void**) &probe->getProbeStruct(level)->quantitiesD, tmpQ) );
    checkCudaErrors( cudaMallocHost((void**) &probe->getProbeStruct(level)->arrayOffsetsH, tmpA) );
    checkCudaErrors( cudaMalloc    ((void**) &probe->getProbeStruct(level)->arrayOffsetsD, tmpA) );
    setMemsizeGPU(tmpA+tmpQ, false);
}

void CudaMemoryManager::cudaCopyProbeQuantitiesAndOffsetsHtoD(Probe* probe, int level)
{
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->quantitiesD, probe->getProbeStruct(level)->quantitiesH, int(Statistic::LAST)*sizeof(bool), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->arrayOffsetsD, probe->getProbeStruct(level)->arrayOffsetsH, int(Statistic::LAST)*sizeof(int), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyProbeQuantitiesAndOffsetsDtoH(Probe* probe, int level)
{
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->quantitiesH, probe->getProbeStruct(level)->quantitiesD, int(Statistic::LAST)*sizeof(bool), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->arrayOffsetsH, probe->getProbeStruct(level)->arrayOffsetsD, int(Statistic::LAST)*sizeof(int), cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeProbeQuantitiesAndOffsets(Probe* probe, int level)
{
    checkCudaErrors( cudaFreeHost(probe->getProbeStruct(level)->quantitiesH) );
    checkCudaErrors( cudaFreeHost(probe->getProbeStruct(level)->arrayOffsetsH) );
    checkCudaErrors( cudaFree    (probe->getProbeStruct(level)->quantitiesD) );
    checkCudaErrors( cudaFree    (probe->getProbeStruct(level)->arrayOffsetsD) );
}

void CudaMemoryManager::cudaAllocPrecursorWriter(PrecursorWriter* writer, int level)
{
    auto prec =  writer->getPrecursorStruct(level);
    size_t indSize = prec->numberOfPointsInBC*sizeof(uint);

    checkCudaErrors( cudaStreamCreate(&prec->stream) );

    checkCudaErrors( cudaMallocHost((void**) &prec->indicesH, indSize));
    checkCudaErrors( cudaMalloc((void**) &prec->indicesD, indSize));

    size_t dataSize  = prec->numberOfPointsInBC*sizeof(real)*prec->numberOfQuantities;
    size_t dataSizeH = dataSize * prec->numberOfTimestepsPerFile;
    
    checkCudaErrors( cudaMallocHost((void**) &prec->dataH, dataSizeH));
    checkCudaErrors( cudaMallocHost((void**) &prec->bufferH, dataSizeH));
    checkCudaErrors( cudaMalloc((void**) &prec->dataD, dataSize));
    checkCudaErrors( cudaMalloc((void**) &prec->bufferD, dataSize));

    setMemsizeGPU(indSize+2*dataSize, false);
}

void CudaMemoryManager::cudaCopyPrecursorWriterIndicesHtoD(PrecursorWriter* writer, int level)
{
    checkCudaErrors( cudaMemcpy(writer->getPrecursorStruct(level)->indicesD, writer->getPrecursorStruct(level)->indicesH, writer->getPrecursorStruct(level)->numberOfPointsInBC*sizeof(uint), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyPrecursorWriterOutputVariablesDtoH(PrecursorWriter* writer, int level)
{
    auto prec =  writer->getPrecursorStruct(level);
    int sizeTimestep = prec->numberOfPointsInBC*prec->numberOfQuantities;

    checkCudaErrors( cudaStreamSynchronize(prec->stream) );
    checkCudaErrors( cudaMemcpyAsync( &prec->bufferH[prec->numberOfTimestepsBuffered*sizeTimestep], prec->bufferD, sizeof(real)*sizeTimestep, cudaMemcpyDeviceToHost, prec->stream));
}

void CudaMemoryManager::cudaFreePrecursorWriter(PrecursorWriter* writer, int level)
{
    checkCudaErrors( cudaFreeHost(writer->getPrecursorStruct(level)->indicesH));
    checkCudaErrors( cudaFree(writer->getPrecursorStruct(level)->indicesD));

    checkCudaErrors( cudaFreeHost(writer->getPrecursorStruct(level)->dataH));
    checkCudaErrors( cudaFreeHost(writer->getPrecursorStruct(level)->bufferH));
    checkCudaErrors( cudaFree(writer->getPrecursorStruct(level)->dataD));
    checkCudaErrors( cudaFree(writer->getPrecursorStruct(level)->bufferD));
}


CudaMemoryManager::CudaMemoryManager(std::shared_ptr<Parameter> parameter) : parameter(parameter)
{

}


void CudaMemoryManager::setMemsizeGPU(double admem, bool reset)
{
    if (reset == true)
    {
        this->memsizeGPU = 0.;
    }
    else
    {
        this->memsizeGPU += admem;
    }
}

double CudaMemoryManager::getMemsizeGPU()
{
    return this->memsizeGPU;
}

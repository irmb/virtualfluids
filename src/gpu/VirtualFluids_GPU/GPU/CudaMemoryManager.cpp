#include "CudaMemoryManager.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Parameter/Parameter.h>
#include "Parameter/CudaStreamManager.h"
#include "PreCollisionInteractor/ActuatorLine.h"
#include "PreCollisionInteractor/Probes/Probe.h"

#include "Calculation/PorousMedia.h"

#include "lbm/constants/NumericConstants.h"

//void CudaMemoryManager::cudaAllocFull(int lev) //DEPRECATED: related to full matrix
//{
//    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geo      ), parameter->getParH(lev)->mem_size_int  ));
//    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->k        ), parameter->getParH(lev)->mem_size_int  ));
//}
//void CudaMemoryManager::cudaFreeFull(int lev) //DEPRECATED: related to full matrix
//{
//    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geo   ));
//    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->k     ));
//}
void CudaMemoryManager::cudaCopyPrint(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityX   , parameter->getParD(lev)->velocityX   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityY   , parameter->getParD(lev)->velocityY   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityZ   , parameter->getParD(lev)->velocityZ   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->rho  , parameter->getParD(lev)->rho  , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->pressure, parameter->getParD(lev)->pressure, parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));

    if(parameter->getIsBodyForce())
    {
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->forceX_SP   , parameter->getParD(lev)->forceX_SP   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->forceY_SP   , parameter->getParD(lev)->forceY_SP   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->forceZ_SP   , parameter->getParD(lev)->forceZ_SP   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    }

    if(parameter->getUseTurbulentViscosity())
    {
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->turbViscosity   , parameter->getParD(lev)->turbViscosity   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    }
}
void CudaMemoryManager::cudaCopyMedianPrint(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vx_SP_Med   , parameter->getParD(lev)->vx_SP_Med   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vy_SP_Med   , parameter->getParD(lev)->vy_SP_Med   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vz_SP_Med   , parameter->getParD(lev)->vz_SP_Med   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->rho_SP_Med  , parameter->getParD(lev)->rho_SP_Med  , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->press_SP_Med, parameter->getParD(lev)->press_SP_Med, parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaAllocCoord(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordinateX      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordinateY      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordinateZ      ), parameter->getParH(lev)->mem_size_real_SP  ));
	//Device (spinning ship + uppsala)
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordinateX      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordinateY      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordinateZ      ), parameter->getParH(lev)->mem_size_real_SP  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (double)parameter->getParH(lev)->mem_size_real_SP;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyCoord(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordinateX,  parameter->getParH(lev)->coordinateX,  parameter->getParH(lev)->mem_size_real_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordinateY,  parameter->getParH(lev)->coordinateY,  parameter->getParH(lev)->mem_size_real_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordinateZ,  parameter->getParH(lev)->coordinateZ,  parameter->getParH(lev)->mem_size_real_SP     , cudaMemcpyHostToDevice));
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
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->forceX_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->forceY_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->forceZ_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->forceX_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->forceY_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->forceZ_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (double)parameter->getParH(lev)->mem_size_real_SP;
	setMemsizeGPU(tmp, false);

}
void CudaMemoryManager::cudaCopyBodyForce(int lev)
{
   	//copy host to device
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->forceX_SP,  parameter->getParH(lev)->forceX_SP,  parameter->getParH(lev)->mem_size_real_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->forceY_SP,  parameter->getParH(lev)->forceY_SP,  parameter->getParH(lev)->mem_size_real_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->forceZ_SP,  parameter->getParH(lev)->forceZ_SP,  parameter->getParH(lev)->mem_size_real_SP     , cudaMemcpyHostToDevice));

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
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityX   , parameter->getParD(lev)->velocityX   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityY   , parameter->getParD(lev)->velocityY   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->velocityZ   , parameter->getParD(lev)->velocityZ   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->rho  , parameter->getParD(lev)->rho  , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->pressure, parameter->getParD(lev)->pressure, parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
//sparse
void CudaMemoryManager::cudaAllocSP(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->typeOfGridNode           ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborX    ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborY    ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborZ    ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->rho          ), parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityX           ), parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityY           ), parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityZ           ), parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->pressure        ), parameter->getParH(lev)->mem_size_real_SP));
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->typeOfGridNode               ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborX        ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborY        ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborZ        ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->rho              ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityX               ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityY               ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityZ               ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressure            ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->distributions.f[0]           ), (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParD(lev)->mem_size_real_SP));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 4. * (double)parameter->getParH(lev)->mem_size_int_SP + 5. * (double)parameter->getParH(lev)->mem_size_real_SP + (double)parameter->getD3Qxx() * (double)parameter->getParH(lev)->mem_size_real_SP;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopySP(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->typeOfGridNode       ,  parameter->getParH(lev)->typeOfGridNode       ,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborX,  parameter->getParH(lev)->neighborX,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborY,  parameter->getParH(lev)->neighborY,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborZ,  parameter->getParH(lev)->neighborZ,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->rho      ,  parameter->getParH(lev)->rho      ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityX       ,  parameter->getParH(lev)->velocityX       ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityY       ,  parameter->getParH(lev)->velocityY       ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityZ       ,  parameter->getParH(lev)->velocityZ       ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->pressure    ,  parameter->getParH(lev)->pressure    ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeSP(int lev)
{
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->typeOfGridNode       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityX       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityY       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->velocityZ       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->rho      ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->pressure    ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborX));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborY));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborZ));
}
void CudaMemoryManager::cudaAllocF3SP(int lev)
{
    //Device
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->g6.g[0]), (unsigned long long)6*(unsigned long long)parameter->getParD(lev)->mem_size_real_SP));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)6 * (double)parameter->getParH(lev)->mem_size_real_SP;
    setMemsizeGPU(tmp, false);
}





















//Velo
void CudaMemoryManager::cudaAllocVeloBC(int lev)
{
	unsigned int mem_size_inflow_Q_k = sizeof(int)*parameter->getParH(lev)->velocityBC.numberOfBCnodes;
	unsigned int mem_size_inflow_Q_q = sizeof(real)*parameter->getParH(lev)->velocityBC.numberOfBCnodes;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.q27[0]),  parameter->getD3Qxx()*mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.k),                  mem_size_inflow_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.Vx),                 mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.Vy),                 mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.Vz),                 mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.deltaVz),            mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->velocityBC.RhoBC),              mem_size_inflow_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.q27[0]),      parameter->getD3Qxx()*mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.k),                      mem_size_inflow_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.Vx),                     mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.Vy),                     mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.Vz),                     mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->velocityBC.deltaVz),                mem_size_inflow_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_inflow_Q_k + 4. * (double)mem_size_inflow_Q_q + (double)parameter->getD3Qxx() * (double)mem_size_inflow_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyVeloBC(int lev)
{
	unsigned int mem_size_inflow_Q_k = sizeof(int)*parameter->getParH(lev)->velocityBC.numberOfBCnodes;
	unsigned int mem_size_inflow_Q_q = sizeof(real)*parameter->getParH(lev)->velocityBC.numberOfBCnodes;

	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.q27[0],  parameter->getParH(lev)->velocityBC.q27[0], parameter->getD3Qxx()* mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.k,       parameter->getParH(lev)->velocityBC.k,                  mem_size_inflow_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.Vx,      parameter->getParH(lev)->velocityBC.Vx,                 mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.Vy,      parameter->getParH(lev)->velocityBC.Vy,                 mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.Vz,      parameter->getParH(lev)->velocityBC.Vz,                 mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->velocityBC.deltaVz, parameter->getParH(lev)->velocityBC.deltaVz,            mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));

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
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBC.k),                 mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBC.kN),                mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->outflowBC.RhoBC),             mem_size_outflow_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBC.q27[0]),     parameter->getD3Qxx()* mem_size_outflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBC.k),                      mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBC.kN),                     mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->outflowBC.RhoBC),                  mem_size_outflow_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_outflow_Q_q + 2. * (double)mem_size_outflow_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_outflow_Q_q;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyOutflowBC(int lev)
{
	unsigned int mem_size_outflow_Q_k = sizeof(int)*parameter->getParH(lev)->outflowBC.numberOfBCnodes;
	unsigned int mem_size_outflow_Q_q = sizeof(real)*parameter->getParH(lev)->outflowBC.numberOfBCnodes;

	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBC.q27[0],  parameter->getParH(lev)->outflowBC.q27[0], parameter->getD3Qxx()* mem_size_outflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBC.k,       parameter->getParH(lev)->outflowBC.k,                  mem_size_outflow_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBC.kN,      parameter->getParH(lev)->outflowBC.kN,                 mem_size_outflow_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->outflowBC.RhoBC,   parameter->getParH(lev)->outflowBC.RhoBC,              mem_size_outflow_Q_q,  cudaMemcpyHostToDevice));
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
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->noSlipBC.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->noSlipBC.qread),             mem_size_Q_q_read ));//Geller
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->noSlipBC.valueQ),            mem_size_Q_value  ));//Geller

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->noSlipBC.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->noSlipBC.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyNoSlipBC(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->noSlipBC.numberOfBCnodes;
	unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->noSlipBC.numberOfBCnodes;

	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->noSlipBC.q27[0], parameter->getParH(lev)->noSlipBC.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->noSlipBC.k,      parameter->getParH(lev)->noSlipBC.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
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
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geometryBC.k),                 mem_size_Q_k      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBC.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geometryBC.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyGeomBC(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->geometryBC.numberOfBCnodes;
	unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->geometryBC.numberOfBCnodes;

	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBC.q27[0], parameter->getParH(lev)->geometryBC.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geometryBC.k,      parameter->getParH(lev)->geometryBC.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
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
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->pressureBC.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->pressureBC.kN),                mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->pressureBC.RhoBC),             mem_size_Q_q      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressureBC.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressureBC.k),                      mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressureBC.kN),                     mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->pressureBC.RhoBC),                  mem_size_Q_q     ));

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
        fx_t *= vf::lbm::constant::c2o1;
        fy_t *= vf::lbm::constant::c2o1;
        fz_t *= vf::lbm::constant::c2o1;
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
                                                     const unsigned int &memsizeFsRecv, int streamIndex)
{
    if (streamIndex == -1)
        checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].f[0],
						 parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].f[0],
						 parameter->getD3Qxx() * memsizeFsRecv,
						 cudaMemcpyHostToDevice));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].f[0],
                         parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].f[0],
                         parameter->getD3Qxx() * memsizeFsRecv,
                         cudaMemcpyHostToDevice,
                         parameter->getStreamManager()->getStream(streamIndex)));
}
void CudaMemoryManager::cudaCopyProcessNeighborXFsDH(int lev, unsigned int processNeighbor,
                                                     const unsigned int &memsizeFsSend, int streamIndex)
{
    if (streamIndex == -1)
    	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].f[0],
    								parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].f[0],
    								parameter->getD3Qxx() * memsizeFsSend,
    								cudaMemcpyDeviceToHost));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].f[0],
    								     parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].f[0],
    								     parameter->getD3Qxx() * memsizeFsSend,
    								     cudaMemcpyDeviceToHost,
                                         parameter->getStreamManager()->getStream(streamIndex)));
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
void CudaMemoryManager::cudaCopyProcessNeighborYFsHD(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsRecv,
                                                     int streamIndex)
{
    if (streamIndex == -1)
	    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].f[0],
								    parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].f[0],
								    parameter->getD3Qxx() * memsizeFsRecv,
								    cudaMemcpyHostToDevice));
    else
        checkCudaErrors(cudaMemcpyAsync(parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].f[0],
                        parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].f[0],
                        parameter->getD3Qxx() * memsizeFsRecv,
                        cudaMemcpyHostToDevice,
                        parameter->getStreamManager()->getStream(streamIndex)));
}
void CudaMemoryManager::cudaCopyProcessNeighborYFsDH(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsSend,
                                                     int streamIndex)
{
    if (streamIndex == -1)
	    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].f[0],
	    							parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].f[0],
	    							parameter->getD3Qxx() * memsizeFsSend,
	    							cudaMemcpyDeviceToHost));
    else
        checkCudaErrors(
            cudaMemcpyAsync(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].f[0],
                            parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].f[0],
                            parameter->getD3Qxx() * memsizeFsSend,
                            cudaMemcpyDeviceToHost, parameter->getStreamManager()->getStream(streamIndex)));
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
                                                     const unsigned int &memsizeFsRecv, int streamIndex)
{
    if (streamIndex == -1)
	    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].f[0],
	    							parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].f[0],
	    							parameter->getD3Qxx() * memsizeFsRecv,
	    							cudaMemcpyHostToDevice));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].f[0],
	    				                 parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].f[0],
	    				                 parameter->getD3Qxx() * memsizeFsRecv,
	    				                 cudaMemcpyHostToDevice,
                                         parameter->getStreamManager()->getStream(streamIndex)));
}
void CudaMemoryManager::cudaCopyProcessNeighborZFsDH(int lev, unsigned int processNeighbor,
                                                     const unsigned int &memsizeFsSend, int streamIndex)
{
    if (streamIndex == -1)
        checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].f[0],
	        					    parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].f[0],
	        					    parameter->getD3Qxx() * memsizeFsSend,
	        					    cudaMemcpyDeviceToHost));
    else
        checkCudaErrors( cudaMemcpyAsync(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].f[0],
	        						     parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].f[0],
	        						     parameter->getD3Qxx() * memsizeFsSend,
	        						     cudaMemcpyDeviceToHost,
                                         parameter->getStreamManager()->getStream(streamIndex)));
}
void CudaMemoryManager::cudaFreeProcessNeighborZ(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].index ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].f[0]     ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].index  ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].f[0]     ));
}
//////////////////////////////////////////////////////////////////////////
//Process Neighbors
//  3D domain decomposition F3
//  X  /////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocProcessNeighborF3X(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].index), parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].g[0]),  parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].memsizeGs));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].index), parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].g[0]),  parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].memsizeGs));

	//Device
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborF3X[processNeighbor].index), parameter->getParD(lev)->sendProcessNeighborF3X[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborF3X[processNeighbor].g[0]),  parameter->getParD(lev)->sendProcessNeighborF3X[processNeighbor].memsizeGs));
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborF3X[processNeighbor].index), parameter->getParD(lev)->recvProcessNeighborF3X[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborF3X[processNeighbor].g[0]),  parameter->getParD(lev)->recvProcessNeighborF3X[processNeighbor].memsizeGs));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp =
		(double)parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].memsizeIndex + (double)parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].memsizeGs +
		(double)parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].memsizeIndex + (double)parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].memsizeGs;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyProcessNeighborF3XIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors(cudaMemcpy(
		parameter->getParD(lev)->sendProcessNeighborF3X[processNeighbor].index,
		parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].index,
		parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].memsizeIndex,
		cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors(cudaMemcpy(
		parameter->getParD(lev)->recvProcessNeighborF3X[processNeighbor].index,
		parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].index,
		parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].memsizeIndex,
		cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborF3XFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors(cudaMemcpy(
		parameter->getParD(lev)->recvProcessNeighborF3X[processNeighbor].g[0],
		parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].g[0],
		parameter->getParD(lev)->recvProcessNeighborF3X[processNeighbor].memsizeGs,
		cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborF3XFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors(cudaMemcpy(
		parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].g[0],
		parameter->getParD(lev)->sendProcessNeighborF3X[processNeighbor].g[0],
		parameter->getParD(lev)->sendProcessNeighborF3X[processNeighbor].memsizeGs,
		cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeProcessNeighborF3X(int lev, unsigned int processNeighbor)
{
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].index));
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborF3X[processNeighbor].g[0]));
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].index));
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborF3X[processNeighbor].g[0]));
}
//  Y  /////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocProcessNeighborF3Y(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].index), parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].g[0]),  parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].memsizeGs));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].index), parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].g[0]),  parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].memsizeGs));

	//Device
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborF3Y[processNeighbor].index), parameter->getParD(lev)->sendProcessNeighborF3Y[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborF3Y[processNeighbor].g[0]),  parameter->getParD(lev)->sendProcessNeighborF3Y[processNeighbor].memsizeGs));
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborF3Y[processNeighbor].index), parameter->getParD(lev)->recvProcessNeighborF3Y[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborF3Y[processNeighbor].g[0]),  parameter->getParD(lev)->recvProcessNeighborF3Y[processNeighbor].memsizeGs));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp =
		(double)parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].memsizeIndex + (double)parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].memsizeGs +
		(double)parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].memsizeIndex + (double)parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].memsizeGs;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyProcessNeighborF3YIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors(cudaMemcpy(
		parameter->getParD(lev)->sendProcessNeighborF3Y[processNeighbor].index,
		parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].index,
		parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].memsizeIndex,
		cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors(cudaMemcpy(
		parameter->getParD(lev)->recvProcessNeighborF3Y[processNeighbor].index,
		parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].index,
		parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].memsizeIndex,
		cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborF3YFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors(cudaMemcpy(
		parameter->getParD(lev)->recvProcessNeighborF3Y[processNeighbor].g[0],
		parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].g[0],
		parameter->getParD(lev)->recvProcessNeighborF3Y[processNeighbor].memsizeGs,
		cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborF3YFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors(cudaMemcpy(
		parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].g[0],
		parameter->getParD(lev)->sendProcessNeighborF3Y[processNeighbor].g[0],
		parameter->getParD(lev)->sendProcessNeighborF3Y[processNeighbor].memsizeGs,
		cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeProcessNeighborF3Y(int lev, unsigned int processNeighbor)
{
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].index));
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborF3Y[processNeighbor].g[0]));
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].index));
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborF3Y[processNeighbor].g[0]));
}
//  Z  /////////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocProcessNeighborF3Z(int lev, unsigned int processNeighbor)
{
	//Host
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].index), parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].g[0]),  parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].memsizeGs));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].index), parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].g[0]),  parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].memsizeGs));

	//Device
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborF3Z[processNeighbor].index), parameter->getParD(lev)->sendProcessNeighborF3Z[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborF3Z[processNeighbor].g[0]),  parameter->getParD(lev)->sendProcessNeighborF3Z[processNeighbor].memsizeGs));
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborF3Z[processNeighbor].index), parameter->getParD(lev)->recvProcessNeighborF3Z[processNeighbor].memsizeIndex));
	checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborF3Z[processNeighbor].g[0]),  parameter->getParD(lev)->recvProcessNeighborF3Z[processNeighbor].memsizeGs));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tmp =
		(double)parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].memsizeIndex + (double)parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].memsizeGs +
		(double)parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].memsizeIndex + (double)parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].memsizeGs;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyProcessNeighborF3ZIndex(int lev, unsigned int processNeighbor)
{
	//copy send Index
	checkCudaErrors(cudaMemcpy(
		parameter->getParD(lev)->sendProcessNeighborF3Z[processNeighbor].index,
		parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].index,
		parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].memsizeIndex,
		cudaMemcpyHostToDevice));
	//copy recv Index
	checkCudaErrors(cudaMemcpy(
		parameter->getParD(lev)->recvProcessNeighborF3Z[processNeighbor].index,
		parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].index,
		parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].memsizeIndex,
		cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborF3ZFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors(cudaMemcpy(
		parameter->getParD(lev)->recvProcessNeighborF3Z[processNeighbor].g[0],
		parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].g[0],
		parameter->getParD(lev)->recvProcessNeighborF3Z[processNeighbor].memsizeGs,
		cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborF3ZFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors(cudaMemcpy(
		parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].g[0],
		parameter->getParD(lev)->sendProcessNeighborF3Z[processNeighbor].g[0],
		parameter->getParD(lev)->sendProcessNeighborF3Z[processNeighbor].memsizeGs,
		cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeProcessNeighborF3Z(int lev, unsigned int processNeighbor)
{
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].index));
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->sendProcessNeighborF3Z[processNeighbor].g[0]));
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].index));
	checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->recvProcessNeighborF3Z[processNeighbor].g[0]));
}
//////////////////////////////////////////////////////////////////////////

void CudaMemoryManager::cudaAllocNeighborWSB(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborInverse    ), parameter->getParH(lev)->mem_size_int_SP    ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborInverse        ), parameter->getParD(lev)->mem_size_int_SP    ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->mem_size_int_SP;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyNeighborWSB(int lev)
{
    //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborInverse,  parameter->getParH(lev)->neighborInverse,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeNeighborWSB(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborInverse));
}
//turbulent viscosity
void CudaMemoryManager::cudaAllocTurbulentViscosity(int lev)
{
    //Host
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->turbViscosity), parameter->getParH(lev)->mem_size_real_SP));
    //Debug
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gSij ), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gSDij), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDxvx), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDyvx), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDzvx), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDxvy), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDyvy), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDzvy), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDxvz), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDyvz), parameter->getParH(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDzvz), parameter->getParH(lev)->mem_size_real_SP));

    //Device
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->turbViscosity), parameter->getParD(lev)->mem_size_real_SP));
    //Debug
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gSij ), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gSDij), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDxvx), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDyvx), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDzvx), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDxvy), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDyvy), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDzvy), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDxvz), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDyvz), parameter->getParD(lev)->mem_size_real_SP));
    // checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDzvz), parameter->getParD(lev)->mem_size_real_SP));
    // //////////////////////////////////////////////////////////////////////////
    // double tmp = (double)parameter->getParH(lev)->mem_size_real_SP * 12.0;
    double tmp = (double)parameter->getParH(lev)->mem_size_real_SP;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyTurbulentViscosityHD(int lev)
{
    //copy host to device
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->turbViscosity, parameter->getParH(lev)->turbViscosity, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    //Debug
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gSij , parameter->getParH(lev)->gSij , parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gSDij, parameter->getParH(lev)->gSDij, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDxvx, parameter->getParH(lev)->gDxvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDyvx, parameter->getParH(lev)->gDyvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDzvx, parameter->getParH(lev)->gDzvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDxvy, parameter->getParH(lev)->gDxvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDyvy, parameter->getParH(lev)->gDyvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDzvy, parameter->getParH(lev)->gDzvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDxvz, parameter->getParH(lev)->gDxvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDyvz, parameter->getParH(lev)->gDyvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    // checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDzvz, parameter->getParH(lev)->gDzvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyTurbulentViscosityDH(int lev)
{
    //copy device to host
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->turbViscosity, parameter->getParD(lev)->turbViscosity, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    //Debug
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gSij , parameter->getParD(lev)->gSij , parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gSDij, parameter->getParD(lev)->gSDij, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDxvx, parameter->getParD(lev)->gDxvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDyvx, parameter->getParD(lev)->gDyvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDzvx, parameter->getParD(lev)->gDzvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDxvy, parameter->getParD(lev)->gDxvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDyvy, parameter->getParD(lev)->gDyvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDzvy, parameter->getParD(lev)->gDzvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDxvz, parameter->getParD(lev)->gDxvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDyvz, parameter->getParD(lev)->gDyvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDzvz, parameter->getParD(lev)->gDzvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeTurbulentViscosity(int lev)
{
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->turbViscosity));
    //Debug
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gSij ));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gSDij));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDxvx));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDyvx));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDzvx));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDxvy));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDyvy));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDzvy));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDxvz));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDyvz));
    // checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDzvz));
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
//median
void CudaMemoryManager::cudaAllocMedianSP(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->rho_SP_Med      ), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vx_SP_Med       ), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vy_SP_Med       ), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vz_SP_Med       ), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->press_SP_Med    ), parameter->getParH(lev)->mem_size_real_SP));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->rho_SP_Med          ), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vx_SP_Med           ), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vy_SP_Med           ), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vz_SP_Med           ), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->press_SP_Med        ), parameter->getParD(lev)->mem_size_real_SP));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 5. * (double)parameter->getParH(lev)->mem_size_real_SP;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyMedianSP(int lev)
{
    //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->rho_SP_Med  ,  parameter->getParH(lev)->rho_SP_Med  ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vx_SP_Med   ,  parameter->getParH(lev)->vx_SP_Med   ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vy_SP_Med   ,  parameter->getParH(lev)->vy_SP_Med   ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vz_SP_Med   ,  parameter->getParH(lev)->vz_SP_Med   ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->press_SP_Med,  parameter->getParH(lev)->press_SP_Med,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeMedianSP(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vx_SP_Med   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vy_SP_Med   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vz_SP_Med   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->rho_SP_Med  ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->press_SP_Med));
}
void CudaMemoryManager::cudaAllocMedianOut(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->rho_SP_Med_Out      ), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vx_SP_Med_Out       ), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vy_SP_Med_Out       ), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vz_SP_Med_Out       ), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->press_SP_Med_Out    ), parameter->getParH(lev)->mem_size_real_SP));
}
void CudaMemoryManager::cudaFreeMedianOut(int lev)
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
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->intCF.ICellCFC), parameter->getParH(lev)->mem_size_kCF  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->intCF.ICellCFF), parameter->getParH(lev)->mem_size_kCF  ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->intCF.ICellCFC), parameter->getParD(lev)->mem_size_kCF  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->intCF.ICellCFF), parameter->getParD(lev)->mem_size_kCF  ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 2. * (double)parameter->getParH(lev)->mem_size_kCF;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInterfaceCF(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->intCF.ICellCFC, parameter->getParH(lev)->intCF.ICellCFC, parameter->getParH(lev)->mem_size_kCF, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->intCF.ICellCFF, parameter->getParH(lev)->intCF.ICellCFF, parameter->getParH(lev)->mem_size_kCF, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeInterfaceCF(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->intCF.ICellCFC));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->intCF.ICellCFF));
}
//Interface FC
void CudaMemoryManager::cudaAllocInterfaceFC(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->intFC.ICellFCF), parameter->getParH(lev)->mem_size_kFC  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->intFC.ICellFCC), parameter->getParH(lev)->mem_size_kFC  ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->intFC.ICellFCF), parameter->getParD(lev)->mem_size_kFC  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->intFC.ICellFCC), parameter->getParD(lev)->mem_size_kFC  ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 2. * (double)parameter->getParH(lev)->mem_size_kFC;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInterfaceFC(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->intFC.ICellFCF, parameter->getParH(lev)->intFC.ICellFCF, parameter->getParH(lev)->mem_size_kFC, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->intFC.ICellFCC, parameter->getParH(lev)->intFC.ICellFCC, parameter->getParH(lev)->mem_size_kFC, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCheckInterfaceFCBulk(int lev)
{
    // only use for testing!
    size_t memsize = sizeof(uint) * parameter->getParH(lev)->intFCBulk.kFC;
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->intFCBulk.ICellFCC, parameter->getParH(lev)->intFCBulk.ICellFCC, memsize, cudaMemcpyDeviceToDevice));
    for (uint i = 0; i < parameter->getParH(lev)->intFCBulk.kFC; i++)
        printf("%d %d\n", i, parameter->getParH(lev)->intFCBulk.ICellFCC[i]);
}
void CudaMemoryManager::cudaFreeInterfaceFC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->intFC.ICellFCF));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->intFC.ICellFCC));
}
//Interface Offset CF
void CudaMemoryManager::cudaAllocInterfaceOffCF(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->offCF.xOffCF),   parameter->getParH(lev)->mem_size_kCF_off  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->offCF.yOffCF),   parameter->getParH(lev)->mem_size_kCF_off  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->offCF.zOffCF),   parameter->getParH(lev)->mem_size_kCF_off  ));
    getLastCudaError("Allocate host memory");
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->offCF.xOffCF),   parameter->getParD(lev)->mem_size_kCF_off  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->offCF.yOffCF),   parameter->getParD(lev)->mem_size_kCF_off  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->offCF.zOffCF),   parameter->getParD(lev)->mem_size_kCF_off  ));
    getLastCudaError("Allocate device memory");
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)parameter->getParH(lev)->mem_size_kCF_off;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInterfaceOffCF(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->offCF.xOffCF,   parameter->getParH(lev)->offCF.xOffCF,   parameter->getParH(lev)->mem_size_kCF_off, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->offCF.yOffCF,   parameter->getParH(lev)->offCF.yOffCF,   parameter->getParH(lev)->mem_size_kCF_off, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->offCF.zOffCF,   parameter->getParH(lev)->offCF.zOffCF,   parameter->getParH(lev)->mem_size_kCF_off, cudaMemcpyHostToDevice));
    getLastCudaError("Copy host memory to device");
}
void CudaMemoryManager::cudaFreeInterfaceOffCF(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->offCF.xOffCF));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->offCF.yOffCF));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->offCF.zOffCF));
}
//Interface Offset FC
void CudaMemoryManager::cudaAllocInterfaceOffFC(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->offFC.xOffFC),   parameter->getParH(lev)->mem_size_kFC_off  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->offFC.yOffFC),   parameter->getParH(lev)->mem_size_kFC_off  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->offFC.zOffFC),   parameter->getParH(lev)->mem_size_kFC_off  ));
    getLastCudaError("Allocate host memory");
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->offFC.xOffFC),   parameter->getParD(lev)->mem_size_kFC_off  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->offFC.yOffFC),   parameter->getParD(lev)->mem_size_kFC_off  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->offFC.zOffFC),   parameter->getParD(lev)->mem_size_kFC_off  ));
    getLastCudaError("Allocate device memory");
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)parameter->getParH(lev)->mem_size_kFC_off;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInterfaceOffFC(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->offFC.xOffFC,   parameter->getParH(lev)->offFC.xOffFC,   parameter->getParH(lev)->mem_size_kFC_off, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->offFC.yOffFC,   parameter->getParH(lev)->offFC.yOffFC,   parameter->getParH(lev)->mem_size_kFC_off, cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->offFC.zOffFC,   parameter->getParH(lev)->offFC.zOffFC,   parameter->getParH(lev)->mem_size_kFC_off, cudaMemcpyHostToDevice));
    getLastCudaError("Copy host memory to device");
}
void CudaMemoryManager::cudaFreeInterfaceOffFC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->offFC.xOffFC));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->offFC.yOffFC));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->offFC.zOffFC));
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
//Propeller Velocity
void CudaMemoryManager::cudaAllocVeloPropeller(int lev)
{
    unsigned int mem_size_Propeller_k = sizeof(int)*parameter->getParH(lev)->propellerBC.numberOfBCnodes;
    unsigned int mem_size_Propeller_q = sizeof(real)*parameter->getParH(lev)->propellerBC.numberOfBCnodes;

    //Host
    //checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->propellerBC.q27[0]),  parameter->getD3Qxx()*mem_size_Propeller_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->propellerBC.k),                  mem_size_Propeller_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->propellerBC.Vx),                 mem_size_Propeller_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->propellerBC.Vy),                 mem_size_Propeller_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->propellerBC.Vz),                 mem_size_Propeller_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->propellerBC.RhoBC),              mem_size_Propeller_q ));

    //Device
    //checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->propellerBC.q27[0]),      parameter->getD3Qxx()*mem_size_Propeller_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->propellerBC.k),                      mem_size_Propeller_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->propellerBC.Vx),                     mem_size_Propeller_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->propellerBC.Vy),                     mem_size_Propeller_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->propellerBC.Vz),                     mem_size_Propeller_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->propellerBC.RhoBC),                  mem_size_Propeller_q ));

    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_Propeller_k + 4. * (double)mem_size_Propeller_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyVeloPropeller(int lev)
{
    unsigned int mem_size_Propeller_k = sizeof(int)*parameter->getParH(lev)->propellerBC.numberOfBCnodes;
    unsigned int mem_size_Propeller_q = sizeof(real)*parameter->getParH(lev)->propellerBC.numberOfBCnodes;

    //checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->propellerBC.q27[0],  parameter->getParH(lev)->propellerBC.q27[0], parameter->getD3Qxx()* mem_size_Propeller_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->propellerBC.k,       parameter->getParH(lev)->propellerBC.k,                  mem_size_Propeller_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->propellerBC.Vx,      parameter->getParH(lev)->propellerBC.Vx,                 mem_size_Propeller_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->propellerBC.Vy,      parameter->getParH(lev)->propellerBC.Vy,                 mem_size_Propeller_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->propellerBC.Vz,      parameter->getParH(lev)->propellerBC.Vz,                 mem_size_Propeller_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->propellerBC.RhoBC,   parameter->getParH(lev)->propellerBC.RhoBC,              mem_size_Propeller_q,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeVeloPropeller(int lev)
{
    //checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->propellerBC.q27[0] ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->propellerBC.k      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->propellerBC.Vx     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->propellerBC.Vy     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->propellerBC.Vz     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->propellerBC.RhoBC  ));
}
//Measure Points
//void CudaMemoryManager::cudaAllocMeasurePoints(int lev, int i)
//{
//    //Host
//    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->MP[i].Vx),                 parameter->getParH(lev)->memSizerealMP ));
//    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->MP[i].Vy),                 parameter->getParH(lev)->memSizerealMP ));
//    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->MP[i].Vz),                 parameter->getParH(lev)->memSizerealMP ));
//    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->MP[i].Rho),                parameter->getParH(lev)->memSizerealMP ));
//
//    //Device
//    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->MP[i].Vx),                     parameter->getParD(lev)->memSizerealMP ));
//    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->MP[i].Vy),                     parameter->getParD(lev)->memSizerealMP ));
//    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->MP[i].Vz),                     parameter->getParD(lev)->memSizerealMP ));
//    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->MP[i].Rho),                    parameter->getParD(lev)->memSizerealMP ));
//}
//void CudaMemoryManager::cudaCopyMeasurePoints(int lev, int i)
//{
//    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->MP[i].Vx,      parameter->getParH(lev)->MP[i].Vx,           parameter->getParH(lev)->memSizerealMP,  cudaMemcpyHostToDevice));
//    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->MP[i].Vy,      parameter->getParH(lev)->MP[i].Vy,           parameter->getParH(lev)->memSizerealMP,  cudaMemcpyHostToDevice));
//    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->MP[i].Vz,      parameter->getParH(lev)->MP[i].Vz,           parameter->getParH(lev)->memSizerealMP,  cudaMemcpyHostToDevice));
//    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->MP[i].Rho,     parameter->getParH(lev)->MP[i].Rho,          parameter->getParH(lev)->memSizerealMP,  cudaMemcpyHostToDevice));
//}
//void CudaMemoryManager::cudaFreeMeasurePoints(int lev, int i)
//{
//    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->MP[i].Vx     ));
//    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->MP[i].Vy     ));
//    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->MP[i].Vz     ));
//    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->MP[i].Rho    ));
//}
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
void CudaMemoryManager::cudaAllocFsForCheckPointAndRestart(int lev)
{
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->distributions.f[0] ),           (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParH(lev)->mem_size_real_SP));
}
void CudaMemoryManager::cudaCopyFsForRestart(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->distributions.f[0],  parameter->getParH(lev)->distributions.f[0],     (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyFsForCheckPoint(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->distributions.f[0],  parameter->getParD(lev)->distributions.f[0],     (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeFsForCheckPointAndRestart(int lev)
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
//particles
void CudaMemoryManager::cudaAllocParticles(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.coordXlocal),        parameter->getParH(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.coordYlocal),        parameter->getParH(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.coordZlocal),        parameter->getParH(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.coordXabsolut),      parameter->getParH(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.coordYabsolut),      parameter->getParH(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.coordZabsolut),      parameter->getParH(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.veloX),              parameter->getParH(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.veloY),              parameter->getParH(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.veloZ),              parameter->getParH(lev)->plp.memSizerealAll  ));
    //checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.randomLocationInit), parameter->getParH(lev)->plp.memSizereal     ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.ID),                 parameter->getParH(lev)->plp.memSizeID          ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.cellBaseID),         parameter->getParH(lev)->plp.memSizeID          ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.timestep),           parameter->getParH(lev)->plp.memSizeTimestep    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.stuck),              parameter->getParH(lev)->plp.memSizeBool        ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->plp.hot),                parameter->getParH(lev)->plp.memSizeBoolBC      ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.coordXlocal),            parameter->getParD(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.coordYlocal),            parameter->getParD(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.coordZlocal),            parameter->getParD(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.coordXabsolut),          parameter->getParD(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.coordYabsolut),          parameter->getParD(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.coordZabsolut),          parameter->getParD(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.veloX),                  parameter->getParD(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.veloY),                  parameter->getParD(lev)->plp.memSizerealAll  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.veloZ),                  parameter->getParD(lev)->plp.memSizerealAll  ));
    //checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.randomLocationInit),     parameter->getParD(lev)->plp.memSizereal     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.ID),                     parameter->getParD(lev)->plp.memSizeID          ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.cellBaseID),             parameter->getParD(lev)->plp.memSizeID          ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.timestep),               parameter->getParD(lev)->plp.memSizeTimestep    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.stuck),                  parameter->getParD(lev)->plp.memSizeBool        ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->plp.hot),                    parameter->getParD(lev)->plp.memSizeBoolBC      ));

    //////////////////////////////////////////////////////////////////////////
    double tmp =
        (double)parameter->getParD(lev)->plp.memSizerealAll * (double)9.0 +
        (double)parameter->getParD(lev)->plp.memSizeID *      (double)2.0 +
        (double)parameter->getParD(lev)->plp.memSizeTimestep              +
        (double)parameter->getParD(lev)->plp.memSizeBool                  +
        (double)parameter->getParD(lev)->plp.memSizeBoolBC;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyParticles(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.coordXlocal,        parameter->getParD(lev)->plp.coordXlocal,        parameter->getParH(lev)->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.coordYlocal,        parameter->getParD(lev)->plp.coordYlocal,        parameter->getParH(lev)->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.coordZlocal,        parameter->getParD(lev)->plp.coordZlocal,        parameter->getParH(lev)->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.coordXabsolut,      parameter->getParD(lev)->plp.coordXabsolut,      parameter->getParH(lev)->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.coordYabsolut,      parameter->getParD(lev)->plp.coordYabsolut,      parameter->getParH(lev)->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.coordZabsolut,      parameter->getParD(lev)->plp.coordZabsolut,      parameter->getParH(lev)->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.veloX,              parameter->getParD(lev)->plp.veloX,              parameter->getParH(lev)->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.veloY,              parameter->getParD(lev)->plp.veloY,              parameter->getParH(lev)->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.veloZ,              parameter->getParD(lev)->plp.veloZ,              parameter->getParH(lev)->plp.memSizerealAll,  cudaMemcpyDeviceToHost));
    //checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.randomLocationInit, parameter->getParD(lev)->plp.randomLocationInit, parameter->getParH(lev)->plp.memSizereal,     cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.ID,                 parameter->getParD(lev)->plp.ID,                 parameter->getParH(lev)->plp.memSizeID,          cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.cellBaseID,         parameter->getParD(lev)->plp.cellBaseID,         parameter->getParH(lev)->plp.memSizeID,          cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->plp.timestep,           parameter->getParD(lev)->plp.timestep,           parameter->getParH(lev)->plp.memSizeTimestep,    cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeParticles(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.coordXlocal)       );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.coordYlocal)       );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.coordZlocal)       );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.coordXabsolut)     );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.coordYabsolut)     );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.coordZabsolut)     );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.veloX)             );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.veloY)             );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.veloZ)             );
    //checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.randomLocationInit));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.ID)                );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.cellBaseID)        );
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->plp.timestep)          );
}
//random values
void CudaMemoryManager::cudaAllocRandomValues()
{
    //Device
    checkCudaErrors( cudaMalloc((void**)(parameter->getRandomState()), (sizeof(curandState)*parameter->getParD(parameter->getFine())->plp.numberOfParticles) ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)(sizeof(curandState) * parameter->getParD(parameter->getFine())->plp.numberOfParticles);
    setMemsizeGPU(tmp, false);
}
//////////////////////////////////////////////////////////////////////////
//porous media
void CudaMemoryManager::cudaAllocPorousMedia(PorousMedia* pm, int lev)
{
    unsigned int mem_size_IDsPM = sizeof(unsigned int)*pm->getSizePM();
    unsigned int *tmpIDHost, *tmpIDDevice;
    //std::cout << "cudaMallocHost" << endl;
    //Host
    checkCudaErrors(cudaMallocHost((void**) &(tmpIDHost), mem_size_IDsPM));

    //std::cout << "cudaMalloc" << endl;
    //Device
    checkCudaErrors(cudaMalloc((void**) &(tmpIDDevice), mem_size_IDsPM));

    //std::cout << "set Host and Device arrays PM" << endl;
    //////////////////////////////////////////////////////////////////////////
    pm->setHostNodeIDsPM(tmpIDHost);
    pm->setDeviceNodeIDsPM(tmpIDDevice);
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_IDsPM;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyPorousMedia(PorousMedia* pm, int lev)
{
    unsigned int mem_size_IDsPM = sizeof(unsigned int)*pm->getSizePM();
    unsigned int *tmpIDHost   = pm->getHostNodeIDsPM();
    unsigned int *tmpIDDevice = pm->getDeviceNodeIDsPM();
    //////////////////////////////////////////////////////////////////////////
    checkCudaErrors(cudaMemcpy(tmpIDDevice, tmpIDHost, mem_size_IDsPM, cudaMemcpyHostToDevice));
    //////////////////////////////////////////////////////////////////////////
    pm->setDeviceNodeIDsPM(tmpIDDevice);
}
void CudaMemoryManager::cudaFreePorousMedia(PorousMedia* pm, int lev)
{
    checkCudaErrors(cudaFreeHost(pm->getHostNodeIDsPM()));
}
//////////////////////////////////////////////////////////////////////////
//advection diffusion
void CudaMemoryManager::cudaAllocConcentration(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Conc), parameter->getParH(lev)->mem_size_real_SP));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Conc), parameter->getParD(lev)->mem_size_real_SP));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->mem_size_real_SP;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyConcentrationDeviceToHost(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->Conc, parameter->getParD(lev)->Conc,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaCopyConcentrationHostToDevice(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Conc, parameter->getParH(lev)->Conc, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeConcentration(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Conc));
}
//////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocTempFs(int lev)
{
    //Device
    if (parameter->getDiffMod() == 7)
    {
        checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->distributionsAD7.f[0]), parameter->getDiffMod()*parameter->getParH(lev)->mem_size_real_SP));
    }
    else if (parameter->getDiffMod() == 27)
    {
        checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->distributionsAD27.f[0]), parameter->getDiffMod()*parameter->getParH(lev)->mem_size_real_SP));
    }
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)(parameter->getDiffMod() * parameter->getParH(lev)->mem_size_real_SP);
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
void CudaMemoryManager::cudaAllocMedianOutAD(int lev)
{
	//Host
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->rho_SP_Med_Out),   parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->vx_SP_Med_Out),    parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->vy_SP_Med_Out),    parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->vz_SP_Med_Out),    parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->press_SP_Med_Out), parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->Conc_Med_Out),     parameter->getParH(lev)->mem_size_real_SP));
}
void CudaMemoryManager::cudaFreeMedianOutAD(int lev)
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
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].index ),                  parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].f[0]  ),	parameter->getDiffMod() * parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].index ),                  parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].f[0]  ),   parameter->getDiffMod() * parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].index ),                      parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].index ),                      parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp =
        (double)parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeIndex +
        (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeFs +
        (double)parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeIndex +
        (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeFs;
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
                                parameter->getDiffMod() * parameter->getParD(lev)->recvProcessNeighborADX[processNeighbor].memsizeFs,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborADXFsDH(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].f[0],
                                parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].f[0],
                                parameter->getDiffMod() * parameter->getParD(lev)->sendProcessNeighborADX[processNeighbor].memsizeFs,
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
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].index ),                             parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].f[0]  ),   parameter->getDiffMod() * parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].index ),                             parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].f[0]  ),   parameter->getDiffMod() * parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].index ),                                 parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].index ),                                 parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp =
        (double)parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex +
        (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs +
        (double)parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex +
        (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs;
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
                                parameter->getDiffMod() * parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborADYFsDH(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].f[0],
                                parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].f[0],
                                parameter->getDiffMod() * parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs,
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
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].index ),                             parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].f[0]  ),   parameter->getDiffMod() * parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].index ),                             parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].f[0]  ),   parameter->getDiffMod() * parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs    ));

    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].index ),                                 parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].index ),                                 parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs    ));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp =
        (double)parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex +
        (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs +
        (double)parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex +
        (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs;
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
                                parameter->getDiffMod() * parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs,
                                cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborADZFsDH(int lev, unsigned int processNeighbor)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].f[0],
                                parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].f[0],
                                parameter->getDiffMod() * parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs,
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
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->dxxUx), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->dyyUy), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->dzzUz), parameter->getParH(lev)->mem_size_real_SP));
    //Device (spinning ship)
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->dxxUx), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->dyyUy), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->dzzUz), parameter->getParH(lev)->mem_size_real_SP));
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)parameter->getParH(lev)->mem_size_real_SP;
    setMemsizeGPU(tmp, false);
    //printf("Coord = %f MB",tmp/1000000.);
}
void CudaMemoryManager::cudaCopy2ndOrderDerivitivesIsoTestDH(int lev)
{
    //copy device to host
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->dxxUx, parameter->getParD(lev)->dxxUx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->dyyUy, parameter->getParD(lev)->dyyUy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->dzzUz, parameter->getParD(lev)->dzzUz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaCopy2ndOrderDerivitivesIsoTestHD(int lev)
{
    //copy host to device
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->dxxUx, parameter->getParH(lev)->dxxUx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->dyyUy, parameter->getParH(lev)->dyyUy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->dzzUz, parameter->getParH(lev)->dzzUz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));

}
void CudaMemoryManager::cudaFree2ndOrderDerivitivesIsoTest(int lev)
{
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->dxxUx));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->dyyUy));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->dzzUz));

}

void CudaMemoryManager::cudaAllocFluidNodeIndices(int lev) {
    uint mem_size_geo_fluid_nodes = sizeof(uint) * parameter->getParH(lev)->numberOfFluidNodes;
    // Host
    checkCudaErrors(cudaMallocHost((void **)&(parameter->getParH(lev)->fluidNodeIndices), mem_size_geo_fluid_nodes));
    // Device
    checkCudaErrors(cudaMalloc((void **)&(parameter->getParD(lev)->fluidNodeIndices), mem_size_geo_fluid_nodes));
    //////////////////////////////////////////////////////////////////////////
    setMemsizeGPU((double)mem_size_geo_fluid_nodes, false);
}

void CudaMemoryManager::cudaCopyFluidNodeIndices(int lev) {
    uint mem_size_geo_fluid_nodes = sizeof(uint) * parameter->getParH(lev)->numberOfFluidNodes;
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->fluidNodeIndices,
                               parameter->getParH(lev)->fluidNodeIndices,
                               mem_size_geo_fluid_nodes, cudaMemcpyHostToDevice));
}

void CudaMemoryManager::cudaFreeFluidNodeIndices(int lev) {
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->fluidNodeIndices));
}

void CudaMemoryManager::cudaAllocFluidNodeIndicesBorder(int lev) {
    uint mem_size_fluid_nodes_border = sizeof(uint) * parameter->getParH(lev)->numberOfFluidNodesBorder;
    // Host
    checkCudaErrors(
        cudaMallocHost((void **)&(parameter->getParH(lev)->fluidNodeIndicesBorder), mem_size_fluid_nodes_border));
    // Device
    checkCudaErrors(
        cudaMalloc((void **)&(parameter->getParD(lev)->fluidNodeIndicesBorder), mem_size_fluid_nodes_border));
    //////////////////////////////////////////////////////////////////////////
    setMemsizeGPU((double)mem_size_fluid_nodes_border, false);
}

void CudaMemoryManager::cudaCopyFluidNodeIndicesBorder(int lev) {
    uint mem_size_fluid_nodes_border = sizeof(uint) * parameter->getParH(lev)->numberOfFluidNodesBorder;
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->fluidNodeIndicesBorder,
                               parameter->getParH(lev)->fluidNodeIndicesBorder,
                               mem_size_fluid_nodes_border, cudaMemcpyHostToDevice));
}

void CudaMemoryManager::cudaFreeFluidNodeIndicesBorder(int lev) {
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->fluidNodeIndicesBorder));
}

////////////////////////////////////////////////////////////////////////////////////
//  ActuatorLine
///////////////////////////////////////////////////////////////////////////////

void CudaMemoryManager::cudaAllocBladeRadii(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeRadiiH, sizeof(real)*actuatorLine->getNBladeNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeRadiiD, sizeof(real)*actuatorLine->getNBladeNodes()) );

    setMemsizeGPU(sizeof(real)*actuatorLine->getNBladeNodes(), false);
}

void CudaMemoryManager::cudaCopyBladeRadiiHtoD(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeRadiiD, actuatorLine->bladeRadiiH, sizeof(real)*actuatorLine->getNBladeNodes(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyBladeRadiiDtoH(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeRadiiH, actuatorLine->bladeRadiiD, sizeof(real)*actuatorLine->getNBladeNodes(), cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeBladeRadii(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaFree(actuatorLine->bladeRadiiD) );

    checkCudaErrors( cudaFreeHost(actuatorLine->bladeRadiiH) );
}

void CudaMemoryManager::cudaAllocBladeCoords(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeCoordsXH, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeCoordsYH, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeCoordsZH, sizeof(real)*actuatorLine->getNNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeCoordsXD, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeCoordsYD, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeCoordsZD, sizeof(real)*actuatorLine->getNNodes()) );

    setMemsizeGPU(3.f*actuatorLine->getNNodes(), false);
}

void CudaMemoryManager::cudaCopyBladeCoordsHtoD(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeCoordsXD, actuatorLine->bladeCoordsXH, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeCoordsYD, actuatorLine->bladeCoordsYH, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeCoordsZD, actuatorLine->bladeCoordsZH, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyBladeCoordsDtoH(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeCoordsXH, actuatorLine->bladeCoordsXD, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeCoordsYH, actuatorLine->bladeCoordsYD, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeCoordsZH, actuatorLine->bladeCoordsZD, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeBladeCoords(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaFree(actuatorLine->bladeCoordsXD) );
    checkCudaErrors( cudaFree(actuatorLine->bladeCoordsYD) );
    checkCudaErrors( cudaFree(actuatorLine->bladeCoordsZD) );

    checkCudaErrors( cudaFreeHost(actuatorLine->bladeCoordsXH) );
    checkCudaErrors( cudaFreeHost(actuatorLine->bladeCoordsYH) );
    checkCudaErrors( cudaFreeHost(actuatorLine->bladeCoordsZH) );
}

void CudaMemoryManager::cudaAllocBladeIndices(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeIndicesH, sizeof(uint)*actuatorLine->getNNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeIndicesD, sizeof(uint)*actuatorLine->getNNodes()) );

    setMemsizeGPU(sizeof(uint)*actuatorLine->getNNodes(), false);
}

void CudaMemoryManager::cudaCopyBladeIndicesHtoD(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeIndicesD, actuatorLine->bladeIndicesH, sizeof(uint)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaFreeBladeIndices(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaFree(actuatorLine->bladeIndicesD) );

    checkCudaErrors( cudaFreeHost(actuatorLine->bladeIndicesH) );
}

void CudaMemoryManager::cudaAllocBladeVelocities(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeVelocitiesXH, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeVelocitiesYH, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeVelocitiesZH, sizeof(real)*actuatorLine->getNNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeVelocitiesXD, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeVelocitiesYD, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeVelocitiesZD, sizeof(real)*actuatorLine->getNNodes()) );

    setMemsizeGPU(3.*sizeof(real)*actuatorLine->getNNodes(), false);
}

void CudaMemoryManager::cudaCopyBladeVelocitiesHtoD(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeVelocitiesXD, actuatorLine->bladeVelocitiesXH, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeVelocitiesYD, actuatorLine->bladeVelocitiesYH, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeVelocitiesZD, actuatorLine->bladeVelocitiesZH, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyBladeVelocitiesDtoH(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeVelocitiesXH, actuatorLine->bladeVelocitiesXD, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeVelocitiesYH, actuatorLine->bladeVelocitiesYD, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeVelocitiesZH, actuatorLine->bladeVelocitiesZD, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeBladeVelocities(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaFree(actuatorLine->bladeVelocitiesXD) );
    checkCudaErrors( cudaFree(actuatorLine->bladeVelocitiesYD) );
    checkCudaErrors( cudaFree(actuatorLine->bladeVelocitiesZD) );

    checkCudaErrors( cudaFreeHost(actuatorLine->bladeVelocitiesXH) );
    checkCudaErrors( cudaFreeHost(actuatorLine->bladeVelocitiesYH) );
    checkCudaErrors( cudaFreeHost(actuatorLine->bladeVelocitiesZH) );
}

void CudaMemoryManager::cudaAllocBladeForces(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeForcesXH, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeForcesYH, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMallocHost((void**) &actuatorLine->bladeForcesZH, sizeof(real)*actuatorLine->getNNodes()) );

    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeForcesXD, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeForcesYD, sizeof(real)*actuatorLine->getNNodes()) );
    checkCudaErrors( cudaMalloc((void**) &actuatorLine->bladeForcesZD, sizeof(real)*actuatorLine->getNNodes()) );

    setMemsizeGPU(3.*sizeof(real)*actuatorLine->getNNodes(), false);
}

void CudaMemoryManager::cudaCopyBladeForcesHtoD(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeForcesXD, actuatorLine->bladeForcesXH, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeForcesYD, actuatorLine->bladeForcesYH, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeForcesZD, actuatorLine->bladeForcesZH, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaCopyBladeForcesDtoH(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeForcesXH, actuatorLine->bladeForcesXD, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeForcesYH, actuatorLine->bladeForcesYD, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(actuatorLine->bladeForcesZH, actuatorLine->bladeForcesZD, sizeof(real)*actuatorLine->getNNodes(), cudaMemcpyDeviceToHost) );
}

void CudaMemoryManager::cudaFreeBladeForces(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaFree(actuatorLine->bladeForcesXD) );
    checkCudaErrors( cudaFree(actuatorLine->bladeForcesYD) );
    checkCudaErrors( cudaFree(actuatorLine->bladeForcesZD) );

    checkCudaErrors( cudaFreeHost(actuatorLine->bladeForcesXH) );
    checkCudaErrors( cudaFreeHost(actuatorLine->bladeForcesYH) );
    checkCudaErrors( cudaFreeHost(actuatorLine->bladeForcesZH) );
}

void CudaMemoryManager::cudaAllocSphereIndices(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMallocHost((void**) &(actuatorLine->boundingSphereIndicesH), sizeof(int)*actuatorLine->getNIndices()));
    checkCudaErrors( cudaMalloc((void**) &(actuatorLine->boundingSphereIndicesD), sizeof(int)*actuatorLine->getNIndices()));
    setMemsizeGPU(sizeof(int)*actuatorLine->getNIndices(), false);
}

void CudaMemoryManager::cudaCopySphereIndicesHtoD(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaMemcpy(actuatorLine->boundingSphereIndicesD, actuatorLine->boundingSphereIndicesH, sizeof(int)*actuatorLine->getNIndices(), cudaMemcpyHostToDevice) );
}

void CudaMemoryManager::cudaFreeSphereIndices(ActuatorLine* actuatorLine)
{
    checkCudaErrors( cudaFreeHost(actuatorLine->boundingSphereIndicesH) );
    checkCudaErrors( cudaFree(actuatorLine->boundingSphereIndicesD) );
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
    size_t tmp = sizeof(real)*probe->getProbeStruct(level)->nArrays*probe->getProbeStruct(level)->nPoints;

    checkCudaErrors( cudaMallocHost((void**) &probe->getProbeStruct(level)->quantitiesArrayH, tmp) );
    if(probe->getHasDeviceQuantityArray())
    {
        checkCudaErrors( cudaMalloc    ((void**) &probe->getProbeStruct(level)->quantitiesArrayD, tmp) );
        setMemsizeGPU(1.f*tmp, false);
    }
}

void CudaMemoryManager::cudaCopyProbeQuantityArrayHtoD(Probe* probe, int level)
{
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->quantitiesArrayD, probe->getProbeStruct(level)->quantitiesArrayH, probe->getProbeStruct(level)->nArrays*sizeof(real)*probe->getProbeStruct(level)->nPoints, cudaMemcpyHostToDevice) );
}
void CudaMemoryManager::cudaCopyProbeQuantityArrayDtoH(Probe* probe, int level)
{
    checkCudaErrors( cudaMemcpy(probe->getProbeStruct(level)->quantitiesArrayH, probe->getProbeStruct(level)->quantitiesArrayD, probe->getProbeStruct(level)->nArrays*sizeof(real)*probe->getProbeStruct(level)->nPoints, cudaMemcpyDeviceToHost) );
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

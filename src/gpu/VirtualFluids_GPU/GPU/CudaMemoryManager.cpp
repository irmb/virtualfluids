#include "CudaMemoryManager.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Parameter/Parameter.h>

void CudaMemoryManager::cudaAllocFull(int lev)
{
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geo      ), parameter->getParH(lev)->mem_size_int  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->k        ), parameter->getParH(lev)->mem_size_int  ));
}
void CudaMemoryManager::cudaFreeFull(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geo   ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->k     ));
}
void CudaMemoryManager::cudaCopyPrint(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vx_SP   , parameter->getParD(lev)->vx_SP   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vy_SP   , parameter->getParD(lev)->vy_SP   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vz_SP   , parameter->getParD(lev)->vz_SP   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->rho_SP  , parameter->getParD(lev)->rho_SP  , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->press_SP, parameter->getParD(lev)->press_SP, parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
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
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordX_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordY_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordZ_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	//Device (spinning ship)
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordX_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordY_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordZ_SP      ), parameter->getParH(lev)->mem_size_real_SP  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (double)parameter->getParH(lev)->mem_size_real_SP;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyCoord(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordX_SP,  parameter->getParH(lev)->coordX_SP,  parameter->getParH(lev)->mem_size_real_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordY_SP,  parameter->getParH(lev)->coordY_SP,  parameter->getParH(lev)->mem_size_real_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordZ_SP,  parameter->getParH(lev)->coordZ_SP,  parameter->getParH(lev)->mem_size_real_SP     , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeCoord(int lev)
{
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->coordX_SP   ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->coordY_SP   ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->coordZ_SP   ));
}
//print
void CudaMemoryManager::cudaCopyDataToHost(int lev)
{
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vx_SP   , parameter->getParD(lev)->vx_SP   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vy_SP   , parameter->getParD(lev)->vy_SP   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vz_SP   , parameter->getParD(lev)->vz_SP   , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->rho_SP  , parameter->getParD(lev)->rho_SP  , parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->press_SP, parameter->getParD(lev)->press_SP, parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
//sparse
void CudaMemoryManager::cudaAllocSP(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geoSP           ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborX_SP    ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborY_SP    ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborZ_SP    ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->rho_SP          ), parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vx_SP           ), parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vy_SP           ), parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vz_SP           ), parameter->getParH(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->press_SP        ), parameter->getParH(lev)->mem_size_real_SP));
	//Device						 
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geoSP               ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborX_SP        ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborY_SP        ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborZ_SP        ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->rho_SP              ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vx_SP               ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vy_SP               ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vz_SP               ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->press_SP            ), parameter->getParD(lev)->mem_size_real_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->d0SP.f[0]           ), (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParD(lev)->mem_size_real_SP));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 4. * (double)parameter->getParH(lev)->mem_size_int_SP + 5. * (double)parameter->getParH(lev)->mem_size_real_SP + (double)parameter->getD3Qxx() * (double)parameter->getParH(lev)->mem_size_real_SP;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopySP(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geoSP       ,  parameter->getParH(lev)->geoSP       ,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborX_SP,  parameter->getParH(lev)->neighborX_SP,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborY_SP,  parameter->getParH(lev)->neighborY_SP,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborZ_SP,  parameter->getParH(lev)->neighborZ_SP,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->rho_SP      ,  parameter->getParH(lev)->rho_SP      ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vx_SP       ,  parameter->getParH(lev)->vx_SP       ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vy_SP       ,  parameter->getParH(lev)->vy_SP       ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vz_SP       ,  parameter->getParH(lev)->vz_SP       ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->press_SP    ,  parameter->getParH(lev)->press_SP    ,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeSP(int lev)
{
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->geoSP       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vx_SP       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vy_SP       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->vz_SP       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->rho_SP      ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->press_SP    ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborX_SP));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborY_SP));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborZ_SP));
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
	unsigned int mem_size_inflow_Q_k = sizeof(int)*parameter->getParH(lev)->Qinflow.kQ;
	unsigned int mem_size_inflow_Q_q = sizeof(real)*parameter->getParH(lev)->Qinflow.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qinflow.q27[0]),  parameter->getD3Qxx()*mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qinflow.k),                  mem_size_inflow_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qinflow.Vx),                 mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qinflow.Vy),                 mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qinflow.Vz),                 mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qinflow.deltaVz),            mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qinflow.RhoBC),              mem_size_inflow_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qinflow.q27[0]),      parameter->getD3Qxx()*mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qinflow.k),                      mem_size_inflow_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qinflow.Vx),                     mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qinflow.Vy),                     mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qinflow.Vz),                     mem_size_inflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qinflow.deltaVz),                mem_size_inflow_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_inflow_Q_k + 4. * (double)mem_size_inflow_Q_q + (double)parameter->getD3Qxx() * (double)mem_size_inflow_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyVeloBC(int lev)
{
	unsigned int mem_size_inflow_Q_k = sizeof(int)*parameter->getParH(lev)->Qinflow.kQ;
	unsigned int mem_size_inflow_Q_q = sizeof(real)*parameter->getParH(lev)->Qinflow.kQ;

	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qinflow.q27[0],  parameter->getParH(lev)->Qinflow.q27[0], parameter->getD3Qxx()* mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qinflow.k,       parameter->getParH(lev)->Qinflow.k,                  mem_size_inflow_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qinflow.Vx,      parameter->getParH(lev)->Qinflow.Vx,                 mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qinflow.Vy,      parameter->getParH(lev)->Qinflow.Vy,                 mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qinflow.Vz,      parameter->getParH(lev)->Qinflow.Vz,                 mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qinflow.deltaVz, parameter->getParH(lev)->Qinflow.deltaVz,            mem_size_inflow_Q_q,  cudaMemcpyHostToDevice));

}
void CudaMemoryManager::cudaFreeVeloBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qinflow.q27[0] ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qinflow.k      ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qinflow.Vx     ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qinflow.Vy     ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qinflow.Vz     ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qinflow.deltaVz));
}
//Press
void CudaMemoryManager::cudaAllocOutflowBC(int lev)
{
	unsigned int mem_size_outflow_Q_k = sizeof(int)*parameter->getParH(lev)->Qoutflow.kQ;
	unsigned int mem_size_outflow_Q_q = sizeof(real)*parameter->getParH(lev)->Qoutflow.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qoutflow.q27[0]), parameter->getD3Qxx()*mem_size_outflow_Q_q ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qoutflow.k),                 mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qoutflow.kN),                mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Qoutflow.RhoBC),             mem_size_outflow_Q_q ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qoutflow.q27[0]),     parameter->getD3Qxx()* mem_size_outflow_Q_q ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qoutflow.k),                      mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qoutflow.kN),                     mem_size_outflow_Q_k ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Qoutflow.RhoBC),                  mem_size_outflow_Q_q ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_outflow_Q_k + 2. * (double)mem_size_outflow_Q_q + (double)parameter->getD3Qxx()*(double)mem_size_outflow_Q_q;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyOutflowBC(int lev)
{
	unsigned int mem_size_outflow_Q_k = sizeof(int)*parameter->getParH(lev)->Qoutflow.kQ;
	unsigned int mem_size_outflow_Q_q = sizeof(real)*parameter->getParH(lev)->Qoutflow.kQ;

	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qoutflow.q27[0],  parameter->getParH(lev)->Qoutflow.q27[0], parameter->getD3Qxx()* mem_size_outflow_Q_q,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qoutflow.k,       parameter->getParH(lev)->Qoutflow.k,                  mem_size_outflow_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qoutflow.kN,      parameter->getParH(lev)->Qoutflow.kN,                 mem_size_outflow_Q_k,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Qoutflow.RhoBC,   parameter->getParH(lev)->Qoutflow.RhoBC,              mem_size_outflow_Q_q,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeOutflowBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qoutflow.q27[0] ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qoutflow.k      ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qoutflow.kN     ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Qoutflow.RhoBC  ));
}
//Wall
void CudaMemoryManager::cudaAllocWallBC(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->QWall.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->QWall.kQ;
	unsigned int mem_size_Q_value  = sizeof(long long)*parameter->getParH(lev)->QWall.kQ; //Geller
	unsigned int mem_size_Q_q_read = sizeof(real)*parameter->getParH(lev)->kQread;     //Geller

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QWall.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QWall.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QWall.qread),             mem_size_Q_q_read ));//Geller
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QWall.valueQ),            mem_size_Q_value  ));//Geller

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QWall.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QWall.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyWallBC(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QWall.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QWall.kQ;

	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QWall.q27[0], parameter->getParH(lev)->QWall.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QWall.k,      parameter->getParH(lev)->QWall.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeWallBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QWall.q27[0]));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QWall.k));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QWall.valueQ));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QWall.qread));
}
//Geometrie
void CudaMemoryManager::cudaAllocGeomBC(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->QGeom.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->QGeom.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeom.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeom.k),                 mem_size_Q_k      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeom.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeom.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyGeomBC(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QGeom.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QGeom.kQ;

	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeom.q27[0], parameter->getParH(lev)->QGeom.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeom.k,      parameter->getParH(lev)->QGeom.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeGeomBC(int lev)
{
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeom.q27[0]));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeom.k));
}
//Press
void CudaMemoryManager::cudaAllocPress(int lev)
{
	unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->QPress.kQ;
	unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->QPress.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPress.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPress.k),                 mem_size_Q_k      ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPress.kN),                mem_size_Q_k      )); 
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPress.RhoBC),             mem_size_Q_q      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPress.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPress.k),                      mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPress.kN),                     mem_size_Q_k     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPress.RhoBC),                  mem_size_Q_q     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = 2. * (double)mem_size_Q_k + (double)mem_size_Q_q + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyPress(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QPress.kQ;
	unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QPress.kQ;

	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPress.q27[0], parameter->getParH(lev)->QPress.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPress.k,      parameter->getParH(lev)->QPress.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPress.kN,     parameter->getParH(lev)->QPress.kN,                 mem_size_Q_k,       cudaMemcpyHostToDevice)); 
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPress.RhoBC,  parameter->getParH(lev)->QPress.RhoBC,              mem_size_Q_q,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreePress(int lev)
{
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPress.q27[0]));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPress.k));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPress.kN));
	checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPress.RhoBC));
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
	double tmp = (real)mem_size;
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
	double tmp = (real)mem_size;
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
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].index ),                  parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].index ),                  parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].index ),                      parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].index ),                      parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].memsizeIndex ));
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
void CudaMemoryManager::cudaCopyProcessNeighborXFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].f[0], 
								parameter->getParH(lev)->recvProcessNeighborX[processNeighbor].f[0], 
								parameter->getD3Qxx() * parameter->getParD(lev)->recvProcessNeighborX[processNeighbor].memsizeFs, 
								cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborXFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborX[processNeighbor].f[0], 
								parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].f[0], 
								parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighborX[processNeighbor].memsizeFs, 
								cudaMemcpyDeviceToHost));
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
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].index ),                  parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].index ),                  parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].index ),                      parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].index ),                      parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].memsizeIndex ));
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
void CudaMemoryManager::cudaCopyProcessNeighborYFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].f[0], 
								parameter->getParH(lev)->recvProcessNeighborY[processNeighbor].f[0], 
								parameter->getD3Qxx() * parameter->getParD(lev)->recvProcessNeighborY[processNeighbor].memsizeFs, 
								cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborYFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborY[processNeighbor].f[0], 
								parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].f[0], 
								parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighborY[processNeighbor].memsizeFs, 
								cudaMemcpyDeviceToHost));
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
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].index ),                  parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].index ),                  parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].memsizeFs    ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].index ),                      parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].memsizeIndex ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].memsizeFs    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].index ),                      parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].memsizeIndex ));
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
void CudaMemoryManager::cudaCopyProcessNeighborZFsHD(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].f[0], 
								parameter->getParH(lev)->recvProcessNeighborZ[processNeighbor].f[0], 
								parameter->getD3Qxx() * parameter->getParD(lev)->recvProcessNeighborZ[processNeighbor].memsizeFs, 
								cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyProcessNeighborZFsDH(int lev, unsigned int processNeighbor)
{
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->sendProcessNeighborZ[processNeighbor].f[0], 
								parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].f[0], 
								parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighborZ[processNeighbor].memsizeFs, 
								cudaMemcpyDeviceToHost));
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
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborWSB_SP    ), parameter->getParH(lev)->mem_size_int_SP    ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborWSB_SP        ), parameter->getParD(lev)->mem_size_int_SP    ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->mem_size_int_SP;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyNeighborWSB(int lev)
{
    //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborWSB_SP,  parameter->getParH(lev)->neighborWSB_SP,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeNeighborWSB(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->neighborWSB_SP));
}
//turbulent viscosity
void CudaMemoryManager::cudaAllocTurbulentViscosity(int lev)
{
    //Host
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->turbViscosity), parameter->getParH(lev)->mem_size_real_SP));
    //Debug
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gSij ), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gSDij), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDxvx), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDyvx), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDzvx), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDxvy), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDyvy), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDzvy), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDxvz), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDyvz), parameter->getParH(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMallocHost((void**) &(parameter->getParH(lev)->gDzvz), parameter->getParH(lev)->mem_size_real_SP));
    
    //Device
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->turbViscosity), parameter->getParD(lev)->mem_size_real_SP));
    //Debug
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gSij ), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gSDij), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDxvx), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDyvx), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDzvx), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDxvy), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDyvy), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDzvy), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDxvz), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDyvz), parameter->getParD(lev)->mem_size_real_SP));
    checkCudaErrors(cudaMalloc((void**) &(parameter->getParD(lev)->gDzvz), parameter->getParD(lev)->mem_size_real_SP));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->mem_size_real_SP;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyTurbulentViscosityHD(int lev)
{
    //copy host to device
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->turbViscosity, parameter->getParH(lev)->turbViscosity, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    //Debug
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gSij , parameter->getParH(lev)->gSij , parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gSDij, parameter->getParH(lev)->gSDij, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDxvx, parameter->getParH(lev)->gDxvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDyvx, parameter->getParH(lev)->gDyvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDzvx, parameter->getParH(lev)->gDzvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDxvy, parameter->getParH(lev)->gDxvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDyvy, parameter->getParH(lev)->gDyvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDzvy, parameter->getParH(lev)->gDzvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDxvz, parameter->getParH(lev)->gDxvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDyvz, parameter->getParH(lev)->gDyvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(parameter->getParD(lev)->gDzvz, parameter->getParH(lev)->gDzvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyTurbulentViscosityDH(int lev)
{
    //copy device to host
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->turbViscosity, parameter->getParD(lev)->turbViscosity, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    //Debug
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gSij , parameter->getParD(lev)->gSij , parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gSDij, parameter->getParD(lev)->gSDij, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDxvx, parameter->getParD(lev)->gDxvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDyvx, parameter->getParD(lev)->gDyvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDzvx, parameter->getParD(lev)->gDzvx, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDxvy, parameter->getParD(lev)->gDxvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDyvy, parameter->getParD(lev)->gDyvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDzvy, parameter->getParD(lev)->gDzvy, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDxvz, parameter->getParD(lev)->gDxvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDyvz, parameter->getParD(lev)->gDyvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(parameter->getParH(lev)->gDzvz, parameter->getParD(lev)->gDzvz, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeTurbulentViscosity(int lev)
{
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->turbViscosity));
    //Debug
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gSij ));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gSDij));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDxvx));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDyvx));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDzvx));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDxvy));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDyvy));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDzvy));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDxvz));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDyvz));
    checkCudaErrors(cudaFreeHost(parameter->getParH(lev)->gDzvz));
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
    unsigned int mem_size_inlet_Q_k = sizeof(int)*parameter->getParH(lev)->QInlet.kQ;
    unsigned int mem_size_inlet_Q_q = sizeof(real)*parameter->getParH(lev)->QInlet.kQ;
    
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
    unsigned int mem_size_inlet_Q_k = sizeof(int)*parameter->getParH(lev)->QInlet.kQ;
    unsigned int mem_size_inlet_Q_q = sizeof(real)*parameter->getParH(lev)->QInlet.kQ;
    
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
    unsigned int mem_size_outlet_Q_k = sizeof(int)*parameter->getParH(lev)->QOutlet.kQ;
    unsigned int mem_size_outlet_Q_q = sizeof(real)*parameter->getParH(lev)->QOutlet.kQ;
    
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
    unsigned int mem_size_outlet_Q_k = sizeof(int)*parameter->getParH(lev)->QOutlet.kQ;
    unsigned int mem_size_outlet_Q_q = sizeof(real)*parameter->getParH(lev)->QOutlet.kQ;
    
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
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->QGeom.kQ;
    
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeom.Vx),  mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeom.Vy),  mem_size_Q_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeom.Vz),  mem_size_Q_q ));
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeom.Vx),      mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeom.Vy),      mem_size_Q_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeom.Vz),      mem_size_Q_q ));
    
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3. * (double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyGeomValuesBC(int lev)
{
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QGeom.kQ;
    
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeom.Vx, parameter->getParH(lev)->QGeom.Vx,  mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeom.Vy, parameter->getParH(lev)->QGeom.Vy,  mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeom.Vz, parameter->getParH(lev)->QGeom.Vz,  mem_size_Q_q,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeGeomValuesBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeom.Vx));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeom.Vy));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeom.Vz));
}
//Geometrie inkl. Normale f¸r Slip
void CudaMemoryManager::cudaAllocGeomNormals(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->QGeomNormalX.kQ;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->QGeomNormalX.kQ;
    
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeomNormalX.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeomNormalX.k),                 mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeomNormalY.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeomNormalY.k),                 mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeomNormalZ.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeomNormalZ.k),                 mem_size_Q_k      ));
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeomNormalX.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeomNormalX.k),                      mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeomNormalY.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeomNormalY.k),                      mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeomNormalZ.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeomNormalZ.k),                      mem_size_Q_k     ));
    
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3.0 * (double)mem_size_Q_k + 3.0 * (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyGeomNormals(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QGeomNormalX.kQ;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QGeomNormalX.kQ;
    
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeomNormalX.q27[0], parameter->getParH(lev)->QGeomNormalX.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeomNormalX.k,      parameter->getParH(lev)->QGeomNormalX.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeomNormalY.q27[0], parameter->getParH(lev)->QGeomNormalY.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeomNormalY.k,      parameter->getParH(lev)->QGeomNormalY.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeomNormalZ.q27[0], parameter->getParH(lev)->QGeomNormalZ.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QGeomNormalZ.k,      parameter->getParH(lev)->QGeomNormalZ.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeGeomNormals(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeomNormalX.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeomNormalX.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeomNormalY.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeomNormalY.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeomNormalZ.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QGeomNormalZ.k));
}
//Geometrie inkl. Normale f¸r Inflow
void CudaMemoryManager::cudaAllocInflowNormals(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->QInflowNormalX.kQ;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->QInflowNormalX.kQ;
    
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInflowNormalX.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInflowNormalX.k),                 mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInflowNormalY.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInflowNormalY.k),                 mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInflowNormalZ.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QInflowNormalZ.k),                 mem_size_Q_k      ));
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInflowNormalX.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInflowNormalX.k),                      mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInflowNormalY.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInflowNormalY.k),                      mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInflowNormalZ.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QInflowNormalZ.k),                      mem_size_Q_k     ));
    
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3.0 * (double)mem_size_Q_k + 3.0 * (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyInflowNormals(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QInflowNormalX.kQ;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QInflowNormalX.kQ;
    
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInflowNormalX.q27[0], parameter->getParH(lev)->QInflowNormalX.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInflowNormalX.k,      parameter->getParH(lev)->QInflowNormalX.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInflowNormalY.q27[0], parameter->getParH(lev)->QInflowNormalY.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInflowNormalY.k,      parameter->getParH(lev)->QInflowNormalY.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInflowNormalZ.q27[0], parameter->getParH(lev)->QInflowNormalZ.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QInflowNormalZ.k,      parameter->getParH(lev)->QInflowNormalZ.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeInflowNormals(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInflowNormalX.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInflowNormalX.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInflowNormalY.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInflowNormalY.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInflowNormalZ.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QInflowNormalZ.k));
}
//Geometrie inkl. Normale f¸r Outflow
void CudaMemoryManager::cudaAllocOutflowNormals(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->QOutflowNormalX.kQ;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->QOutflowNormalX.kQ;
    
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutflowNormalX.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutflowNormalX.k),                 mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutflowNormalY.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutflowNormalY.k),                 mem_size_Q_k      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutflowNormalZ.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QOutflowNormalZ.k),                 mem_size_Q_k      ));
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutflowNormalX.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutflowNormalX.k),                      mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutflowNormalY.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutflowNormalY.k),                      mem_size_Q_k     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutflowNormalZ.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QOutflowNormalZ.k),                      mem_size_Q_k     ));
    
    //////////////////////////////////////////////////////////////////////////
    double tmp = 3.0 * (double)mem_size_Q_k + 3.0 * (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyOutflowNormals(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QOutflowNormalX.kQ;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QOutflowNormalX.kQ;
    
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutflowNormalX.q27[0], parameter->getParH(lev)->QOutflowNormalX.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutflowNormalX.k,      parameter->getParH(lev)->QOutflowNormalX.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutflowNormalY.q27[0], parameter->getParH(lev)->QOutflowNormalY.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutflowNormalY.k,      parameter->getParH(lev)->QOutflowNormalY.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutflowNormalZ.q27[0], parameter->getParH(lev)->QOutflowNormalZ.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QOutflowNormalZ.k,      parameter->getParH(lev)->QOutflowNormalZ.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeOutflowNormals(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutflowNormalX.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutflowNormalX.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutflowNormalY.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutflowNormalY.k));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutflowNormalZ.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QOutflowNormalZ.k));
}
//Slip
void CudaMemoryManager::cudaAllocSlipBC(int lev)
{
    unsigned int mem_size_Q_k      = sizeof(int)*parameter->getParH(lev)->QSlip.kQ;
    unsigned int mem_size_Q_q      = sizeof(real)*parameter->getParH(lev)->QSlip.kQ;
    //unsigned int mem_size_Q_value  = sizeof(long long)*parameter->getParH(lev)->QSlip.kQ; //Geller
    //unsigned int mem_size_Q_q_read = sizeof(real)*parameter->getParH(lev)->kSlipQread;     //Geller
    
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QSlip.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QSlip.k),                 mem_size_Q_k      ));
    //checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QSlip.qread),             mem_size_Q_q_read ));//Geller
    //checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QSlip.valueQ),            mem_size_Q_value  ));//Geller
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QSlip.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QSlip.k),                      mem_size_Q_k     ));
    
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopySlipBC(int lev)
{
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QSlip.kQ;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QSlip.kQ;
    
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QSlip.q27[0], parameter->getParH(lev)->QSlip.q27[0], parameter->getD3Qxx()* mem_size_Q_q,       cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QSlip.k,      parameter->getParH(lev)->QSlip.k,                  mem_size_Q_k,       cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeSlipBC(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QSlip.q27[0]));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QSlip.k));
    //checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QSlip.valueQ));
    //checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QSlip.qread));
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
    double tmp = (double)size;
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
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QpressX0.kQ;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QpressX0.kQ;
    
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
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QpressX0.kQ;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QpressX0.kQ;
    
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
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QpressX1.kQ;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QpressX1.kQ;
    
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
    unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QpressX1.kQ;
    unsigned int mem_size_Q_q = sizeof(real)*parameter->getParH(lev)->QpressX1.kQ;
    
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
    unsigned int mem_size_Propeller_k = sizeof(int)*parameter->getParH(lev)->QPropeller.kQ;
    unsigned int mem_size_Propeller_q = sizeof(real)*parameter->getParH(lev)->QPropeller.kQ;
    
    //Host
    //checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPropeller.q27[0]),  parameter->getD3Qxx()*mem_size_Propeller_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPropeller.k),                  mem_size_Propeller_k ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPropeller.Vx),                 mem_size_Propeller_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPropeller.Vy),                 mem_size_Propeller_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPropeller.Vz),                 mem_size_Propeller_q ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QPropeller.RhoBC),              mem_size_Propeller_q ));
    
    //Device
    //checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPropeller.q27[0]),      parameter->getD3Qxx()*mem_size_Propeller_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPropeller.k),                      mem_size_Propeller_k ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPropeller.Vx),                     mem_size_Propeller_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPropeller.Vy),                     mem_size_Propeller_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPropeller.Vz),                     mem_size_Propeller_q ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QPropeller.RhoBC),                  mem_size_Propeller_q ));
    
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)mem_size_Propeller_k + 4. * (double)mem_size_Propeller_q;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyVeloPropeller(int lev)
{
    unsigned int mem_size_Propeller_k = sizeof(int)*parameter->getParH(lev)->QPropeller.kQ;
    unsigned int mem_size_Propeller_q = sizeof(real)*parameter->getParH(lev)->QPropeller.kQ;
    
    //checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPropeller.q27[0],  parameter->getParH(lev)->QPropeller.q27[0], parameter->getD3Qxx()* mem_size_Propeller_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPropeller.k,       parameter->getParH(lev)->QPropeller.k,                  mem_size_Propeller_k,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPropeller.Vx,      parameter->getParH(lev)->QPropeller.Vx,                 mem_size_Propeller_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPropeller.Vy,      parameter->getParH(lev)->QPropeller.Vy,                 mem_size_Propeller_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPropeller.Vz,      parameter->getParH(lev)->QPropeller.Vz,                 mem_size_Propeller_q,  cudaMemcpyHostToDevice));
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->QPropeller.RhoBC,   parameter->getParH(lev)->QPropeller.RhoBC,              mem_size_Propeller_q,  cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeVeloPropeller(int lev)
{
    //checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPropeller.q27[0] ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPropeller.k      ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPropeller.Vx     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPropeller.Vy     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPropeller.Vz     ));
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->QPropeller.RhoBC  ));
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
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->kMP),                        parameter->getParH(lev)->memSizeIntkMP     ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->VxMP),                    parameter->getParH(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->VyMP),                    parameter->getParH(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->VzMP),                    parameter->getParH(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->RhoMP),                    parameter->getParH(lev)->memSizerealkMP ));
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->kMP),                            parameter->getParD(lev)->memSizeIntkMP     ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->VxMP),                        parameter->getParD(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->VyMP),                        parameter->getParD(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->VzMP),                        parameter->getParD(lev)->memSizerealkMP ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->RhoMP),                        parameter->getParD(lev)->memSizerealkMP ));
    
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
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->d0SP.f[0] ),           (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParH(lev)->mem_size_real_SP));
}
void CudaMemoryManager::cudaCopyFsForRestart(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->d0SP.f[0],  parameter->getParH(lev)->d0SP.f[0],     (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyFsForCheckPoint(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->d0SP.f[0],  parameter->getParD(lev)->d0SP.f[0],     (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeFsForCheckPointAndRestart(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->d0SP.f[0]));
}
//DragLift
void CudaMemoryManager::cudaAllocDragLift(int lev, int numofelem)
{
    unsigned int mem_size = sizeof(double)*numofelem;
    
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPreX), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPreY), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPreZ), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPostX), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPostY), mem_size  ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->DragPostZ), mem_size  ));
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->DragPreX), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->DragPreY), mem_size  ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->DragPreZ), mem_size  ));
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
    double tmp = 5. * (real)mem_size;
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
    double tmp = 7. * (real)mem_size;
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
    double tmp = 7. * (real)mem_size;
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
    double tmp = 3. * (real)mem_size;
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
    double tmp = (double)parameter->getParD(lev)->plp.memSizerealAll * (double)9.0 + (double)parameter->getParD(lev)->plp.memSizeID * (double)2.0 + (double)parameter->getParD(lev)->plp.memSizeTimestep
    + (double)parameter->getParD(lev)->plp.memSizeBool + (double)parameter->getParD(lev)->plp.memSizeBoolBC;
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
void CudaMemoryManager::cudaAllocConc(int lev)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->Conc), parameter->getParH(lev)->mem_size_real_SP));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->Conc), parameter->getParD(lev)->mem_size_real_SP));
}
void CudaMemoryManager::cudaCopyConcDH(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->Conc, parameter->getParD(lev)->Conc,  parameter->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaCopyConcHD(int lev)
{
    checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->Conc, parameter->getParH(lev)->Conc, parameter->getParH(lev)->mem_size_real_SP, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeConc(int lev)
{
    checkCudaErrors( cudaFreeHost(parameter->getParH(lev)->Conc));
}
//////////////////////////////////////////////////////////////////////////
void CudaMemoryManager::cudaAllocTempFs(int lev)
{
    //Device
    if (parameter->getDiffMod() == 7)
    {
        checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->d7.f[0]), parameter->getDiffMod()*parameter->getParH(lev)->mem_size_real_SP));
    }
    else if (parameter->getDiffMod() == 27)
    {
        checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->d27.f[0]), parameter->getDiffMod()*parameter->getParH(lev)->mem_size_real_SP));
    }
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

//////////////////////////////////////////////////////////////////////////
//Process Neighbors
//1D domain decomposition
void CudaMemoryManager::cudaAllocProcessNeighbor(int lev, unsigned int processNeighbor)
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].index ),                  parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].index ),                  parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].f[0]  ),     parameter->getD3Qxx() * parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].memsizeFs    ));
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].index ),                      parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->sendProcessNeighbor[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].index ),                      parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].f[0]  ),         parameter->getD3Qxx() * parameter->getParD(lev)->recvProcessNeighbor[processNeighbor].memsizeFs    ));
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].memsizeIndex + (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->sendProcessNeighbor[processNeighbor].memsizeFs +
    (double)parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].memsizeIndex + (double)parameter->getD3Qxx()*(double)parameter->getParH(lev)->recvProcessNeighbor[processNeighbor].memsizeFs;
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
    double tmp = (double)parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeIndex + (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->sendProcessNeighborADX[processNeighbor].memsizeFs +
    (double)parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeIndex + (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->recvProcessNeighborADX[processNeighbor].memsizeFs;
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
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].index ),                  parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].f[0]  ),   parameter->getDiffMod() * parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].index ),                  parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].f[0]  ),   parameter->getDiffMod() * parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs    ));
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].index ),                      parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].index ),                      parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs    ));
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeIndex + (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->sendProcessNeighborADY[processNeighbor].memsizeFs +
    (double)parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeIndex + (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->recvProcessNeighborADY[processNeighbor].memsizeFs;
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
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].index ),                  parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].f[0]  ),   parameter->getDiffMod() * parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].index ),                  parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].f[0]  ),   parameter->getDiffMod() * parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs    ));
    
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].index ),                      parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs    ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].index ),                      parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex ));
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].f[0]  ),       parameter->getDiffMod() * parameter->getParD(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs    ));
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeIndex + (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->sendProcessNeighborADZ[processNeighbor].memsizeFs +
    (double)parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeIndex + (double)parameter->getDiffMod()*(double)parameter->getParH(lev)->recvProcessNeighborADZ[processNeighbor].memsizeFs;
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
























std::shared_ptr<CudaMemoryManager> CudaMemoryManager::make(std::shared_ptr<Parameter> parameter)
{
    return std::shared_ptr<CudaMemoryManager>(new CudaMemoryManager(parameter));
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
	return memsizeGPU;
}

CudaMemoryManager::CudaMemoryManager(std::shared_ptr<Parameter> parameter)
{
    this->parameter = parameter;
}

CudaMemoryManager::CudaMemoryManager(const CudaMemoryManager&)
{

}
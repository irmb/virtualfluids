#include "CudaMemoryManager.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Parameter/Parameter.h>

void CudaMemoryManager::cudaAllocCoord(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordX_SP      ), parameter->getParH(lev)->mem_size_doubflo_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordY_SP      ), parameter->getParH(lev)->mem_size_doubflo_SP  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->coordZ_SP      ), parameter->getParH(lev)->mem_size_doubflo_SP  ));
	//Device (spinning ship)
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordX_SP      ), parameter->getParH(lev)->mem_size_doubflo_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordY_SP      ), parameter->getParH(lev)->mem_size_doubflo_SP  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->coordZ_SP      ), parameter->getParH(lev)->mem_size_doubflo_SP  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (double)parameter->getParH(lev)->mem_size_doubflo_SP;
	parameter->setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyCoord(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordX_SP,  parameter->getParH(lev)->coordX_SP,  parameter->getParH(lev)->mem_size_doubflo_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordY_SP,  parameter->getParH(lev)->coordY_SP,  parameter->getParH(lev)->mem_size_doubflo_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->coordZ_SP,  parameter->getParH(lev)->coordZ_SP,  parameter->getParH(lev)->mem_size_doubflo_SP     , cudaMemcpyHostToDevice));
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
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vx_SP   , parameter->getParD(lev)->vx_SP   , parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vy_SP   , parameter->getParD(lev)->vy_SP   , parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->vz_SP   , parameter->getParD(lev)->vz_SP   , parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->rho_SP  , parameter->getParD(lev)->rho_SP  , parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH(lev)->press_SP, parameter->getParD(lev)->press_SP, parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyDeviceToHost));
}
//sparse
void CudaMemoryManager::cudaAllocSP(int lev)
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->geoSP           ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborX_SP    ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborY_SP    ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->neighborZ_SP    ), parameter->getParH(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->rho_SP          ), parameter->getParH(lev)->mem_size_doubflo_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vx_SP           ), parameter->getParH(lev)->mem_size_doubflo_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vy_SP           ), parameter->getParH(lev)->mem_size_doubflo_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->vz_SP           ), parameter->getParH(lev)->mem_size_doubflo_SP));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->press_SP        ), parameter->getParH(lev)->mem_size_doubflo_SP));
	//Device						 
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->geoSP               ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborX_SP        ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborY_SP        ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->neighborZ_SP        ), parameter->getParD(lev)->mem_size_int_SP    ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->rho_SP              ), parameter->getParD(lev)->mem_size_doubflo_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vx_SP               ), parameter->getParD(lev)->mem_size_doubflo_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vy_SP               ), parameter->getParD(lev)->mem_size_doubflo_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->vz_SP               ), parameter->getParD(lev)->mem_size_doubflo_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->press_SP            ), parameter->getParD(lev)->mem_size_doubflo_SP));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->d0SP.f[0]           ), (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParD(lev)->mem_size_doubflo_SP));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 4. * (double)parameter->getParH(lev)->mem_size_int_SP + 5. * (double)parameter->getParH(lev)->mem_size_doubflo_SP + (double)parameter->getD3Qxx() * (double)parameter->getParH(lev)->mem_size_doubflo_SP;
	parameter->setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopySP(int lev)
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->geoSP       ,  parameter->getParH(lev)->geoSP       ,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborX_SP,  parameter->getParH(lev)->neighborX_SP,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborY_SP,  parameter->getParH(lev)->neighborY_SP,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->neighborZ_SP,  parameter->getParH(lev)->neighborZ_SP,  parameter->getParH(lev)->mem_size_int_SP     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->rho_SP      ,  parameter->getParH(lev)->rho_SP      ,  parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vx_SP       ,  parameter->getParH(lev)->vx_SP       ,  parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vy_SP       ,  parameter->getParH(lev)->vy_SP       ,  parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->vz_SP       ,  parameter->getParH(lev)->vz_SP       ,  parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD(lev)->press_SP    ,  parameter->getParH(lev)->press_SP    ,  parameter->getParH(lev)->mem_size_doubflo_SP , cudaMemcpyHostToDevice));
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
//Velo
void CudaMemoryManager::cudaAllocVeloBC(int lev)
{
	unsigned int mem_size_inflow_Q_k = sizeof(int)*parameter->getParH(lev)->Qinflow.kQ;
	unsigned int mem_size_inflow_Q_q = sizeof(doubflo)*parameter->getParH(lev)->Qinflow.kQ;

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
    parameter->setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyVeloBC(int lev)
{
	unsigned int mem_size_inflow_Q_k = sizeof(int)*parameter->getParH(lev)->Qinflow.kQ;
	unsigned int mem_size_inflow_Q_q = sizeof(doubflo)*parameter->getParH(lev)->Qinflow.kQ;

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
	unsigned int mem_size_outflow_Q_q = sizeof(doubflo)*parameter->getParH(lev)->Qoutflow.kQ;

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
	parameter->setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyOutflowBC(int lev)
{
	unsigned int mem_size_outflow_Q_k = sizeof(int)*parameter->getParH(lev)->Qoutflow.kQ;
	unsigned int mem_size_outflow_Q_q = sizeof(doubflo)*parameter->getParH(lev)->Qoutflow.kQ;

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
	unsigned int mem_size_Q_q      = sizeof(doubflo)*parameter->getParH(lev)->QWall.kQ;
	unsigned int mem_size_Q_value  = sizeof(long long)*parameter->getParH(lev)->QWall.kQ; //Geller
	unsigned int mem_size_Q_q_read = sizeof(doubflo)*parameter->getParH(lev)->kQread;     //Geller

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
	parameter->setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyWallBC(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QWall.kQ;
	unsigned int mem_size_Q_q = sizeof(doubflo)*parameter->getParH(lev)->QWall.kQ;

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
	unsigned int mem_size_Q_q      = sizeof(doubflo)*parameter->getParH(lev)->QGeom.kQ;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeom.q27[0]), parameter->getD3Qxx()*mem_size_Q_q      ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH(lev)->QGeom.k),                 mem_size_Q_k      ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeom.q27[0]),     parameter->getD3Qxx()* mem_size_Q_q     ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD(lev)->QGeom.k),                      mem_size_Q_k     ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_Q_k + (double)parameter->getD3Qxx()*(double)mem_size_Q_q;
	parameter->setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyGeomBC(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QGeom.kQ;
	unsigned int mem_size_Q_q = sizeof(doubflo)*parameter->getParH(lev)->QGeom.kQ;

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
	unsigned int mem_size_Q_q      = sizeof(doubflo)*parameter->getParH(lev)->QPress.kQ;

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
	parameter->setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyPress(int lev)
{
	unsigned int mem_size_Q_k = sizeof(int)*parameter->getParH(lev)->QPress.kQ;
	unsigned int mem_size_Q_q = sizeof(doubflo)*parameter->getParH(lev)->QPress.kQ;

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
	unsigned int mem_size = sizeof(doubflo) * 3;
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->forcingH), mem_size));
    parameter->forcingH[0] = parameter->getForcesDouble()[0];
    parameter->forcingH[1] = parameter->getForcesDouble()[1];
    parameter->forcingH[2] = parameter->getForcesDouble()[2];
	//Device
	checkCudaErrors( cudaMalloc((void**) &parameter->forcingD, mem_size));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (doubflo)mem_size;
	parameter->setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyForcingToDevice()
{
	unsigned int mem_size = sizeof(doubflo) * 3;
	checkCudaErrors( cudaMemcpy(parameter->forcingD, parameter->forcingH, mem_size, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyForcingToHost()
{
	unsigned int mem_size = sizeof(doubflo) * 3;
	checkCudaErrors( cudaMemcpy(parameter->forcingH, parameter->forcingD, mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeForcing()
{
	checkCudaErrors( cudaFreeHost(parameter->getForcesHost()));
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
	parameter->setMemsizeGPU(tmp, false);
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
	parameter->setMemsizeGPU(tmp, false);
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
	parameter->setMemsizeGPU(tmp, false);
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

std::shared_ptr<CudaMemoryManager> CudaMemoryManager::make(std::shared_ptr<Parameter> parameter)
{
    return std::shared_ptr<CudaMemoryManager>(new CudaMemoryManager(parameter));
}

CudaMemoryManager::CudaMemoryManager(std::shared_ptr<Parameter> parameter)
{
    this->parameter = parameter;
}

CudaMemoryManager::CudaMemoryManager(const CudaMemoryManager&)
{

}

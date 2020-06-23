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
//! \file CudaMemoryManager.cpp
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "CudaMemoryManager.h"
#include <Parameter/Parameter.h>

//coordinates
void CudaMemoryManager::cudaAllocCoord()
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->coordinateX ), parameter->getParH()->mem_size_real  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->coordinateY ), parameter->getParH()->mem_size_real  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->coordinateZ ), parameter->getParH()->mem_size_real  ));
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->coordinateX     ), parameter->getParH()->mem_size_real  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->coordinateY     ), parameter->getParH()->mem_size_real  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->coordinateZ     ), parameter->getParH()->mem_size_real  ));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 3. * (double)parameter->getParH()->mem_size_real;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyCoord()
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parameter->getParD()->coordinateX,  parameter->getParH()->coordinateX,  parameter->getParH()->mem_size_real     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->coordinateY,  parameter->getParH()->coordinateY,  parameter->getParH()->mem_size_real     , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->coordinateZ,  parameter->getParH()->coordinateZ,  parameter->getParH()->mem_size_real     , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeCoord()
{
	checkCudaErrors( cudaFreeHost(parameter->getParH()->coordinateX   ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->coordinateY   ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->coordinateZ   ));
}





//print
void CudaMemoryManager::cudaCopyDataToHost()
{
	checkCudaErrors( cudaMemcpy(parameter->getParH()->velocityX   , parameter->getParD()->velocityX   , parameter->getParH()->mem_size_real , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH()->velocityY   , parameter->getParD()->velocityY   , parameter->getParH()->mem_size_real , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH()->velocityZ   , parameter->getParD()->velocityZ   , parameter->getParH()->mem_size_real , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH()->rho         , parameter->getParD()->rho         , parameter->getParH()->mem_size_real , cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(parameter->getParH()->pressure    , parameter->getParD()->pressure    , parameter->getParH()->mem_size_real , cudaMemcpyDeviceToHost));
}





//sparse
void CudaMemoryManager::cudaAllocSP()
{
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->typeOfGridNode ), parameter->getParH()->mem_size_int ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->neighborX      ), parameter->getParH()->mem_size_int ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->neighborY      ), parameter->getParH()->mem_size_int ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->neighborZ      ), parameter->getParH()->mem_size_int ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->rho            ), parameter->getParH()->mem_size_real));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->velocityX      ), parameter->getParH()->mem_size_real));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->velocityY      ), parameter->getParH()->mem_size_real));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->velocityZ      ), parameter->getParH()->mem_size_real));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->pressure       ), parameter->getParH()->mem_size_real));
	//Device						 
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->typeOfGridNode     ), parameter->getParD()->mem_size_int ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->neighborX          ), parameter->getParD()->mem_size_int ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->neighborY          ), parameter->getParD()->mem_size_int ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->neighborZ          ), parameter->getParD()->mem_size_int ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->rho                ), parameter->getParD()->mem_size_real));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->velocityX          ), parameter->getParD()->mem_size_real));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->velocityY          ), parameter->getParD()->mem_size_real));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->velocityZ          ), parameter->getParD()->mem_size_real));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->pressure           ), parameter->getParD()->mem_size_real));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->distributions.f[0] ), (unsigned long long)parameter->getD3Qxx()*(unsigned long long)parameter->getParD()->mem_size_real));
	//////////////////////////////////////////////////////////////////////////
	double tmp = 4. * (double)parameter->getParH()->mem_size_int + 5. * (double)parameter->getParH()->mem_size_real + (double)parameter->getD3Qxx() * (double)parameter->getParH()->mem_size_real;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopySP()
{
	//copy host to device
	checkCudaErrors( cudaMemcpy(parameter->getParD()->typeOfGridNode ,  parameter->getParH()->typeOfGridNode ,  parameter->getParH()->mem_size_int  , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->neighborX      ,  parameter->getParH()->neighborX      ,  parameter->getParH()->mem_size_int  , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->neighborY      ,  parameter->getParH()->neighborY      ,  parameter->getParH()->mem_size_int  , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->neighborZ      ,  parameter->getParH()->neighborZ      ,  parameter->getParH()->mem_size_int  , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->rho            ,  parameter->getParH()->rho            ,  parameter->getParH()->mem_size_real , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->velocityX      ,  parameter->getParH()->velocityX      ,  parameter->getParH()->mem_size_real , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->velocityY      ,  parameter->getParH()->velocityY      ,  parameter->getParH()->mem_size_real , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->velocityZ      ,  parameter->getParH()->velocityZ      ,  parameter->getParH()->mem_size_real , cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->pressure       ,  parameter->getParH()->pressure       ,  parameter->getParH()->mem_size_real , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeSP()
{
	checkCudaErrors( cudaFreeHost(parameter->getParH()->typeOfGridNode  ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->velocityX       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->velocityY       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->velocityZ       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->rho             ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->pressure        ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->neighborX       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->neighborY       ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->neighborZ       ));
}





//velocity boundary condition
void CudaMemoryManager::cudaAllocVeloBC()
{
	unsigned int mem_size_inflow_BC_INT = sizeof(int)*parameter->getParH()->numberOfInflowBCnodes;
	unsigned int mem_size_inflow_BC_REAL = sizeof(real)*parameter->getParH()->numberOfInflowBCnodes;

	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->inflowBC.q27[0]), parameter->getD3Qxx() * mem_size_inflow_BC_REAL ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->inflowBC.k),								 mem_size_inflow_BC_INT  ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->inflowBC.Vx),				             mem_size_inflow_BC_REAL ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->inflowBC.Vy),				             mem_size_inflow_BC_REAL ));
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->inflowBC.Vz),				             mem_size_inflow_BC_REAL ));

	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->inflowBC.q27[0]),     parameter->getD3Qxx() * mem_size_inflow_BC_REAL ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->inflowBC.k),									 mem_size_inflow_BC_INT  ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->inflowBC.Vx),				                 mem_size_inflow_BC_REAL ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->inflowBC.Vy),				                 mem_size_inflow_BC_REAL ));
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->inflowBC.Vz),				                 mem_size_inflow_BC_REAL ));

	//////////////////////////////////////////////////////////////////////////
	double tmp = (double)mem_size_inflow_BC_INT + 4. * (double)mem_size_inflow_BC_REAL + (double)parameter->getD3Qxx() * (double)mem_size_inflow_BC_REAL;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyVeloBC()
{
	unsigned int mem_size_inflow_BC_INT = sizeof(int)*parameter->getParH()->numberOfInflowBCnodes;
	unsigned int mem_size_inflow_BC_REAL = sizeof(real)*parameter->getParH()->numberOfInflowBCnodes;

	checkCudaErrors( cudaMemcpy(parameter->getParD()->inflowBC.q27[0],  parameter->getParH()->inflowBC.q27[0], parameter->getD3Qxx() *	mem_size_inflow_BC_REAL,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->inflowBC.k,       parameter->getParH()->inflowBC.k,							    mem_size_inflow_BC_INT ,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->inflowBC.Vx,      parameter->getParH()->inflowBC.Vx,								mem_size_inflow_BC_REAL,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->inflowBC.Vy,      parameter->getParH()->inflowBC.Vy,								mem_size_inflow_BC_REAL,  cudaMemcpyHostToDevice));
	checkCudaErrors( cudaMemcpy(parameter->getParD()->inflowBC.Vz,      parameter->getParH()->inflowBC.Vz,								mem_size_inflow_BC_REAL,  cudaMemcpyHostToDevice));

}
void CudaMemoryManager::cudaFreeVeloBC()
{
	checkCudaErrors( cudaFreeHost(parameter->getParH()->inflowBC.q27[0] ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->inflowBC.k      ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->inflowBC.Vx     ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->inflowBC.Vy     ));
	checkCudaErrors( cudaFreeHost(parameter->getParH()->inflowBC.Vz     ));
}





//forcing
void CudaMemoryManager::cudaAllocForcing()
{
	unsigned int mem_size = sizeof(real) * 3;
	//Host
	checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->forcing), mem_size));
	parameter->getParH()->forcing[0] = (real)0.0;
	parameter->getParH()->forcing[1] = (real)0.0;
	parameter->getParH()->forcing[2] = (real)0.0;
	//Device
	checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->forcing), mem_size));
	//////////////////////////////////////////////////////////////////////////
	double tmp = (real)mem_size;
	setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyForcingToDevice()
{
	unsigned int mem_size = sizeof(real) * 3;
	checkCudaErrors( cudaMemcpy(parameter->getParD()->forcing, parameter->getParH()->forcing, mem_size, cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaCopyForcingToHost()
{
	unsigned int mem_size = sizeof(real) * 3;
	checkCudaErrors( cudaMemcpy(parameter->getParH()->forcing, parameter->getParD()->forcing, mem_size, cudaMemcpyDeviceToHost));
}
void CudaMemoryManager::cudaFreeForcing()
{
	checkCudaErrors( cudaFreeHost(parameter->getParH()->forcing));
}





//Neighbor WSB
void CudaMemoryManager::cudaAllocNeighborWSB()
{
    //Host
    checkCudaErrors( cudaMallocHost((void**) &(parameter->getParH()->neighborInverse    ), parameter->getParH()->mem_size_int    ));
    //Device
    checkCudaErrors( cudaMalloc((void**) &(parameter->getParD()->neighborInverse        ), parameter->getParD()->mem_size_int    ));
    //////////////////////////////////////////////////////////////////////////
    double tmp = (double)parameter->getParH()->mem_size_int;
    setMemsizeGPU(tmp, false);
}
void CudaMemoryManager::cudaCopyNeighborWSB()
{
    //copy host to device
    checkCudaErrors( cudaMemcpy(parameter->getParD()->neighborInverse,  parameter->getParH()->neighborInverse,  parameter->getParH()->mem_size_int     , cudaMemcpyHostToDevice));
}
void CudaMemoryManager::cudaFreeNeighborWSB()
{
    checkCudaErrors( cudaFreeHost(parameter->getParH()->neighborInverse));
}





std::shared_ptr<CudaMemoryManager> CudaMemoryManager::make(SPtr<Parameter> parameter)
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

CudaMemoryManager::CudaMemoryManager(SPtr<Parameter> parameter)
{
    this->parameter = parameter;
}

CudaMemoryManager::CudaMemoryManager(const CudaMemoryManager&)
{

}

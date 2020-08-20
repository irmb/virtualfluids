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
//! \file InitializerKernel.cu
//! \ingroup Initializer
//! \author Stephan Lenz
//=======================================================================================
#include "Initializer.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "CudaUtility/CudaRunKernel.hpp"

//! \brief This is a CUDA Kernel that computes the cell index and calls \ref initializeDataUpdateFunction for this index
__global__                 void initializeDataUpdateKernel  ( DataBaseStruct dataBase, uint numberOfEntities );

//! \brief Sets \ref DataBase::dataUpdate for one cell to zero 
__host__ __device__ inline void initializeDataUpdateFunction( DataBaseStruct dataBase, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Initializer::initializeDataUpdate( SPtr<DataBase> dataBase )
{
    CudaUtility::CudaGrid grid( dataBase->numberOfCells, 32 );

    runKernel( initializeDataUpdateKernel,
               initializeDataUpdateFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct() );

    cudaDeviceSynchronize();

    getLastCudaError("Initializer::initializeDataUpdate( SPtr<DataBase> dataBase )");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void initializeDataUpdateKernel(DataBaseStruct dataBase, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    initializeDataUpdateFunction( dataBase, index );
}

__host__ __device__ inline void initializeDataUpdateFunction(DataBaseStruct dataBase, uint index)
{
    dataBase.dataUpdate[ RHO__(index, dataBase.numberOfCells) ] = c0o1;
    dataBase.dataUpdate[ RHO_U(index, dataBase.numberOfCells) ] = c0o1;
    dataBase.dataUpdate[ RHO_V(index, dataBase.numberOfCells) ] = c0o1;
    dataBase.dataUpdate[ RHO_W(index, dataBase.numberOfCells) ] = c0o1;
    dataBase.dataUpdate[ RHO_E(index, dataBase.numberOfCells) ] = c0o1;
#ifdef USE_PASSIVE_SCALAR
	dataBase.dataUpdate[ RHO_S_1(index, dataBase.numberOfCells) ] = c0o1;
	dataBase.dataUpdate[ RHO_S_2(index, dataBase.numberOfCells) ] = c0o1;
#endif // USE_PASSIVE_SCALAR

    dataBase.massFlux[ VEC_X(index, dataBase.numberOfCells) ]   = c0o1;
    dataBase.massFlux[ VEC_Y(index, dataBase.numberOfCells) ]   = c0o1;
    dataBase.massFlux[ VEC_Z(index, dataBase.numberOfCells) ]   = c0o1;

    dataBase.diffusivity[ index ] = c1o1;
}

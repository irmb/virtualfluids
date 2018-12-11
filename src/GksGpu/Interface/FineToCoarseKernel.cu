#include "Interface.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"

#include "DataBase/DataBaseStruct.h"

#include "FlowStateData/FlowStateData.cuh"

#include "Definitions/PassiveScalar.h"
#include "Definitions/MemoryAccessPattern.h"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void fineToCoarseKernel  ( DataBaseStruct dataBase, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void fineToCoarseFunction( DataBaseStruct dataBase, uint startIndex, uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void Interface::runFineToCoarse( SPtr<DataBase> dataBase, uint level )
{    
    CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFineToCoarse, 128);

    runKernel(fineToCoarseKernel,
              fineToCoarseFunction,
              dataBase->getDeviceType(), grid,
              dataBase->toStruct(),
              dataBase->perLevelCount[level].startOfFineToCoarse);

    cudaDeviceSynchronize();

    getLastCudaError("Interface::runFineToCoarse( SPtr<DataBase> dataBase, uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void fineToCoarseKernel( DataBaseStruct dataBase, uint startIndex, uint numberOfEntities )
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    fineToCoarseFunction( dataBase, startIndex, index );
}

__host__ __device__ inline void fineToCoarseFunction( DataBaseStruct dataBase, uint startIndex, uint index )
{
    index += startIndex;

    ConservedVariables cons;

#pragma unroll
    for( uint childIdx = 1; childIdx < LENGTH_FINE_TO_COARSE; childIdx++ ){

        uint cellIdx = dataBase.fineToCoarse[ FINE_TO_COARSE( index, childIdx, dataBase.numberOfCoarseGhostCells ) ];

        cons.rho  += c1o8 * dataBase.data[ RHO__(cellIdx, dataBase.numberOfCells) ];
        cons.rhoU += c1o8 * dataBase.data[ RHO_U(cellIdx, dataBase.numberOfCells) ];
        cons.rhoV += c1o8 * dataBase.data[ RHO_V(cellIdx, dataBase.numberOfCells) ];
        cons.rhoW += c1o8 * dataBase.data[ RHO_W(cellIdx, dataBase.numberOfCells) ];
        cons.rhoE += c1o8 * dataBase.data[ RHO_E(cellIdx, dataBase.numberOfCells) ];
    #ifdef USE_PASSIVE_SCALAR
	    cons.rhoS += c1o8 * dataBase.data[ RHO_S(cellIdx, dataBase.numberOfCells) ];
    #endif // USE_PASSIVE_SCALAR
    }

    uint cellIdx = dataBase.fineToCoarse[ FINE_TO_COARSE( index, 0, dataBase.numberOfCoarseGhostCells ) ];

    dataBase.data[ RHO__(cellIdx, dataBase.numberOfCells) ] = cons.rho ;
    dataBase.data[ RHO_U(cellIdx, dataBase.numberOfCells) ] = cons.rhoU;
    dataBase.data[ RHO_V(cellIdx, dataBase.numberOfCells) ] = cons.rhoV;
    dataBase.data[ RHO_W(cellIdx, dataBase.numberOfCells) ] = cons.rhoW;
    dataBase.data[ RHO_E(cellIdx, dataBase.numberOfCells) ] = cons.rhoE;
#ifdef USE_PASSIVE_SCALAR
	dataBase.data[ RHO_S(cellIdx, dataBase.numberOfCells) ] = cons.rhoS;
#endif // USE_PASSIVE_SCALAR
}

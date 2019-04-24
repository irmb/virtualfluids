#include "Interface.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"

#include "DataBase/DataBaseStruct.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"

#include "FlowStateData/AccessDeviceData.cuh"

#include "Definitions/PassiveScalar.h"
#include "Definitions/MemoryAccessPattern.h"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void fineToCoarseKernel  ( DataBaseStruct dataBase, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void fineToCoarseFunction                      ( DataBaseStruct dataBase, uint startIndex, uint index );
__host__ __device__ inline void fineToCoarseFunctionPrimitiveInterpolation( DataBaseStruct dataBase, uint startIndex, uint index );

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
    //fineToCoarseFunctionPrimitiveInterpolation( dataBase, startIndex, index );
}

__host__ __device__ inline void fineToCoarseFunction( DataBaseStruct dataBase, uint startIndex, uint index )
{
    index += startIndex;

    ConservedVariables parentCons;

#pragma unroll
    for( uint childIdx = 1; childIdx < LENGTH_FINE_TO_COARSE; childIdx++ ){

        uint cellIdx = dataBase.fineToCoarse[ FINE_TO_COARSE( index, childIdx, dataBase.numberOfCoarseGhostCells ) ];

        ConservedVariables cons;

        readCellData( cellIdx, dataBase, cons );

        parentCons = parentCons + c1o8 * cons;
    }

    {
        uint cellIdx = dataBase.fineToCoarse[FINE_TO_COARSE(index, 0, dataBase.numberOfCoarseGhostCells)];

        writeCellData(cellIdx, dataBase, parentCons);
    }
}

__host__ __device__ inline void fineToCoarseFunctionPrimitiveInterpolation( DataBaseStruct dataBase, uint startIndex, uint index )
{
    index += startIndex;

    PrimitiveVariables parentPrim;

#pragma unroll
    for( uint childIdx = 1; childIdx < LENGTH_FINE_TO_COARSE; childIdx++ ){

        uint cellIdx = dataBase.fineToCoarse[ FINE_TO_COARSE( index, childIdx, dataBase.numberOfCoarseGhostCells ) ];

        ConservedVariables cons;

        readCellData( cellIdx, dataBase, cons );

        parentPrim = parentPrim + c1o8 * toPrimitiveVariables(cons, two);
    }

    {
        uint cellIdx = dataBase.fineToCoarse[FINE_TO_COARSE(index, 0, dataBase.numberOfCoarseGhostCells)];

        ConservedVariables parentCons = toConservedVariables(parentPrim, two);

        writeCellData(cellIdx, dataBase, parentCons);
    }
}

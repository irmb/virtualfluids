#include "FluxComputation.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

__global__                 void fluxKernel  ( DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void fluxFunction( DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, uint level )
{
    {
        CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFacesX, 32);

        runKernel(fluxKernel,
            fluxFunction,
            dataBase->getDeviceType(), grid,
            dataBase->toStruct(),
            parameters,
            'x',
            dataBase->perLevelCount[level].startOfFacesX);

        getLastCudaError("FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, 'x', uint level )");
    }
    {
        CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFacesY, 32);

        runKernel(fluxKernel,
            fluxFunction,
            dataBase->getDeviceType(), grid,
            dataBase->toStruct(),
            parameters,
            'y',
            dataBase->perLevelCount[level].startOfFacesY);

        getLastCudaError("FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, 'y', uint level )");
    }
    {
        CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFacesZ, 32);

        runKernel(fluxKernel,
            fluxFunction,
            dataBase->getDeviceType(), grid,
            dataBase->toStruct(),
            parameters,
            'z',
            dataBase->perLevelCount[level].startOfFacesZ);

        getLastCudaError("FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, 'z', uint level )");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void fluxKernel(DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index > numberOfEntities ) return;

    fluxFunction( dataBase, parameters, direction, startIndex, index );
}

__host__ __device__ inline void fluxFunction(DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint index)
{
    uint faceIndex = startIndex + index;

    //////////////////////////////////////////////////////////////////////////



}

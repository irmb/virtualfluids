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

__global__                 void initializeDataUpdateKernel  ( DataBaseStruct dataBase, uint numberOfEntities );

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

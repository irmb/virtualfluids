#include "Communicator.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseAllocator.h"
#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void sendBufferKernel  ( const DataBaseStruct dataBase, 
                                                    const uint numberOfSendNodes,
                                                    const uint* sendIndices,
                                                    real* sendBuffer,
                                                    const uint startIndex,
                                                    const uint numberOfEntities );

__host__ __device__ inline void sendBufferFunction( const DataBaseStruct dataBase, 
                                                    const uint numberOfSendNodes,
                                                    const uint* sendIndices,
                                                    real* sendBuffer,
                                                    const uint startIndex,
                                                    const uint index );

__global__                 void recvBufferKernel  ( const DataBaseStruct dataBase, 
                                                    const uint numberOfRecvNodes,
                                                    const uint* recvIndices,
                                                    real* recvBuffer,
                                                    const uint startIndex,
                                                    const uint numberOfEntities );

__host__ __device__ inline void recvBufferFunction( const DataBaseStruct dataBase, 
                                                    const uint numberOfRecvNodes,
                                                    const uint* recvIndices,
                                                    real* recvBuffer,
                                                    const uint startIndex,
                                                    const uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Communicator::copyFromMeshToSendBuffer(const SPtr<DataBase> dataBase)
{    
    CudaUtility::CudaGrid grid( this->numberOfSendNodes, 32 );

    runKernel( sendBufferKernel,
               sendBufferFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               this->numberOfSendNodes,
               this->sendIndices,
               this->sendBuffer,
               0 );

    cudaDeviceSynchronize();

    getLastCudaError("Communicator::copyFromMeshToSendBuffer(const SPtr<DataBase> dataBase)");
}

void Communicator::copyFromRecvBufferToMesh(const SPtr<DataBase> dataBase)
{    
    CudaUtility::CudaGrid grid( this->numberOfRecvNodes, 32 );

    runKernel( recvBufferKernel,
               recvBufferFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               this->numberOfRecvNodes,
               this->recvIndices,
               this->recvBuffer,
               0 );

    cudaDeviceSynchronize();

    getLastCudaError("Communicator::copyFromRecvBufferToMesh(const SPtr<DataBase> dataBase)");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void sendBufferKernel( const DataBaseStruct dataBase, 
                                  const uint numberOfSendNodes,
                                  const uint* sendIndices,
                                  real* sendBuffer,
                                  const uint startIndex,
                                  const uint numberOfEntities )
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    sendBufferFunction( dataBase, numberOfSendNodes, sendIndices, sendBuffer, startIndex, index );
}

__host__ __device__ inline void sendBufferFunction( const DataBaseStruct dataBase, 
                                                    const uint numberOfSendNodes,
                                                    const uint* sendIndices,
                                                    real* sendBuffer,
                                                    const uint startIndex,
                                                    const uint index )
{
    

    uint cellIdx  = sendIndices [ index ];

    sendBuffer[ RHO__(index, numberOfSendNodes) ] = dataBase.data[ RHO__(cellIdx, dataBase.numberOfCells) ];
    sendBuffer[ RHO_U(index, numberOfSendNodes) ] = dataBase.data[ RHO_U(cellIdx, dataBase.numberOfCells) ];
    sendBuffer[ RHO_V(index, numberOfSendNodes) ] = dataBase.data[ RHO_V(cellIdx, dataBase.numberOfCells) ];
    sendBuffer[ RHO_W(index, numberOfSendNodes) ] = dataBase.data[ RHO_W(cellIdx, dataBase.numberOfCells) ];
    sendBuffer[ RHO_E(index, numberOfSendNodes) ] = dataBase.data[ RHO_E(cellIdx, dataBase.numberOfCells) ];
#ifdef USE_PASSIVE_SCALAR
    sendBuffer[ RHO_S_1(index, numberOfSendNodes) ] = dataBase.data[ RHO_S_1(cellIdx, dataBase.numberOfCells) ];
    sendBuffer[ RHO_S_2(index, numberOfSendNodes) ] = dataBase.data[ RHO_S_2(cellIdx, dataBase.numberOfCells) ];
#endif // USE_PASSIVE_SCALAR
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void recvBufferKernel( const DataBaseStruct dataBase, 
                                  const uint numberOfRecvNodes,
                                  const uint* recvIndices,
                                  real* recvBuffer,
                                  const uint startIndex,
                                  const uint numberOfEntities )
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    recvBufferFunction( dataBase, numberOfRecvNodes, recvIndices, recvBuffer, startIndex, index );
}

__host__ __device__ inline void recvBufferFunction( const DataBaseStruct dataBase, 
                                                    const uint numberOfRecvNodes,
                                                    const uint* recvIndices,
                                                    real* recvBuffer,
                                                    const uint startIndex,
                                                    const uint index )
{
    

    uint cellIdx  = recvIndices [ index ];

    dataBase.data[ RHO__(cellIdx, dataBase.numberOfCells) ] = recvBuffer[ RHO__(index, numberOfRecvNodes) ] ;
    dataBase.data[ RHO_U(cellIdx, dataBase.numberOfCells) ] = recvBuffer[ RHO_U(index, numberOfRecvNodes) ] ;
    dataBase.data[ RHO_V(cellIdx, dataBase.numberOfCells) ] = recvBuffer[ RHO_V(index, numberOfRecvNodes) ] ;
    dataBase.data[ RHO_W(cellIdx, dataBase.numberOfCells) ] = recvBuffer[ RHO_W(index, numberOfRecvNodes) ] ;
    dataBase.data[ RHO_E(cellIdx, dataBase.numberOfCells) ] = recvBuffer[ RHO_E(index, numberOfRecvNodes) ] ;
#ifdef USE_PASSIVE_SCALAR
    dataBase.data[ RHO_S_1(cellIdx, dataBase.numberOfCells) ] = recvBuffer[ RHO_S_1(index, numberOfRecvNodes) ] ;
    dataBase.data[ RHO_S_2(cellIdx, dataBase.numberOfCells) ] = recvBuffer[ RHO_S_2(index, numberOfRecvNodes) ] ;
#endif // USE_PASSIVE_SCALAR
}
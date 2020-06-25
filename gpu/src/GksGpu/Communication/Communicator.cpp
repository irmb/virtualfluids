#include "Communicator.h"

#ifdef VF_DOUBLE_ACCURACY
#define MPI_TYPE_GPU  MPI_DOUBLE
#else
#define MPI_TYPE_GPU  MPI_FLOAT
#endif

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

#include "Core/PointerDefinitions.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseAllocator.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/CudaAwareMpi.h"

#include "CudaUtility/CudaUtility.h"

namespace GksGpu {

int Communicator::tagSendPositive = 0;
int Communicator::tagSendNegative = 1;

Communicator::Communicator( SPtr<DataBase> dataBase )
    : myAllocator ( dataBase->myAllocator )
{
    this->numberOfSendNodes = INVALID_INDEX;
    this->numberOfRecvNodes = INVALID_INDEX;

    this->sendIndices    = nullptr;
    this->recvIndices    = nullptr;
    this->sendBuffer     = nullptr;
    this->recvBuffer     = nullptr;
    this->sendBufferHost = nullptr;
    this->recvBufferHost = nullptr;

    this->sendBufferIsReady = MPI_REQUEST_NULL;
}

void Communicator::initialize(GksMeshAdapter & adapter, uint level, uint direction)
{
    this->myAllocator->freeMemory( *this );

    this->numberOfSendNodes = adapter.communicationIndices[level].sendIndices[direction].size();
    this->numberOfRecvNodes = adapter.communicationIndices[level].recvIndices[direction].size();

    this->myAllocator->allocateMemory( *this, adapter.communicationIndices[level].sendIndices[direction], 
                                              adapter.communicationIndices[level].recvIndices[direction] );

    this->opposingRank = adapter.communicationProcesses[direction];
}

void Communicator::sendData( SPtr<DataBase> dataBase, int tag )
{
#ifdef USE_CUDA_AWARE_MPI

    this->copyFromMeshToSendBuffer( dataBase );
    
    MPI_Isend( this->sendBuffer, this->numberOfSendNodes * LENGTH_CELL_DATA, MPI_TYPE_GPU, this->opposingRank, tag, MPI_COMM_WORLD, &this->sendBufferIsReady );

#else // USE_CUDA_AWARE_MPI

    this->copyFromMeshToSendBuffer( dataBase );

    MPI_Wait(&this->sendBufferIsReady, MPI_STATUSES_IGNORE);

    this->myAllocator->copyBuffersDeviceToHost( shared_from_this() );
    
    CudaUtility::synchronizeCudaStream( CudaUtility::communicationStream );

    MPI_Isend( this->sendBufferHost, this->numberOfSendNodes * LENGTH_CELL_DATA, MPI_TYPE_GPU, this->opposingRank, tag, MPI_COMM_WORLD, &this->sendBufferIsReady );

#endif // USE_CUDA_AWARE_MPI
}

void Communicator::recvData( SPtr<DataBase> dataBase, int tag )
{
#ifdef USE_CUDA_AWARE_MPI
    
    MPI_Recv ( this->recvBuffer, this->numberOfRecvNodes * LENGTH_CELL_DATA, MPI_TYPE_GPU, this->opposingRank, tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE );
    
    this->copyFromRecvBufferToMesh( dataBase );

#else // USE_CUDA_AWARE_MPI
    
    MPI_Recv ( this->recvBufferHost, this->numberOfRecvNodes * LENGTH_CELL_DATA, MPI_TYPE_GPU, this->opposingRank, tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE );
    
    this->myAllocator->copyBuffersHostToDevice( shared_from_this() );

    this->copyFromRecvBufferToMesh( dataBase );

#endif // USE_CUDA_AWARE_MPI
}

} // namespace GksGpu

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

Communicator::Communicator( SPtr<DataBase> dataBase )
    : myAllocator ( dataBase->myAllocator )
{
    this->numberOfSendNodes = INVALID_INDEX;
    this->numberOfRecvNodes = INVALID_INDEX;

    this->sendIndices = nullptr;
    this->recvIndices = nullptr;
    this->sendBuffer  = nullptr;
    this->recvBuffer  = nullptr;

    this->sendBufferIsFresh = false;
}

void Communicator::initialize(GksMeshAdapter & adapter, uint direction)
{
    this->myAllocator->freeMemory( *this );

    this->numberOfSendNodes = adapter.sendIndices[direction].size();
    this->numberOfRecvNodes = adapter.recvIndices[direction].size();

    this->myAllocator->allocateMemory( *this, adapter.sendIndices[direction], adapter.recvIndices[direction] );
}

void Communicator::exchangeData( SPtr<DataBase> dataBase )
{
    while( this->sendBufferIsFresh ) /* do nothing */;

    this->copyFromMeshToSendBuffer( dataBase );

    this->sendBufferIsFresh = true;

    //////////////////////////////////////////////////////////////////////////

    this->myAllocator->copyDataDeviceToDevice(shared_from_this(), opposingCommunicator);

    this->opposingCommunicator->sendBufferIsFresh = false;

    //////////////////////////////////////////////////////////////////////////

    this->copyFromRecvBufferToMesh( dataBase );
}

#ifndef Communicator_H
#define Communicator_H

#include <memory>
//#include <mutex>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

class  GksMeshAdapter;
class  DataBaseAllocator;
struct DataBase;

struct VF_PUBLIC Communicator : public std::enable_shared_from_this<Communicator>
{
    SPtr<DataBaseAllocator> myAllocator;

    uint numberOfSendNodes;
    uint numberOfRecvNodes;

    uint* sendIndices;
    uint* recvIndices;

    real* sendBuffer;
    real* recvBuffer;

    SPtr<Communicator> opposingCommunicator;

    bool sendBufferIsFresh;

    Communicator( SPtr<DataBase> dataBase );

    void initialize( GksMeshAdapter& adapter, uint direction );

    void exchangeData( SPtr<DataBase> dataBase );

    void copyFromMeshToSendBuffer( SPtr<DataBase> dataBase );

    void copyFromRecvBufferToMesh( SPtr<DataBase> dataBase );
};

#endif

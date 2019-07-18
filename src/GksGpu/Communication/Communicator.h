#ifndef Communicator_H
#define Communicator_H

#include <memory>
#include <vector>
#include <mpi.h>
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

    uint* sendIndices; // device
    uint* recvIndices; // device

    real* sendBuffer; // device
    real* recvBuffer; // device

    real* sendBufferHost; // pinned memory
    real* recvBufferHost; // pinned memory

    uint rank;
    uint opposingRank;

    MPI_Request sendBufferIsReady;

    static int tagSendPositive;
    static int tagSendNegative;

    //////////////////////////////////////////////////////////////////////////

    Communicator( SPtr<DataBase> dataBase );

    void initialize( GksMeshAdapter& adapter, uint level, uint direction );

    void copyFromMeshToSendBuffer( SPtr<DataBase> dataBase );

    void copyFromRecvBufferToMesh( SPtr<DataBase> dataBase );

    void sendData( SPtr<DataBase> dataBase, int tag );
    void recvData( SPtr<DataBase> dataBase, int tag );
};

#endif

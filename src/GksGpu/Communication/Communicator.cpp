#include "Communicator.h"

#ifdef zero
#undef zero
#endif

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <mutex>
#include <condition_variable>

#include "Core/PointerDefinitions.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseAllocator.h"

// from https://stackoverflow.com/questions/24465533/implementing-boostbarrier-in-c11
class Barrier
{
private:
    std::mutex m_mutex;
    std::condition_variable m_cv;

    size_t m_count;
    const size_t m_initial;

    enum State : unsigned char {
        Up, Down
    };
    State m_state;

public:
    explicit Barrier(std::size_t count) : m_count{ count }, m_initial{ count }, m_state{ State::Down } { }

    /// Blocks until all N threads reach here
    void Sync()
    {
        std::unique_lock<std::mutex> lock{ m_mutex };

        if (m_state == State::Down)
        {
            // Counting down the number of syncing threads
            if (--m_count == 0) {
                m_state = State::Up;
                m_cv.notify_all();
            }
            else {
                m_cv.wait(lock, [this] { return m_state == State::Up; });
            }
        }

        else // (m_state == State::Up)
        {
            // Counting back up for Auto reset
            if (++m_count == m_initial) {
                m_state = State::Down;
                m_cv.notify_all();
            }
            else {
                m_cv.wait(lock, [this] { return m_state == State::Down; });
            }
        }
    }
};  

Barrier barrier(2);

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
    while( this->sendBufferIsFresh );

    this->copyFromMeshToSendBuffer( dataBase );

    this->sendBufferIsFresh = true;

    barrier.Sync();

    //////////////////////////////////////////////////////////////////////////

    while( !this->opposingCommunicator->sendBufferIsFresh );

    this->myAllocator->copyDataDeviceToDevice(shared_from_this(), opposingCommunicator);

    this->opposingCommunicator->sendBufferIsFresh = false;

    //////////////////////////////////////////////////////////////////////////

    this->copyFromRecvBufferToMesh( dataBase );
}

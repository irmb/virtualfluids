#ifndef VF_GPU_COMMUNICATIONROUTINEMOCKS_H
#define VF_GPU_COMMUNICATIONROUTINEMOCKS_H

#include "CommunicationRoutine.h"

namespace vf::gpu::test 
{

class CommunicationRoutineTestDouble : public vf::gpu::CommunicationRoutine
{
public:
    void receive_send(uint *buffer_receive, int size_buffer_recv, int neighbor_rank_recv, uint *buffer_send,
                              int size_buffer_send, int neighbor_rank_send) const override { } 
    int getPID() const override { return 0; }
};

}



#endif

#ifndef INDEX_EXCHANGE
#define INDEX_EXCHANGE

#include <basics/Core/DataTypes.h>

namespace vf::gpu
{
class IndexExchange
{
public:
    virtual void exchangeIndices(uint *buffer_receive, int size_buffer_recv, int neighbor_rank_recv, uint *buffer_send,
                                 int size_buffer_send, int neighbor_rank_send) const = 0;
    virtual int getPID() const = 0;
};
} // namespace vf::gpu

#endif
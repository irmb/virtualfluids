#ifndef COMMUNICATOR_GPU_H
#define COMMUNICATOR_GPU_H

#include <vector>
#include <basics/DataTypes.h>

#include "VirtualFluids_GPU_export.h"
#include "CommunicationRoutine.h"

namespace vf::gpu
{

class VIRTUALFLUIDS_GPU_EXPORT Communicator : public CommunicationRoutine
{
public:
    virtual void waitAll() = 0;
    virtual int getPID() const override = 0;
    virtual int getNumberOfProcess() const = 0;
    virtual void exchngData(float *sbuf_t, float *rbuf_t, float *sbuf_b, float *rbuf_b, int count) = 0;
    //////////////////////////////////////////////////////////////////////////
    virtual void exchngDataGPU(real *sbuf, int count_s, real *rbuf, int count_r, int nb_rank) = 0;
    virtual void nbRecvDataGPU(real *rbuf, int count_r, int nb_rank) = 0;
    virtual void nbSendDataGPU(real *sbuf, int count_s, int nb_rank) = 0;
    virtual void waitallGPU() = 0;
    virtual void sendDataGPU(real *sbuf, int count_s, int nb_rank) = 0;
    virtual void waitGPU(int id) = 0;
    virtual void resetRequest() = 0;
    //////////////////////////////////////////////////////////////////////////
    virtual int mapCudaDevice(const int &rank, const int &size, const std::vector<unsigned int> &devices, const int &maxdev) = 0;
    virtual std::vector<double> gatherNUPS(double processNups) = 0;
    virtual double sumNups(double processNups) = 0;
    //////////////////////////////////////////////////////////////////////////
    virtual void receive_send(uint *buffer_receive, int size_buffer_recv, int neighbor_rank_recv, uint *buffer_send,
                              int size_buffer_send, int neighbor_rank_send) const override = 0;

};

} // namespace vf::gpu

#endif

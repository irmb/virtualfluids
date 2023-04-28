#ifndef MPIMpiCommunicator_GPU_H
#define MPIMpiCommunicator_GPU_H

#include <vector>

#include <mpi.h>

#include "VirtualFluids_GPU_export.h"

#include "Communicator.h"
#include <basics/DataTypes.h>

//////////////////////////////////
#ifdef VF_DOUBLE_ACCURACY
#define MPI_Type_GPU MPI_DOUBLE
#else
#define MPI_Type_GPU MPI_FLOAT
#endif
//////////////////////////////////

namespace vf::gpu
{

class VIRTUALFLUIDS_GPU_EXPORT MpiCommunicator : public Communicator
{
public:
    static MpiCommunicator &getInstance();
    MpiCommunicator(const MpiCommunicator &) = delete;
    MpiCommunicator &operator=(const MpiCommunicator &) = delete;
    ~MpiCommunicator() override;

    void exchngBottomToTop(float *sbuf, float *rbuf, int count);
    void exchngTopToBottom(float *sbuf, float *rbuf, int count);
    void waitAll() override;
    void distributeGeometry(unsigned int *dataRoot, unsigned int *dataNode, int dataSizePerNode);
    int getPID() const override;
    int getNumberOfProcess() const override;
    int getNeighbourTop();
    int getNeighbourBottom();
    void exchngData(float *sbuf_t, float *rbuf_t, float *sbuf_b, float *rbuf_b, int count) override;
    void exchngDataNB(float *sbuf_t, int count_st, float *rbuf_t, int count_rt, float *sbuf_b, int count_sb,
                      float *rbuf_b, int count_rb);
    //////////////////////////////////////////////////////////////////////////
    void exchngDataGPU(real *sbuf, int count_s, real *rbuf, int count_r, int nb_rank) override;
    void sendRecvGPU(real *sbuf, int count_s, real *rbuf, int count_r, int nb_rank);
    void nbRecvDataGPU(real *rbuf, int count_r, int nb_rank) override;
    void nbSendDataGPU(real *sbuf, int count_s, int nb_rank) override;
    void waitallGPU() override;
    void sendDataGPU(real *sbuf, int count_s, int nb_rank) override;
    void waitGPU(int id) override;
    void resetRequest() override;
    void barrierGPU();
    void barrier();
    //////////////////////////////////////////////////////////////////////////
    void exchngDataGeo(int *sbuf_t, int *rbuf_t, int *sbuf_b, int *rbuf_b, int count);
    MPI_Comm getMpiCommunicator();
    int mapCudaDevice(const int &rank, const int &size, const std::vector<unsigned int> &devices, const int &maxdev) override;
    double reduceSum(double quantityPerProcess) override;
    //////////////////////////////////////////////////////////////////////////
    void receive_send(uint *buffer_receive, int size_buffer_recv, int neighbor_rank_recv, uint *buffer_send,
                      int size_buffer_send, int neighbor_rank_send) const override;

private:
    int numprocs, PID;
    int nbrbottom, nbrtop;
    MPI_Comm comm1d, commGPU;
    MPI_Status status[4];
    MPI_Request request[4];
    //////////////////////////////////////////////////////////////////////////
    std::vector<MPI_Request> requestGPU;
    int rcount;
    //////////////////////////////////////////////////////////////////////////
    double starttime;
    double endtime;
    MpiCommunicator();
};

} // namespace vf::gpu

#endif

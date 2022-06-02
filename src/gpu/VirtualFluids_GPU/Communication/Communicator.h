#ifndef COMMUNICATOR_GPU_H
#define COMMUNICATOR_GPU_H

#include <vector>

#include <mpi.h>

#include "VirtualFluids_GPU_export.h"

#include <basics/Core/DataTypes.h>

//////////////////////////////////
#ifdef VF_DOUBLE_ACCURACY
#define MPI_Type_GPU  MPI_DOUBLE
#else
#define MPI_Type_GPU  MPI_FLOAT
#endif
//////////////////////////////////


namespace vf::gpu
{


class VIRTUALFLUIDS_GPU_EXPORT Communicator
{
public:
    static Communicator& getInstance();
    Communicator(const Communicator&) = delete;
    Communicator& operator=(const Communicator&) = delete;

    void exchngBottomToTop(float* sbuf, float* rbuf, int count);
    void exchngTopToBottom(float* sbuf, float* rbuf, int count);
    void waitAll();
    void distributeGeometry(unsigned int* dataRoot, unsigned int* dataNode, int dataSizePerNode);
    int getPID() const;
    int getNummberOfProcess() const;
    int getNeighbourTop();
    int getNeighbourBottom();
    void exchngData(float* sbuf_t, float* rbuf_t, float* sbuf_b, float* rbuf_b, int count);
    void exchngDataNB(float* sbuf_t, int count_st, float* rbuf_t, int count_rt, float* sbuf_b, int count_sb, float* rbuf_b, int count_rb);
    //////////////////////////////////////////////////////////////////////////
    void exchngDataGPU(real* sbuf, int count_s, real* rbuf, int count_r, int nb_rank);
    void sendRecvGPU(real* sbuf, int count_s, real* rbuf, int count_r, int nb_rank);
    void nbRecvDataGPU( real* rbuf, int count_r, int nb_rank );
    void nbSendDataGPU( real* sbuf, int count_s, int nb_rank );
    void waitallGPU();
    void sendDataGPU( real* sbuf, int count_s, int nb_rank );
    void waitGPU(int id);
    void resetRequest();
    void barrierGPU();
    void barrier();
    //////////////////////////////////////////////////////////////////////////
    void exchngDataGeo(int* sbuf_t, int* rbuf_t, int* sbuf_b, int* rbuf_b, int count);
    MPI_Comm getCommunicator();
    void startTimer();
    void stopTimer();
    double getTime();
    int mapCudaDevice(const int &rank, const int &size, const std::vector<unsigned int> &devices, const int &maxdev);
    std::vector<double> gatherNUPS(double processNups);
    double sumNups(double processNups);
    //////////////////////////////////////////////////////////////////////////
    void exchangeIndices(uint *rbuf, int count_r, int nb_rank_r, uint *sbuf, int count_s, int nb_rank_s);
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
   Communicator();
   ~Communicator();
};

}

#endif


#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <mpi.h>
#include <vector>
#include <memory>
#include "LBM/LB.h"

//////////////////////////////////
#ifdef ISDOUBLE
#define MPI_Type_GPU  MPI_DOUBLE
#endif
//////////////////////////////////
#ifdef ISFLOAT
#define MPI_Type_GPU  MPI_FLOAT
#endif
//////////////////////////////////



class Communicator
{
public:
    static std::shared_ptr<Communicator> getInstance();
    void exchngBottomToTop(float* sbuf, float* rbuf, int count);
    void exchngTopToBottom(float* sbuf, float* rbuf, int count);
    void waitAll();
    void distributeGeometry(unsigned int* dataRoot, unsigned int* dataNode, int dataSizePerNode);
    int getPID();
    int getNummberOfProcess();
    int getNeighbourTop();
    int getNeighbourBottom();
    void exchngData(float* sbuf_t, float* rbuf_t, float* sbuf_b, float* rbuf_b, int count);
    void exchngDataNB(float* sbuf_t, int count_st, float* rbuf_t, int count_rt, float* sbuf_b, int count_sb, float* rbuf_b, int count_rb);
    //////////////////////////////////////////////////////////////////////////
    void exchngDataGPU(doubflo* sbuf, int count_s, doubflo* rbuf, int count_r, int nb_rank);
    void sendRecvGPU(doubflo* sbuf, int count_s, doubflo* rbuf, int count_r, int nb_rank);
    void nbRecvDataGPU(doubflo* rbuf, int count_r, int nb_rank);
    void nbSendDataGPU(doubflo* sbuf, int count_s, int nb_rank);
    void waitallGPU();
    void sendDataGPU(doubflo* sbuf, int count_s, int nb_rank);
    void waitGPU(int id);
    void resetRequest();
    void barrierGPU();
    //////////////////////////////////////////////////////////////////////////
    void exchngDataGeo(int* sbuf_t, int* rbuf_t, int* sbuf_b, int* rbuf_b, int count);
    MPI_Comm getCommunicator();
    void startTimer();
    void stopTimer();
    double getTime();
    int mapCudaDevice(const int &rank, const int &size, const std::vector<int> &devices, const int &maxdev);
protected:
private:
    static std::shared_ptr<Communicator> instance;
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
    Communicator(const int numberOfProcs);
    Communicator(const Communicator&);
};

#endif


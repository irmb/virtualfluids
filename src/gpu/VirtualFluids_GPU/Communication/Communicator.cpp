#include "Communicator.h"

#include <mpi.h>
#include <vector>

#include <logger/Logger.h>

#if defined (_WIN32) || defined (_WIN64)
   #include <Winsock2.h>
#elif defined (__unix__)
   #include <unistd.h>
#endif
//lib for windows Ws2_32.lib

namespace vf::gpu
{


Communicator::Communicator()
{
    int mpiInitialized = 0; // false
    MPI_Initialized(&mpiInitialized);
    if (!mpiInitialized) {
        MPI_Init(NULL, NULL);
        VF_LOG_TRACE("vf::gpu::Communicator(): MPI_Init");
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &PID);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    // Get a new communicator for a decomposition of the domain
    int isperiodic[1] = { 0 };
    MPI_Cart_create(MPI_COMM_WORLD, 1, &numprocs, isperiodic, 1, &comm1d);

    // Get my position in this communicator, and my neighbors
    MPI_Cart_shift(comm1d, 0, 1, &nbrbottom, &nbrtop);
}

Communicator::~Communicator()
{
    // proof if MPI is finalized
    int _mpiFinalized = 0; // false
    MPI_Finalized(&_mpiFinalized);
    if (!_mpiFinalized) {
        MPI_Finalize();
        VF_LOG_TRACE("vf::gpu::~Communicator(): MPI_Finalize");
    }
}


// C++11 thread safe singelton implementation:
// https://stackoverflow.com/questions/1661529/is-meyers-implementation-of-the-singleton-pattern-thread-safe
Communicator& Communicator::getInstance()
{
    static Communicator comm;
    return comm;
}

void Communicator::exchngBottomToTop(float *sbuf, float *rbuf, int count)
{
    MPI_Sendrecv(sbuf, count, MPI_FLOAT, nbrtop, 0, rbuf, count, MPI_FLOAT, nbrbottom, 0, comm1d, status);
}
void Communicator::exchngTopToBottom(float *sbuf, float *rbuf, int count)
{
    MPI_Sendrecv(sbuf, count, MPI_FLOAT, nbrbottom, 0, rbuf, count, MPI_FLOAT, nbrtop, 0, comm1d, status);
}
void Communicator::waitAll() { MPI_Waitall(4, request, status); }
void Communicator::exchngData(float *sbuf_t, float *rbuf_t, float *sbuf_b, float *rbuf_b, int count)
{
    MPI_Sendrecv(sbuf_t, count, MPI_FLOAT, nbrtop, 0, rbuf_t, count, MPI_FLOAT, nbrbottom, 0, comm1d, status);
    MPI_Sendrecv(sbuf_b, count, MPI_FLOAT, nbrbottom, 0, rbuf_b, count, MPI_FLOAT, nbrtop, 0, comm1d, status);
}
void Communicator::exchngDataNB(float *sbuf_t, int count_st, float *rbuf_t, int count_rt, float *sbuf_b, int count_sb,
                                float *rbuf_b, int count_rb)
{
    MPI_Irecv(rbuf_t, count_rt, MPI_FLOAT, nbrbottom, 0, comm1d, &request[0]);
    MPI_Irecv(rbuf_b, count_rb, MPI_FLOAT, nbrtop, 0, comm1d, &request[1]);
    MPI_Isend(sbuf_t, count_st, MPI_FLOAT, nbrtop, 0, comm1d, &request[2]);
    MPI_Isend(sbuf_b, count_sb, MPI_FLOAT, nbrbottom, 0, comm1d, &request[3]);
    MPI_Waitall(4, request, status);
}
//////////////////////////////////////////////////////////////////////////
// Crap by Martin Sch.
void Communicator::exchngDataGPU(real *sbuf, int count_s, real *rbuf, int count_r, int nb_rank)
{
    MPI_Status MSstatus;
    MPI_Send(sbuf, count_s, MPI_Type_GPU, nb_rank, 0, commGPU);
    MPI_Recv(rbuf, count_r, MPI_Type_GPU, nb_rank, 0, commGPU, &MSstatus);
    ////test only - please don't use
    // MPI_Sendrecv(sbuf, count_s, MPI_Type_GPU, nb_rank, 0, rbuf, count_r, MPI_Type_GPU, nb_rank, 0, comm1d,
    // MPI_STATUSES_IGNORE);
}
void Communicator::sendRecvGPU(real *sbuf, int count_s, real *rbuf, int count_r, int nb_rank)
{
    // test only - please don't use
    MPI_Sendrecv(sbuf, count_s, MPI_Type_GPU, nb_rank, 0, rbuf, count_r, MPI_Type_GPU, nb_rank, 0, commGPU,
                 MPI_STATUSES_IGNORE);
}
void Communicator::nbRecvDataGPU(real *rbuf, int count_r, int nb_rank)
{
    // printf("\n Start Recv Rank: %d, neighbor Rank: %d, request = %d \n", PID, nb_rank, (int)requestGPU.size());
    // fflush(stdout);

    requestGPU.push_back(0);
    MPI_Irecv(rbuf, count_r, MPI_Type_GPU, nb_rank, 0, commGPU, &requestGPU[rcount]);
    rcount++;

    // printf("\n End Recv - Rank: %d , neighbor Rank: %d \n", PID, nb_rank);
    // fflush(stdout);
}
void Communicator::nbSendDataGPU(real *sbuf, int count_s, int nb_rank)
{
    // printf("\n Start Send Rank: %d, neighbor Rank: %d, request = %d \n", PID, nb_rank, (int)requestGPU.size());
    // fflush(stdout);

    requestGPU.push_back(0);
    MPI_Isend(sbuf, count_s, MPI_Type_GPU, nb_rank, 0, commGPU, &requestGPU[rcount]);
    rcount++;

    // printf("\n End Send - Rank: %d , neighbor Rank: %d \n", PID, nb_rank);
    // fflush(stdout);
}
void Communicator::waitallGPU()
{
    // printf("\n Start Waitall Rank: %d, request = %d \n", PID, (int)requestGPU.size());
    // fflush(stdout);
    if (requestGPU.size() > 0) {
        MPI_Waitall(static_cast<int>(requestGPU.size()), &requestGPU[0], MPI_STATUSES_IGNORE);
        requestGPU.resize(0);
        rcount = 0;
    }
    // printf("\n End Waitall \n");
    // fflush(stdout);
}
void Communicator::sendDataGPU(real *sbuf, int count_s, int nb_rank)
{
    MPI_Send(sbuf, count_s, MPI_Type_GPU, nb_rank, 0, commGPU);
}
void Communicator::waitGPU(int id) { MPI_Wait(&requestGPU[id], MPI_STATUSES_IGNORE); }
void Communicator::resetRequest()
{
    if (requestGPU.size() > 0) {
        requestGPU.resize(0);
        rcount = 0;
    }
}
void Communicator::barrierGPU()
{
    // printf("\n Start Waitall Rank: %d, request = %d \n", PID, (int)requestGPU.size());
    // fflush(stdout);
    if (requestGPU.size() > 0) {
        MPI_Barrier(commGPU);
    }
    // printf("\n End Waitall \n");
    // fflush(stdout);
}
void Communicator::barrier() { MPI_Barrier(commGPU); }

//////////////////////////////////////////////////////////////////////////
void Communicator::exchngDataGeo(int *sbuf_t, int *rbuf_t, int *sbuf_b, int *rbuf_b, int count)
{
    MPI_Irecv(rbuf_t, count, MPI_INT, nbrbottom, 0, comm1d, &request[0]);
    MPI_Irecv(rbuf_b, count, MPI_INT, nbrtop, 0, comm1d, &request[1]);
    MPI_Isend(sbuf_t, count, MPI_INT, nbrtop, 0, comm1d, &request[2]);
    MPI_Isend(sbuf_b, count, MPI_INT, nbrbottom, 0, comm1d, &request[3]);
    MPI_Waitall(4, request, status);
}
int Communicator::getPID() const { return PID; }
int Communicator::getNummberOfProcess() const { return numprocs; }
int Communicator::getNeighbourTop() { return nbrtop; }
int Communicator::getNeighbourBottom() { return nbrbottom; }
MPI_Comm Communicator::getCommunicator() { return comm1d; }
void Communicator::startTimer() { starttime = MPI_Wtime(); }
void Communicator::stopTimer() { endtime = MPI_Wtime(); }
double Communicator::getTime() { return endtime - starttime; }
void Communicator::distributeGeometry(unsigned int *dataRoot, unsigned int *dataNode, int dataSizePerNode)
{
    MPI_Scatter(dataRoot, dataSizePerNode, MPI_UNSIGNED, dataNode, dataSizePerNode, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
}
int Communicator::mapCudaDevice(const int &rank, const int &size, const std::vector<unsigned int> &devices,
                                const int &maxdev)
{
    int device        = -1;
    char *host        = (char *)malloc(sizeof(char) * size * 255);
    unsigned int *map = (unsigned int *)malloc(sizeof(unsigned int) * size);

    char hostname[255];
    gethostname(hostname, 254);
    hostname[254] = 0;

    MPI_Gather(hostname, 255, MPI_BYTE, host, 255, MPI_BYTE, 0, MPI_COMM_WORLD);

    int i, j;
    if (rank == 0) {
        for (i = 0; i < size; i++) {
            int counter = 0;
            for (j = 0; j < i; j++) {
                if (strcmp(&host[i * 255], &host[j * 255]) == 0)
                    counter++;
            }
            if (counter >= maxdev) {
                VF_LOG_CRITICAL("More processes than GPUs!");
                exit(1);
            }
            map[i] = devices[counter];
        }
    }

    MPI_Scatter(map, 1, MPI_UNSIGNED, &device, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    VF_LOG_INFO("Rank: {} runs on host: {} with GPU: {}", rank, hostname, device);

    free(map);
    free(host);
    return device;
}

}

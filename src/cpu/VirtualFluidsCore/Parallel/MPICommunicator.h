#if defined VF_MPI

#ifndef MPICOMMUNICATOR_H
#define MPICOMMUNICATOR_H

#include "Communicator.h"
#include <PointerDefinitions.h>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbLogger.h>
#include <mpi.h>
#include <string>
#include <vector>

//! \brief A class uses MPI library to communication.
//! \details Support MPI communication. Implements singleton pattern.
//! \author K. Kutscher

class MPICommunicator : public Communicator
{
private:
    MPICommunicator();
    MPICommunicator(const MPICommunicator &) {}

public:
    ~MPICommunicator() override;
    static SPtr<Communicator> getInstance();
    int getBundleID() override;
    int getNumberOfBundles() override;
    int getProcessID() override;
    int getProcessID(int bundle, int rank) override;
    int getNumberOfProcesses() override;
    void *getNativeCommunicator() override;
    int getRoot() override;
    int getBundleRoot() override;
    int getProcessRoot() override;
    int getNumberOfProcessesInBundle(int bundle) override;
    bool isRoot() override;
    void abort(int errorcode) override;

    void sendSerializedObject(std::stringstream &ss, int target) override;
    void receiveSerializedObject(std::stringstream &ss, int source) override;

    void barrier() override;

    std::vector<std::string> gather(const std::string &str) override;
    std::vector<int> gather(std::vector<int> &values) override;
    std::vector<float> gather(std::vector<float> &values) override;
    std::vector<double> gather(std::vector<double> &values) override;
    std::vector<unsigned long long> gather(std::vector<unsigned long long> &values) override;

    void allGather(std::vector<int> &svalues, std::vector<int> &rvalues) override;
    void allGather(std::vector<float> &svalues, std::vector<float> &rvalues) override;
    void allGather(std::vector<double> &svalues, std::vector<double> &rvalues) override;
    void allGather(std::vector<unsigned long long> &svalues, std::vector<unsigned long long> &rvalues) override;

    void broadcast(int &value) override;
    void broadcast(float &value) override;
    void broadcast(double &value) override;
    void broadcast(long int &value) override;
    void broadcast(std::vector<int> &values) override;
    void broadcast(std::vector<float> &values) override;
    void broadcast(std::vector<double> &values) override;
    void broadcast(std::vector<long int> &values) override;

    template <class T>
    std::vector<T> gather(std::vector<T> &values);

    template <class T>
    void allGather(std::vector<T> &svalues, std::vector<T> &rvalues);

    template <class T>
    void broadcast(std::vector<T> &values);

    template <class T>
    void broadcast(T &value);

private:
    int numprocs, PID;
    MPI_Comm comm;
    int root;
};

//////////////////////////////////////////////////////////////////////////
template <class T>
std::vector<T> MPICommunicator::gather(std::vector<T> &values)
{
    MPI_Datatype mpiDataType;
    if ((std::string) typeid(T).name() == (std::string) typeid(double).name())
        mpiDataType = MPI_DOUBLE;
    else if ((std::string) typeid(T).name() == (std::string) typeid(float).name())
        mpiDataType = MPI_FLOAT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(int).name())
        mpiDataType = MPI_INT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(unsigned long long).name())
        mpiDataType = MPI_UNSIGNED_LONG_LONG;
    else
        throw UbException(UB_EXARGS, "no MpiDataType for T" + (std::string) typeid(T).name());

    int count = static_cast<int>(values.size());
    std::vector<T> rvalues(1);

    if (PID == root) {
        rvalues.resize(numprocs * count);
    }

    MPI_Gather(&values[0], count, mpiDataType, &rvalues[0], count, mpiDataType, root, comm);

    return rvalues;
}
//////////////////////////////////////////////////////////////////////////
template <class T>
void MPICommunicator::allGather(std::vector<T> &svalues, std::vector<T> &rvalues)
{
    MPI_Datatype mpiDataType;
    if ((std::string) typeid(T).name() == (std::string) typeid(double).name())
        mpiDataType = MPI_DOUBLE;
    else if ((std::string) typeid(T).name() == (std::string) typeid(float).name())
        mpiDataType = MPI_FLOAT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(int).name())
        mpiDataType = MPI_INT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(unsigned long long).name())
        mpiDataType = MPI_UNSIGNED_LONG_LONG;
    else
        throw UbException(UB_EXARGS, "no MpiDataType for T" + (std::string) typeid(T).name());

    int scount;
    std::vector<int> displs, rcounts;

    scount = (int)(svalues.size());

    rcounts.resize(numprocs);
    MPI_Allgather(&scount, 1, MPI_INT, &rcounts[0], 1, MPI_INT, comm);
    displs.resize(numprocs);

    displs[0] = 0;
    for (int i = 1; i < numprocs; ++i) {
        displs[i] = displs[i - 1] + rcounts[i - 1];
    }

    rvalues.resize(displs[numprocs - 1] + rcounts[numprocs - 1]);

    if (rvalues.size() == 0) {
        rvalues.resize(1);
        rvalues[0] = 999;
    }
    if (scount == 0) {
        svalues.resize(1);
        svalues[0] = 999;
    }

    MPI_Allgatherv(&svalues[0], scount, mpiDataType, &rvalues[0], &rcounts[0], &displs[0], mpiDataType, comm);
}
//////////////////////////////////////////////////////////////////////////
template <class T>
void MPICommunicator::broadcast(std::vector<T> &values)
{
    MPI_Datatype mpiDataType;
    if ((std::string) typeid(T).name() == (std::string) typeid(double).name())
        mpiDataType = MPI_DOUBLE;
    else if ((std::string) typeid(T).name() == (std::string) typeid(float).name())
        mpiDataType = MPI_FLOAT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(int).name())
        mpiDataType = MPI_INT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(long int).name())
        mpiDataType = MPI_LONG_INT;
    else
        throw UbException(UB_EXARGS, "no MpiDataType for T" + (std::string) typeid(T).name());

    int rcount;
    if (this->PID == this->root) {
        rcount = (int)values.size();
    }

    MPI_Bcast(&rcount, 1, MPI_INT, this->root, comm);

    if (this->PID != this->root) {
        values.resize(rcount);
    }

    MPI_Bcast(&values[0], (int)values.size(), mpiDataType, this->root, comm);
}
//////////////////////////////////////////////////////////////////////////
template <class T>
void MPICommunicator::broadcast(T &value)
{
    MPI_Datatype mpiDataType;
    if ((std::string) typeid(T).name() == (std::string) typeid(double).name())
        mpiDataType = MPI_DOUBLE;
    else if ((std::string) typeid(T).name() == (std::string) typeid(float).name())
        mpiDataType = MPI_FLOAT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(int).name())
        mpiDataType = MPI_INT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(long int).name())
        mpiDataType = MPI_LONG_INT;
    else
        throw UbException(UB_EXARGS, "no MpiDataType for T" + (std::string) typeid(T).name());

    MPI_Bcast(&value, 1, mpiDataType, this->root, comm);
}
//////////////////////////////////////////////////////////////////////////

#endif

#endif

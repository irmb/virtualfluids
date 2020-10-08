#if defined VF_MPI

#include "MPICommunicator.h"
#include <mpi.h>

#include <sstream>
using namespace std;
//////////////////////////////////////////////////////////////////////////
MPICommunicator::MPICommunicator()
{
   //proof if MPI is initialized 
   int mpiInitialized = (int)false;
   MPI_Initialized(&mpiInitialized);
   if (!mpiInitialized)
   {
      MPI_Init(NULL, NULL);	
      //MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, NULL);
   }
   MPI_Comm_rank(MPI_COMM_WORLD, &PID);
   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
   //numprocs = 1000;
   comm = MPI_COMM_WORLD;
   root = 0;
}
//////////////////////////////////////////////////////////////////////////
MPICommunicator::~MPICommunicator()
{
   //proof if MPI is finalized
   int _mpiFinalized = (int)false;
   MPI_Finalized(&_mpiFinalized);
   if (!_mpiFinalized)
   {
      MPI_Finalize();
      //UBLOG(logINFO, "MPI_Finalize()");
   }
 }
//////////////////////////////////////////////////////////////////////////
SPtr<Communicator> MPICommunicator::getInstance()
{
   if( !Communicator::instance )
      Communicator::instance = SPtr<Communicator>(new MPICommunicator());
   return Communicator::instance;
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::abort(int errorcode)
{
   MPI_Abort(comm, errorcode);
}
////////////////////////////////////////////////////////////////////////////
vector<string> MPICommunicator::gather(const string& str)
{
   vector<string> parts;
   vector<string> strings;
   int scount;
   vector<char> rbuf(1);
   vector<int> rcounts(1);
   MPI_Status status;

   if (PID == root)
   {
      rcounts.resize(numprocs - 1);
      strings.push_back(str);

      for (int i = 1; i < numprocs; i++)
      {
         MPI_Recv(&rcounts[i-1], 1, MPI_INT, i, 0, comm, &status);
      }
      for (int i = 1; i < numprocs; i++)
      {
         rbuf.resize(rcounts[i-1]);
         MPI_Recv(&rbuf[0], rcounts[i-1], MPI_CHAR, i, 0, comm, &status);
         string s(&rbuf[0], rcounts[i-1]);
         if (s != "") strings.push_back(s);
      }
   }
   else
   {
      scount = (int)str.length();
      MPI_Send(&scount, 1, MPI_INT, root, 0, comm);
      MPI_Send((char *)str.c_str(), scount, MPI_CHAR, root, 0, comm);
   }
   return strings;
}
//////////////////////////////////////////////////////////////////////////
vector<int> MPICommunicator::gather(vector<int>& values)
{
   return gather<int>(values);
}
//////////////////////////////////////////////////////////////////////////
vector<float> MPICommunicator::gather(vector<float>& values)
{
   return gather<float>(values);
}
//////////////////////////////////////////////////////////////////////////
vector<double> MPICommunicator::gather(vector<double>& values)
{
   return gather<double>(values);
}
//////////////////////////////////////////////////////////////////////////
std::vector<unsigned long long> MPICommunicator::gather(std::vector<unsigned long long>& values)
{
   return gather<unsigned long long>(values);
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getProcessID()
{
   return PID;
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getProcessID(int  /*bundle*/, int  /*rank*/)
{
   return PID;
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getNumberOfProcesses()
{
   return numprocs;
}
//////////////////////////////////////////////////////////////////////////
void* MPICommunicator::getNativeCommunicator()
{
   return &comm;
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getBundleID()
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getNumberOfBundles()
{
   return 1;
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getRoot() 
{
   return root;
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getBundleRoot() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getProcessRoot() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getNumberOfProcessesInBundle(int  /*bundle*/)
{
   return numprocs;
}
//////////////////////////////////////////////////////////////////////////
bool MPICommunicator::isRoot()
{
   return PID == root;
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::sendSerializedObject( std::stringstream& ss, int target)
{
   string str = ss.str();
   int scount = static_cast<int> (str.length());
   MPI_Send(&scount,1,MPI_INT,target,0,comm);
   MPI_Send((char *)str.c_str(),scount,MPI_CHAR,target,0,comm);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::receiveSerializedObject( std::stringstream& ss, int source )
{
   vector<char> rbuf;
   int rcount;
   MPI_Status status;
   MPI_Recv(&rcount,1,MPI_INT,source,0,comm,&status);
   rbuf.resize(rcount);
   MPI_Recv(&rbuf[0],rcount,MPI_CHAR,source,0,comm,&status);
   ss.rdbuf()->pubsetbuf(&rbuf[0],rcount);
   string str (&rbuf[0]);
   ss.str(str);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::barrier()
{
   MPI_Barrier(comm);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::allGather(std::vector<int>& svalues, std::vector<int>& rvalues)
{
   allGather<int>(svalues, rvalues);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::allGather(std::vector<float>& svalues, std::vector<float>& rvalues)
{
   allGather<float>(svalues, rvalues);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::allGather(std::vector<double>& svalues, std::vector<double>& rvalues)
{
   allGather<double>(svalues, rvalues);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::allGather(std::vector<unsigned long long>& svalues, std::vector<unsigned long long>& rvalues)
{
   allGather<unsigned long long>(svalues, rvalues);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(std::vector<int>& values)
{
   broadcast<int>(values);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(std::vector<float>& values)
{
   broadcast<float>(values);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(std::vector<double>& values)
{
   broadcast<double>(values);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(std::vector<long int>& values)
{
   broadcast<long int>(values);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(int& value)
{
   broadcast<int>(value);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(float& value)
{
   broadcast<float>(value);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(double& value)
{
   broadcast<double>(value);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(long int& value)
{
   broadcast<long int>(value);
}
#endif

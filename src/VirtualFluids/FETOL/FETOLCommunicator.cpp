#if defined VF_MPI && defined VF_FETOL

#include "FETOLCommunicator.h"
#include <MPICommunicator.h>
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbException.h>

#include <boost/foreach.hpp>

#include <fbond.h>
#include <Version.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////
FETOLCommunicator::FETOLCommunicator()
{
   bond::init();

   mpi_comm = MPICommunicator::getInstance();
   mpi_rank = mpi_comm->getProcessID();
   mpi_size = mpi_comm->getNumberOfProcesses();

   BID = bond::bundleID();
   numberOfBundles = bond::bundleCount();
   numberOfProcesses = bond::processCount();
   PID = bond::processID();
   bundleSizes.resize(numberOfBundles);
   bond::bundleSizes(numberOfBundles, &bundleSizes[0]);

   root = 0;
}
//////////////////////////////////////////////////////////////////////////
FETOLCommunicator::~FETOLCommunicator()
{
   bond::finalize();
}
//////////////////////////////////////////////////////////////////////////
CommunicatorPtr FETOLCommunicator::getInstance()
{
   if( !Communicator::instance )
      Communicator::instance = CommunicatorPtr(new FETOLCommunicator());
   return Communicator::instance;
}//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getBundleID()
{
   return bond::bundleID();
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getNumberOfBundles()
{
   return bond::bundleCount();
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getBONDProcessID()
{
   return bond::processID();
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getBONDProcessCount()
{
   return bond::processCount();
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getProcessID()
{
   return PID;
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getProcessID(int bundle, int rank)
{
   return bond::processID_for_bundleID_and_MPIRank(bundle, rank);
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getNumberOfProcesses()
{
   return numberOfProcesses;
}
//////////////////////////////////////////////////////////////////////////
void* FETOLCommunicator::getNativeCommunicator()
{
   return mpi_comm->getNativeCommunicator();
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getRoot() 
{
   return root;
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getBundleRoot() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getProcessRoot() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getNumberOfProcessesInBundle(int bundle)
{
   return bundleSizes[bundle];
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getMPIRank()
{
   return mpi_rank;
}
//////////////////////////////////////////////////////////////////////////
int FETOLCommunicator::getMPICommSize()
{
   return mpi_size;
}
//////////////////////////////////////////////////////////////////////////
void FETOLCommunicator::broadcastInts(std::vector<int>& values)
{
   int rcount;

   if (PID == root) 
   {
      UBLOG(logDEBUG5, "bond::send-start: "<<"Bundle ID="<<BID<<", PID="<<PID<<", MPI rank="<<mpi_rank);
      rcount = (int)values.size();

      for(int i = 1; i < numberOfProcesses; i++)
      {
         bond::sendComplete(&rcount, 1, MPI_INT, i, 0);
      }
      if (rcount > 0)
      {
         for(int i = 1; i < numberOfProcesses; i++)
         {
            bond::sendComplete(&values[0], rcount, MPI_INT, i, 0);
         }
      }
      UBLOG(logDEBUG5, "bond::send-end: "<<"Bundle ID="<<BID<<", PID="<<PID<<", MPI rank="<<mpi_rank);
   }
   else 
   {
      UBLOG(logDEBUG5, "bond::receive-start: "<<"Bundle ID="<<BID<<", PID="<<PID<<", MPI rank="<<mpi_rank);
      bond::receiveComplete(&rcount, 1, MPI_INT, root, 0);
      if (rcount > 0)
      {
         values.resize(rcount);
         bond::receiveComplete(&values[0], (int)values.size(), MPI_INT, root, 0);
      }
      UBLOG(logDEBUG5, "bond::receive-end: "<<"Bundle ID="<<BID<<", PID="<<PID<<", MPI rank="<<mpi_rank);
   }
}
//////////////////////////////////////////////////////////////////////////
std::vector<std::string> FETOLCommunicator::gatherStrings(const std::string& str)
{
   vector<string> parts;
   vector<string> strings;
   int scount;
   vector<char> rbuf(1);
   vector<int> rcounts(1);

   if (PID == root)
   {
      rcounts.resize(numberOfProcesses - 1);
      strings.push_back(str);

      for (int i = 1; i < numberOfProcesses; i++)
      {
         bond::receiveComplete(&rcounts[i-1], 1, MPI_INT, i, 0);
      }
      for (int i = 1; i < numberOfProcesses; i++)
      {
         rbuf.resize(rcounts[i-1]);
         bond::receiveComplete(&rbuf[0], rcounts[i-1], MPI_CHAR, i, 0);
         string s(&rbuf[0], rcounts[i-1]);
         if (s != "") strings.push_back(s);
      }
   }
   else
   {
      scount = (int)str.length();
      bond::sendComplete(&scount, 1, MPI_INT, root, 0);
      bond::sendComplete((char *)str.c_str(), scount, MPI_CHAR, root, 0);
   }
   return strings;
}
//////////////////////////////////////////////////////////////////////////
std::vector<double> FETOLCommunicator::gatherDoubles(std::vector<double>& values) 
{
   return mpi_comm->gatherDoubles(values);
}
//////////////////////////////////////////////////////////////////////////
void FETOLCommunicator::allGatherInts(std::vector<int>& svalues, std::vector<int>& rvalues)
{
   int scount;
   vector<int> rbuf;
   vector<int> rcounts;

   if (PID == root)
   {
      rvalues.resize(0);
      rvalues.insert(rvalues.end(), svalues.begin(), svalues.end());
      rcounts.resize(numberOfProcesses - 1);

      for (int i = 1; i < numberOfProcesses; i++)
      {
         bond::receiveComplete(&rcounts[i-1], 1, MPI_INT, i, 0);
      }
      for (int i = 1; i < numberOfProcesses; i++)
      {
         if(rcounts[i-1] > 0)
         {
            rbuf.resize(rcounts[i-1]);
            bond::receiveComplete(&rbuf[0], rcounts[i-1], MPI_INT, i, 0);
            rvalues.insert(rvalues.end(), rbuf.begin(), rbuf.end());
         }
      }
   }
   else
   {
      scount = (int)svalues.size();
      bond::sendComplete(&scount, 1, MPI_INT, root, 0);
      if(scount > 0)
         bond::sendComplete(&svalues[0], scount, MPI_INT, root, 0);
   }

   if(rvalues.size() > 0)
      broadcastInts(rvalues);
}
//////////////////////////////////////////////////////////////////////////
void FETOLCommunicator::sendSerializedObject(std::stringstream& ss, int target) 
{

}
//////////////////////////////////////////////////////////////////////////
void FETOLCommunicator::receiveSerializedObject(std::stringstream& ss, int source) 
{

}
//////////////////////////////////////////////////////////////////////////
void FETOLCommunicator::barrier()
{
   int flag = 0;

   if (PID == root)
   {
      for (int i = 1; i < numberOfProcesses; i++)
      {
         bond::receiveComplete(&flag, 1, MPI_INT, i, 0);
      }

      for (int i = 1; i < numberOfProcesses; i++)
      {
         bond::sendComplete(&flag, 1, MPI_INT, i, 0);
      }
   }
   else
   {
      bond::sendComplete(&flag, 1, MPI_INT, root, 0);
      bond::receiveComplete(&flag, 1, MPI_INT, root, 0);
   }
}

#endif

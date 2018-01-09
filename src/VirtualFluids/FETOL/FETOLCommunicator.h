#if defined VF_MPI && defined VF_FETOL

#ifndef FETOLCommunicator_H
#define FETOLCommunicator_H

#include <iostream>
#include <vector>
#include <string>

#include <MPICommunicator.h>

#include <memory>
class FETOLCommunicator;
typedef std::shared_ptr<FETOLCommunicator> BondCommunicatorPtr;

///\brief Communicator for agent-based framework BOND

class FETOLCommunicator : public Communicator
{
private:
   FETOLCommunicator();
public:
   ~FETOLCommunicator();
   static CommunicatorPtr getInstance();
   int getBundleID();
   int getNumberOfBundles();
   int getProcessID();
   int getProcessID(int bundle, int rank);
   int getNumberOfProcesses();
   int getBONDProcessID();
   int getBONDProcessCount();
   int getMPIRank();
   int getMPICommSize();
   int getRoot();
   int getBundleRoot();
   int getProcessRoot();
   int getNumberOfProcessesInBundle(int bundle);
   void* getNativeCommunicator();
   void broadcastInts(std::vector<int>& values);
   std::vector<std::string> gatherStrings(const std::string& str);
   std::vector<double> gatherDoubles(std::vector<double>& values); 
   void allGatherInts(std::vector<int>& svalues, std::vector<int>& rvalues);
   void sendSerializedObject(std::stringstream& ss, int target);
   void receiveSerializedObject(std::stringstream& ss, int source);
   void barrier();
private:
   CommunicatorPtr mpi_comm;
   int mpi_rank, mpi_size;
   int numberOfBundles, BID, numberOfProcesses, PID;
   int root;
   std::vector<int> bundleSizes;
};

#endif

#endif

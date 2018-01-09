#ifndef NullCommunicator_H
#define NullCommunicator_H

#include "Communicator.h"

#include <memory>

class NullCommunicator;
typedef std::shared_ptr<NullCommunicator> NullCommunicatorPtr;

class NullCommunicator : public Communicator
{
public:
   NullCommunicator();
   ~NullCommunicator();
   int getBundleID();
   int getNumberOfBundles();
   int getProcessID();
   int getNumberOfProcesses();
   void* getNativeCommunicator();
   int getRoot();
   int getBundleRoot();
   int getProcessRoot();
   std::vector<std::string> gather(const std::string& str);
   std::vector<double> gatherDoubles(std::vector<double>& values); 
   void allGatherInts(std::vector<int>& svalues, std::vector<int>& rvalues);
   void sendSerializedObject(std::stringstream& ss, int target);
   void receiveSerializedObject(std::stringstream& ss, int source);
protected:
private:
};

#endif

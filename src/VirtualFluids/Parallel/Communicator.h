#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>
#include <string>
#include <memory>

#include <VirtualFluidsDefinitions.h>

class Communicator;
typedef std::shared_ptr<Communicator> CommunicatorPtr;

class VF_PUBLIC Communicator
{
public:
   virtual ~Communicator(){}
   static CommunicatorPtr getInstance();
   virtual int getBundleID() = 0;
   virtual int getNumberOfBundles() = 0;
   virtual int getProcessID() = 0;
   virtual int getProcessID(int bundle, int rank) = 0;
   virtual int getNumberOfProcesses() = 0;
   virtual bool isRoot() = 0;
   virtual void* getNativeCommunicator() = 0;

   virtual void sendSerializedObject(std::stringstream& ss, int target) = 0;
   virtual void receiveSerializedObject(std::stringstream& ss, int source) = 0;

   virtual int getRoot() = 0;
   virtual int getBundleRoot() = 0;
   virtual int getProcessRoot() = 0;
   virtual int getNumberOfProcessesInBundle(int bundle) = 0;
   virtual void barrier() = 0;
   virtual void abort(int errorcode) = 0;

   virtual std::vector<std::string> gather(const std::string& str) = 0;
   virtual std::vector<int> gather(std::vector<int>& values) = 0;
   virtual std::vector<float> gather(std::vector<float>& values) = 0;
   virtual std::vector<double> gather(std::vector<double>& values) = 0;

   virtual void allGather(std::vector<int>& svalues, std::vector<int>& rvalues) = 0;
   virtual void allGather(std::vector<float>& svalues, std::vector<float>& rvalues) = 0;
   virtual void allGather(std::vector<double>& svalues, std::vector<double>& rvalues) = 0;
   
   virtual void broadcast(int& value) = 0;
   virtual void broadcast(float& value) = 0;
   virtual void broadcast(double& value) = 0;
   virtual void broadcast(std::vector<int>& values) = 0;
   virtual void broadcast(std::vector<float>& values) = 0;
   virtual void broadcast(std::vector<double>& values) = 0;
protected:
   Communicator(){}
   static CommunicatorPtr instance;
private:
};

#endif


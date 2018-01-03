#ifndef LOADBALANCER_H
#define LOADBALANCER_H

#if defined VF_MPI

#include <vector>
#include <map>
#include "Block3D.h"
#include "Communicator.h"
#include "MPICommunicator.h"
#include "Grid3D.h"
#include "LBMKernel.h"

struct NeighbourProcess
{
   int processID;
   double load;
};

struct ltweight
{
   bool operator()(const int l1, const int l2) const
   {
      return (l1 > l2);
   }
};

struct TransferBuffer
{
   int comFlag; //communication on = 1, communication off = 0
   int blocksSize;
   int sbufferSize;
   std::string sbuffer;
};

struct  TransferBlock
{
   double load;
   int weight;
   Block3DPtr block;
};


class LoadBalancer;
typedef std::shared_ptr<LoadBalancer> LoadBalancerPtr;

class LoadBalancer
{
public:
   // boundary blocks:    block, load
   //typedef std::map< Block3DPtr, int> BBlocks;
   typedef std::map< Block3DPtr, TransferBlock> BBlocks;
   //typedef std::set< TransferBlock > BBlocks;
   //           neighbor, boundary blocks
   typedef std::map< int, BBlocks > BBlocksMap;
   //transfer blocks: weight, transfer block, sort function
   typedef std::multimap< int, TransferBlock, ltweight> TBlocks;
   //          neighbor,buffer    
   typedef std::map<int, TransferBuffer> TransferBufferMap;
public:
   LoadBalancer(Grid3DPtr grid, CommunicatorPtr comm, const int& endDir);
   virtual ~LoadBalancer();
   bool balance();

protected:
private:
   void callculateLoad();
   void collectData();
   void collectNeighboursLoad();
   void prepareToSendLoad();
   void prepareToSendLoad(const int& nProcessID, const double& transLoad, int& conv);
   void prepareToRecieveLoad();
   void sendRecieveLoad();
   void saveRecievedLoad();
   void sendRecieveChanges();
   bool isConverged();
   void setConverged(int value);
   void addTransferBlock(Block3DPtr block, int neighBlockRank);
   double load;
   int processID;
   double loadThreshold;
   int loadThresholdP;
   int loadP;
   std::vector<NeighbourProcess> neighbourProcesses;
   Grid3DPtr grid;
   CommunicatorPtr comm;
   BBlocksMap neighboursMap;
   BBlocksMap::iterator itnMap;
   BBlocks::iterator itbBlocks;
   int endDir;

   MPI_Comm mpi_comm;

   TransferBufferMap sendBuffer;
   TransferBufferMap receiveBuffer;

   std::vector<int> removeBlocksSend;
   
   std::vector<int> converged;
};

#endif

#endif 

#ifndef Calculator2_H
#define Calculator2_H

#include "Grid3D.h"
#include "Block3D.h"
#include "Synchronizer.h"
#include "MathUtil.hpp"
#include "basics/utilities/UbScheduler.h"
#include "basics/utilities/UbTiming.h"
#include "LoadBalancer.h"


class Calculator2;
typedef boost::shared_ptr<Calculator2> Calculator2Ptr;

#include "CalculationManager.h"

class Calculator2 
{
public:
   Calculator2(Grid3DPtr grid, SynchronizerPtr sync, bool mainThread = true);
   virtual ~Calculator2(){}
   void calculate(const double& endTime, CalculationManagerPtr cm, boost::exception_ptr& error);
   void addBlock(Block3DPtr block);
   void initConnectors();
   void setVisScheduler(UbSchedulerPtr s);
   //double getCallculationTime();
   std::vector< std::vector< Block3DPtr > > getBlocks(); 
   void deleteBlocks();
   void deleteConnectors();
protected:
   void calculateBlocks(int startLevel, int maxInitLevel);
   void calculateBlocks(int minInitLevel, int maxInitLevel, int staggeredStep);
   void initConnectors(std::vector<Block3DConnectorPtr>& connectors);
   void initRemoteConnectors();
   void swapDistributions(int startLevel, int maxInitLevel);
   void exchangeBlockData(int startLevel, int maxInitLevel, bool invStep);
   void exchangeInterfaceBlockData(int startLevel, int maxInitLevel, bool invStep);
   void connectorsPrepare(std::vector< Block3DConnectorPtr >& connectors);
   void connectorsSend(std::vector< Block3DConnectorPtr >& connectors, bool invStep);
   void connectorsReceive(std::vector< Block3DConnectorPtr >& connectors, bool invStep);
   void connectorsSetInvStep(std::vector< Block3DConnectorPtr >& connectors, bool invStep);
   void interpolation(int startLevel, int maxInitLevel);
   void deleteConnectors(std::vector< std::vector< Block3DConnectorPtr > >& conns);
   void applyBCs(int startLevel, int maxInitLevel);
private:
   boost::barrier* bar;
   //double time;
   SynchronizerPtr sync;
   int minLevel, maxLevel;
   bool mainThread;
   bool refinement;
   Grid3DPtr grid;
   UbSchedulerPtr visScheduler;
   int calcStep;

   std::vector< std::vector<Block3DPtr> > blocks;
   std::vector< std::vector< Block3DConnectorPtr > > localConns;
   std::vector< std::vector< Block3DConnectorPtr > > remoteConns;

   std::vector< std::vector< Block3DConnectorPtr > > localInterfaceBlockConns;
   std::vector< std::vector< Block3DConnectorPtr > > remoteInterfaceBlockConns;

   //localInterConns and remoteInterConns save interpolation connectors 
   //every element save CF connectors for current level and FC connectors for next level
   //e.g. 
   //localInterConns[0] = CF(0), FC(1)
   //localInterConns[1] = CF(1), FC(2)
   //localInterConns[2] 
   std::vector< std::vector< Block3DConnectorPtr > > localInterConns;
   std::vector< std::vector< Block3DConnectorPtr > > remoteInterConns;

   //UbTimer timer, timer2, timer3;
   bool loadBalancingComp;

};

#endif

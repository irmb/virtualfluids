#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <memory>
#include <vector>

class Grid3D;
class Synchronizer;
class UbScheduler;
class CalculationManager;
class Block3D;
class Block3DConnector;
class TimeAveragedValuesCoProcessor;


#include <boost/thread/barrier.hpp>
#include <boost/exception_ptr.hpp>

class Calculator;
typedef std::shared_ptr<Calculator> CalculatorPtr;


class Calculator 
{
public:
   Calculator();
   Calculator(std::shared_ptr<Grid3D> grid, std::shared_ptr<Synchronizer> sync, bool mainThread = true);
   virtual ~Calculator(){}
   virtual void calculate(const double& endTime, std::shared_ptr<CalculationManager> cm, boost::exception_ptr& error);
   void addBlock(std::shared_ptr<Block3D> block);
   void initConnectors();
   void setVisScheduler(std::shared_ptr<UbScheduler> s);
   //double getCalculationTime();
   std::vector< std::vector< std::shared_ptr<Block3D> > > getBlocks() const;
   void deleteBlocks();
   void deleteConnectors();

   void setTimeAveragedValuesCoProcessor(std::shared_ptr<TimeAveragedValuesCoProcessor> coProcessor);

protected:
   void calculateBlocks(int startLevel, int maxInitLevel);
   void calculateBlocks(int minInitLevel, int maxInitLevel, int staggeredStep);
   void initConnectors(std::vector<std::shared_ptr<Block3DConnector> >& connectors);
   virtual void initRemoteConnectors();
   void swapDistributions(int startLevel, int maxInitLevel);
   virtual void exchangeBlockData(int startLevel, int maxInitLevel);
   void exchangeInterfaceBlockData(int startLevel, int maxInitLevel);
   virtual void connectorsPrepare(std::vector< std::shared_ptr<Block3DConnector> >& connectors);
   virtual void connectorsSend(std::vector< std::shared_ptr<Block3DConnector> >& connectors);
   virtual void connectorsReceive(std::vector< std::shared_ptr<Block3DConnector> >& connectors);
   void interpolation(int startLevel, int maxInitLevel);
   void deleteConnectors(std::vector< std::vector< std::shared_ptr<Block3DConnector> > >& conns);
   //void applyBCs(int startLevel, int maxInitLevel);
   void applyPreCollisionBC(int startLevel, int maxInitLevel);
   void applyPostCollisionBC(int startLevel, int maxInitLevel);
   int minLevel, maxLevel;
   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > localConns;
   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > remoteConns;
   std::shared_ptr<Synchronizer> sync;

   boost::barrier* bar;
   //double time;

   bool mainThread;
   bool refinement;
   std::shared_ptr<Grid3D> grid;
   std::shared_ptr<UbScheduler> visScheduler;
   int calcStep;
   std::vector< std::vector<std::shared_ptr<Block3D> > > blocks;


   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > localInterfaceBlockConns;
   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > remoteInterfaceBlockConns;

   //localInterConns and remoteInterConns save interpolation connectors 
   //every element save CF connectors for current level and FC connectors for next level
   //e.g. 
   //localInterConns[0] = CF(0), FC(1)
   //localInterConns[1] = CF(1), FC(2)
   //localInterConns[2] 
   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > localInterConns;
   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > remoteInterConns;

   //UbTimer timer, timer2, timer3;
   bool loadBalancingComp;

   std::shared_ptr<TimeAveragedValuesCoProcessor> taValuesCoProcessor;

};

#endif

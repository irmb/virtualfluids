#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <memory>
#include <vector>

class Grid3D;
class UbScheduler;
class Block3D;
class Block3DConnector;
class CoProcessor;

class Calculator;
typedef std::shared_ptr<Calculator> CalculatorPtr;


class Calculator 
{
public:
   Calculator();
   virtual ~Calculator();
   void setGrid(std::shared_ptr<Grid3D> grid);
   void setLastTimeStep(int t);
   void setVisScheduler(std::shared_ptr<UbScheduler> s);
   void addCoProcessor(std::shared_ptr<CoProcessor> coProcessor);
   void coProcess(double step);

   virtual void calculate()=0;
protected:
   virtual void initLocalConnectors();
   virtual void initRemoteConnectors();
   void initConnectors(std::vector<std::shared_ptr<Block3DConnector> >& connectors);
   void deleteBlocks();
   void deleteConnectors();
   void deleteConnectors(std::vector< std::vector< std::shared_ptr<Block3DConnector> > >& conns);

   int minLevel, maxLevel;
   int startTimeStep;
   int lastTimeStep;
   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > localConns;
   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > remoteConns;

   bool refinement;
   std::shared_ptr<Grid3D> grid;
   std::shared_ptr<UbScheduler> visScheduler;
   std::vector< std::vector<std::shared_ptr<Block3D> > > blocks;

   //localInterConns and remoteInterConns save interpolation connectors 
   //every element save CF connectors for current level and FC connectors for next level
   //e.g. 
   //localInterConns[0] = CF(0), FC(1)
   //localInterConns[1] = CF(1), FC(2)
   //localInterConns[2] 
   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > localInterConns;
   std::vector< std::vector< std::shared_ptr<Block3DConnector> > > remoteInterConns;

   std::vector< std::shared_ptr<CoProcessor> > coProcessors;
};

#endif

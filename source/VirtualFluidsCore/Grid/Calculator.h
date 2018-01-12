#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <PointerDefinitions.h>
#include <vector>

class Grid3D;
class UbScheduler;
class Block3D;
class Block3DConnector;
class CoProcessor;

class Calculator 
{
public:
   Calculator();
   virtual ~Calculator();
   void setGrid(SPtr<Grid3D> grid);
   void setLastTimeStep(int t);
   void setVisScheduler(SPtr<UbScheduler> s);
   void addCoProcessor(SPtr<CoProcessor> coProcessor);
   void coProcess(double step);

   virtual void calculate()=0;
protected:
   virtual void initLocalConnectors();
   virtual void initRemoteConnectors();
   void initConnectors(std::vector<SPtr<Block3DConnector> >& connectors);
   void deleteBlocks();
   void deleteConnectors();
   void deleteConnectors(std::vector< std::vector< SPtr<Block3DConnector> > >& conns);

   int minLevel, maxLevel;
   int startTimeStep;
   int lastTimeStep;
   std::vector< std::vector< SPtr<Block3DConnector> > > localConns;
   std::vector< std::vector< SPtr<Block3DConnector> > > remoteConns;

   bool refinement;
   SPtr<Grid3D> grid;
   SPtr<UbScheduler> visScheduler;
   std::vector< std::vector<SPtr<Block3D> > > blocks;

   //localInterConns and remoteInterConns save interpolation connectors 
   //every element save CF connectors for current level and FC connectors for next level
   //e.g. 
   //localInterConns[0] = CF(0), FC(1)
   //localInterConns[1] = CF(1), FC(2)
   //localInterConns[2] 
   std::vector< std::vector< SPtr<Block3DConnector> > > localInterConns;
   std::vector< std::vector< SPtr<Block3DConnector> > > remoteInterConns;

   std::vector< SPtr<CoProcessor> > coProcessors;
};

#endif

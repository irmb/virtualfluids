#ifndef OMPCALCULATOR_H
#define OMPCALCULATOR_H

#include "Calculator.h"

class Block3DConnector;

class OMPCalculator : public Calculator
{
public:
   OMPCalculator();
   virtual ~OMPCalculator();
   virtual void calculate();

protected:
   void calculateBlocks(int startLevel, int maxInitLevel, int calcStep);
   void swapDistributions(int startLevel, int maxInitLevel);
   void exchangeBlockData(int startLevel, int maxInitLevel);
   void connectorsPrepareLocal(std::vector< SPtr<Block3DConnector> >& connectors);
   void connectorsSendLocal(std::vector< SPtr<Block3DConnector> >& connectors);
   void connectorsReceiveLocal(std::vector< SPtr<Block3DConnector> >& connectors);
   void connectorsPrepareRemote(std::vector< SPtr<Block3DConnector> >& connectors);
   void connectorsSendRemote(std::vector< SPtr<Block3DConnector> >& connectors);
   void connectorsReceiveRemote(std::vector< SPtr<Block3DConnector> >& connectors);
   void interpolation(int startLevel, int maxInitLevel);
   void applyPreCollisionBC(int startLevel, int maxInitLevel);
   void applyPostCollisionBC(int startLevel, int maxInitLevel);
private:
};

#endif


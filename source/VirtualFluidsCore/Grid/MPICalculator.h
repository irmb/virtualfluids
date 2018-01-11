#ifndef MPICALCULATOR_H
#define MPICALCULATOR_H

#include "Calculator.h"

class MPICalculator;
typedef std::shared_ptr<MPICalculator> MPICalculatorPtr;

class Block3DConnector;

class MPICalculator  : public Calculator
{
public:
   MPICalculator();
   virtual ~MPICalculator();
   virtual void calculate();

protected:
   void calculateBlocks(int startLevel, int maxInitLevel, int calcStep);
   void swapDistributions(int startLevel, int maxInitLevel);
   virtual void exchangeBlockData(int startLevel, int maxInitLevel);
   virtual void connectorsPrepare(std::vector< std::shared_ptr<Block3DConnector> >& connectors);
   virtual void connectorsSend(std::vector< std::shared_ptr<Block3DConnector> >& connectors);
   virtual void connectorsReceive(std::vector< std::shared_ptr<Block3DConnector> >& connectors);
   void interpolation(int startLevel, int maxInitLevel);
   void applyPreCollisionBC(int startLevel, int maxInitLevel);
   void applyPostCollisionBC(int startLevel, int maxInitLevel);
private:


};

#endif


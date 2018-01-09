#ifndef MPICALCULATOR_H
#define MPICALCULATOR_H

#include "Grid3D.h"
#include "Block3D.h"
#include "Synchronizer.h"
#include "MathUtil.hpp"
#include "basics/utilities/UbScheduler.h"
#include "basics/utilities/UbTiming.h"
#include "LoadBalancer.h"
#include "TimeAveragedValuesCoProcessor.h"

#include "Calculator.h"

class MPICalculator;
typedef std::shared_ptr<MPICalculator> MPICalculatorPtr;

#include "CalculationManager.h"

class Block3DConnector;

class MPICalculator  : public Calculator
{
public:
   MPICalculator();
   MPICalculator(Grid3DPtr grid);
   virtual ~MPICalculator() {}
   virtual void calculate(const double& endTime, CalculationManagerPtr cm);

protected:
   void calculateBlocks(int startLevel, int maxInitLevel);
   void calculateBlocks(int minInitLevel, int maxInitLevel, int staggeredStep);
   void swapDistributions(int startLevel, int maxInitLevel);
   virtual void exchangeBlockData(int startLevel, int maxInitLevel);
   void exchangeInterfaceBlockData(int startLevel, int maxInitLevel);
   virtual void connectorsPrepare(std::vector< std::shared_ptr<Block3DConnector> >& connectors);
   virtual void connectorsSend(std::vector< std::shared_ptr<Block3DConnector> >& connectors);
   virtual void connectorsReceive(std::vector< std::shared_ptr<Block3DConnector> >& connectors);
   void interpolation(int startLevel, int maxInitLevel);
   void deleteConnectors(std::vector< std::vector< std::shared_ptr<Block3DConnector> > >& conns);
   void applyPreCollisionBC(int startLevel, int maxInitLevel);
   void applyPostCollisionBC(int startLevel, int maxInitLevel);
private:


};

#endif


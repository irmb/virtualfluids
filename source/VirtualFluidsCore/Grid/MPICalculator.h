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


class MPICalculator;
typedef boost::shared_ptr<MPICalculator> MPICalculatorPtr;

#include "CalculationManager.h"

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
   virtual void connectorsPrepare(std::vector< Block3DConnectorPtr >& connectors);
   virtual void connectorsSend(std::vector< Block3DConnectorPtr >& connectors);
   virtual void connectorsReceive(std::vector< Block3DConnectorPtr >& connectors);
   void interpolation(int startLevel, int maxInitLevel);
   void deleteConnectors(std::vector< std::vector< Block3DConnectorPtr > >& conns);
   void applyPostCollisionBC(int startLevel, int maxInitLevel);
private:


};

#endif


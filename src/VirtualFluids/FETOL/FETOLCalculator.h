#if defined VF_FETOL

#ifndef FETOLCALCULATOR_H
#define FETOLCALCULATOR_H

#include "Calculator.h"


class FETOLCalculator;
typedef boost::shared_ptr<FETOLCalculator> FETOLCalculator2Ptr;

#include "CalculationManager.h"

class FETOLCalculator : public Calculator
{
public:
   FETOLCalculator();
   FETOLCalculator(Grid3DPtr grid, SynchronizerPtr sync, bool mainThread = true);
   virtual ~FETOLCalculator(){}
   void calculate(const double& endTime, CalculationManagerPtr cm, boost::exception_ptr& error);
protected:
   void initRemoteConnectors();
   void exchangeFETOLBlockData(int startLevel, int maxInitLevel, bool invStep);
   //void connectorsPrepare(std::vector< Block3DConnectorPtr >& connectors);
   //void connectorsSend(std::vector< Block3DConnectorPtr >& connectors, bool invStep);
   //void connectorsReceive(std::vector< Block3DConnectorPtr >& connectors, bool invStep);
   //void bondConnectorsPrepare(std::vector< Block3DConnectorPtr >& connectors, bool invStep);
   //void bondConnectorsSend(std::vector< Block3DConnectorPtr >& connectors, bool invStep);
   //void bondConnectorsReceive(std::vector< Block3DConnectorPtr >& connectors, bool invStep);
   void initFETOLConnectors();
   void ifRestart(int startLevel, int maxInitLevel, bool invStep);
private:
   std::vector< std::vector< Block3DConnectorPtr > > remoteMPIConns;
   std::vector< std::vector< Block3DConnectorPtr > > remoteBondConns;
};

#endif

#endif

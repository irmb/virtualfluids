#ifndef PrePostBcCalculator_h
#define PrePostBcCalculator_h

#include "Grid3D.h"
#include "Block3D.h"
#include "Synchronizer.h"
#include "MathUtil.hpp"
#include "basics/utilities/UbScheduler.h"
#include "basics/utilities/UbTiming.h"
#include "LoadBalancer.h"

#include "Calculator.h"

class PrePostBcCalculator;
typedef std::shared_ptr<PrePostBcCalculator> PrePostBcCalculatorPtr;

#include "CalculationManager.h"

class PrePostBcCalculator : public Calculator
{
public:
   PrePostBcCalculator();
   PrePostBcCalculator(Grid3DPtr grid, SynchronizerPtr sync, bool mainThread = true);
   virtual ~PrePostBcCalculator(){}
   virtual void calculate(const double& endTime, CalculationManagerPtr cm, boost::exception_ptr& error);
protected:
   void applyPreCollisionBC(int startLevel, int maxInitLevel);
   void applyPostCollisionBC(int startLevel, int maxInitLevel);
private:

};

#endif

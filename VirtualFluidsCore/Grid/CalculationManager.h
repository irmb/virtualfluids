#ifndef CALCULATORMANAGER_H
#define CALCULATORMANAGER_H

#include "Grid3D.h"
#include "Communicator.h"

class CalculationManager;
typedef boost::shared_ptr<CalculationManager> CalculationManagerPtr;

#include "Calculator.h"

class CalculationManager : public boost::enable_shared_from_this<CalculationManager>
{
public:
   enum CalculatorType{MPI, FETOL, PrePostBc};
public:
   CalculationManager(Grid3DPtr grid, int numOfThreads, double endTime, UbSchedulerPtr visScheduler, CalculatorType calcType = CalculationManager::MPI);
   CalculationManager(Grid3DPtr grid, int numOfThreads, double endTime, UbSchedulerPtr visScheduler, 
                      CommunicatorPtr comm, int endDir, LBMReal nu, CalculatorType calcType = CalculationManager::MPI);
   CalculationManager(Grid3DPtr grid, int numOfThreads, double endTime, UbSchedulerPtr visScheduler, BoundaryConditionProcessorPtr bcProcessor, CalculatorType calcType = CalculationManager::MPI);
   virtual ~CalculationManager();
   void calculate();
   void setVisScheduler(UbSchedulerPtr s);
   bool balance();
private:
   void init();
   void calculateMain();
   void initCalcThreads();
   void reinitCalcThreads();
   void addBlocksToCalcThreads();
   CalculatorPtr createCalculator(Grid3DPtr grid, SynchronizerPtr sync, bool mainThread);
   Grid3DPtr grid;
   int numOfThreads;
   //boost::exception_ptr error;
   CalculatorType calcType;
   std::vector<CalculatorPtr> calcThreads;
   double endTime;
   UbSchedulerPtr visScheduler;
   CommunicatorPtr comm;   
   int endDir;
   LoadBalancerPtr loadBalancer;
   int rank;
   LBMReal nu;
   BoundaryConditionProcessorPtr bcProcessor;
}; 

#endif 


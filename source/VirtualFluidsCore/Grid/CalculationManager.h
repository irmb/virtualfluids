#ifndef CALCULATORMANAGER_H
#define CALCULATORMANAGER_H


#include <memory>
#include <vector>

class CalculationManager;
typedef std::shared_ptr<CalculationManager> CalculationManagerPtr;

class LoadBalancer;
class Communicator;
class Grid3D;
class Calculator;
class CalculatorFactory;
class TimeAveragedValuesCoProcessor;
enum class CalculatorType;


class CalculationManager : public std::enable_shared_from_this<CalculationManager>
{

public:
   CalculationManager(std::shared_ptr<Grid3D> grid, int numOfThreads, double endTime, std::shared_ptr<CalculatorFactory> calculatorFactory, CalculatorType type);
   CalculationManager(std::shared_ptr<Grid3D> grid, int numOfThreads, double endTime, std::shared_ptr<Communicator> comm, int endDir, std::shared_ptr<CalculatorFactory> calculatorFactory);
   virtual ~CalculationManager();

   void calculate();
   bool balance();
   void setTimeAveragedValuesCoProcessor(std::shared_ptr<TimeAveragedValuesCoProcessor> coProcessor);

private:
   void initCalcThreads();
   void reinitCalcThreads();
   void addBlocksToCalcThreads();

   std::shared_ptr<Grid3D> grid;
   int numOfThreads;
   double endTime;

   std::vector<std::shared_ptr<Calculator> > calcThreads;

   std::shared_ptr<LoadBalancer> loadBalancer;

   std::shared_ptr<CalculatorFactory> calculatorFactory;
   CalculatorType type;

}; 

#endif 


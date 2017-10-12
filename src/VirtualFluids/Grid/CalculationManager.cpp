#include "CalculationManager.h"

#include <iomanip>
#include <list>
#include <vector>

#include <boost/thread.hpp>
#include <boost/foreach.hpp>

#include <Calculator.h>
#include <MPICalculator.h>
#include <PrePostBcCalculator.h>
#include <Communicator.h>

#if defined VF_FETOL
#include <FETOLCalculator.h>
#endif

using namespace std;
//////////////////////////////////////////////////////////////////////////
CalculationManager::CalculationManager(Grid3DPtr grid, int numOfThreads, double endTime, UbSchedulerPtr visScheduler, CalculationManager::CalculatorType calcType)
                                       : grid(grid),
                                         numOfThreads(numOfThreads),
                                         endTime(endTime),
                                         visScheduler(visScheduler),
                                         calcType(calcType)
{
   init();
}
//////////////////////////////////////////////////////////////////////////
CalculationManager::CalculationManager(Grid3DPtr grid, int numOfThreads, double endTime, UbSchedulerPtr visScheduler, 
                                       CommunicatorPtr comm, int endDir, LBMReal nu, CalculatorType calcType)
                                       : grid(grid),
                                       numOfThreads(numOfThreads),
                                       endTime(endTime),
                                       visScheduler(visScheduler),
                                       calcType(calcType), 
                                       comm(comm),
                                       endDir(endDir),
                                       nu(nu)
{
   init();
   loadBalancer = LoadBalancerPtr(new LoadBalancer(grid, comm, endDir));
}
//////////////////////////////////////////////////////////////////////////
CalculationManager::~CalculationManager()
{
}
//////////////////////////////////////////////////////////////////////////
void CalculationManager::init()
{
   this->rank = grid->getRank();
   initCalcThreads();
}
//////////////////////////////////////////////////////////////////////////
void CalculationManager::calculate()
{
   try
   {
      boost::shared_ptr<CalculationManager> this_ = shared_from_this();
      if (calcType == CalculationManager::Hybrid || calcType == CalculationManager::PrePostBc)
      {
         boost::thread_group threads;
         boost::exception_ptr error;

         for (int i = 1; i<calcThreads.size(); i++)
         {
            threads.create_thread(boost::bind(&Calculator::calculate, calcThreads[i], endTime, this_, boost::ref(error)));
         }

         calcThreads[0]->calculate(endTime, this_, boost::ref(error));

         threads.join_all();
      } 
      else
      {
         boost::dynamic_pointer_cast<MPICalculator>(calcThreads[0])->calculate(endTime, this_);
      }


      //if( error )
      //{
      //   boost::rethrow_exception(error);
      //}
   }
   catch(std::exception& e)
   {
      UBLOG(logERROR, e.what());
      //throw e;
      exit (EXIT_FAILURE);
   }
}
//////////////////////////////////////////////////////////////////////////
void CalculationManager::initCalcThreads()
{
   UBLOG(logDEBUG1, "CalculationManager::initCalcThreads() - started");
   
   SynchronizerPtr sync( new Synchronizer(numOfThreads));
   
   for (int i =0; i<numOfThreads; i++)
   {
      //CalculatorPtr calc = CalculatorPtr (new Calculator(grid, sync, i == 0 ? true : false));
      CalculatorPtr calc = createCalculator(grid, sync, i == 0 ? true : false);
      calc->setVisScheduler(visScheduler);
      calcThreads.push_back(calc);
   }

   UBLOG(logDEBUG5, "calcThreads - initialized");

   addBlocksToCalcThreads();

   UBLOG(logDEBUG5, "calcThreads - filled with Blocks");

   //BOOST_FOREACH(CalculatorPtr calc,  calcThreads)
   //{
   //   calc->initConnectors();
   //}
   UBLOG(logDEBUG1, "CalculationManager::initCalcThreads() - stoped");
}
//////////////////////////////////////////////////////////////////////////
void CalculationManager::setVisScheduler(UbSchedulerPtr s)
{
   visScheduler = s;
}
//////////////////////////////////////////////////////////////////////////
bool CalculationManager::balance()
{
   //if(loadBalancer->balance())
   //   return true;
   //else
   //{
   //   grid->deleteConnectors();

   //   BOOST_FOREACH(CalculatorPtr calc,  calcThreads)
   //   {
   //      calc->deleteConnectors();
   //   }

   //   TbCbVectorMpiPool<LBMReal>::eraseMap();

   //   D3Q27SetInterfaceBlocksPatchVisitor setInterfaceBlocksPatchVisitor(nue, comm);
   //   grid->accept(setInterfaceBlocksPatchVisitor);

   //   D3Q27SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR);
   //   grid->accept( setConnsVisitor );

   //   PQueuePartitioningPatchVisitor pqPartVisitor(numOfThreads);
   //   grid->accept(pqPartVisitor);

   //   reinitCalcThreads();
      return false;
   //}
}
//////////////////////////////////////////////////////////////////////////
void CalculationManager::reinitCalcThreads()
{
   BOOST_FOREACH(CalculatorPtr c, calcThreads)
   {
      c->deleteBlocks();
   }

   addBlocksToCalcThreads();

   BOOST_FOREACH(CalculatorPtr calc,  calcThreads)
   {
      calc->initConnectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void CalculationManager::addBlocksToCalcThreads()
{
   int gridRank = grid->getRank();
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, true, blockVector);
      BOOST_FOREACH(Block3DPtr block, blockVector)
      {
         if (block)
         {
            calcThreads[block->getPart()]->addBlock(block);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
CalculatorPtr CalculationManager::createCalculator(Grid3DPtr grid, SynchronizerPtr sync, bool mainThread)
{
   switch (calcType)
   {
   case CalculationManager::Hybrid:
      return CalculatorPtr(new Calculator(grid, sync, mainThread));
   case CalculationManager::MPI:
      return CalculatorPtr (new MPICalculator(grid));
#if defined VF_FETOL
   case CalculationManager::FETOL:
      return CalculatorPtr (new FETOLCalculator(grid, sync, mainThread));
#endif
   case CalculationManager::PrePostBc:
      return CalculatorPtr(new PrePostBcCalculator(grid, sync, mainThread));
   default:
      UB_THROW(UbException(UB_EXARGS,"This Calculator is not defined!"));
   }
}
//////////////////////////////////////////////////////////////////////////
void CalculationManager::setTimeAveragedValuesCoProcessor(TimeAveragedValuesCoProcessorPtr coProcessor)
{
   calcThreads[0]->setTimeAveragedValuesCoProcessor(coProcessor);
}



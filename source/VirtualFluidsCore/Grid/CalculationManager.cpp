#include "CalculationManager.h"

#include <boost/thread.hpp>

#include "CalculatorFactory.h"

#include <Calculator.h>
#include <MPICalculator.h>
#if defined VF_FETOL
#include <FETOLCalculator.h>
#endif

#include <Communicator.h>
#include "TimeAveragedValuesCoProcessor.h"
#include "Grid3D.h"
#include "LoadBalancer.h"


//////////////////////////////////////////////////////////////////////////
CalculationManager::CalculationManager(Grid3DPtr grid, int numOfThreads, double endTime, std::shared_ptr<CalculatorFactory> calculatorFactory, CalculatorType type)
    : grid(grid),
    numOfThreads(numOfThreads),
    endTime(endTime),
    calculatorFactory(calculatorFactory),
    type(type)
{
    this->initCalcThreads();
}
//////////////////////////////////////////////////////////////////////////
CalculationManager::CalculationManager(Grid3DPtr grid, int numOfThreads, double endTime, CommunicatorPtr comm, int endDir, std::shared_ptr<CalculatorFactory> calculatorFactory)
    : grid(grid),
    numOfThreads(numOfThreads),
    endTime(endTime),
    calculatorFactory(calculatorFactory),
    type(type)
{
    this->initCalcThreads();
    loadBalancer = LoadBalancerPtr(new LoadBalancer(grid, comm, endDir));
}
//////////////////////////////////////////////////////////////////////////
CalculationManager::~CalculationManager()
{

}

//////////////////////////////////////////////////////////////////////////
void CalculationManager::calculate()
{
    if (type == CalculatorType::MPI)
    {
        try
        {
            std::dynamic_pointer_cast<MPICalculator>(calcThreads[0])->calculate(endTime, shared_from_this());
        }
        catch (std::exception& e)
        {
            UBLOG(logERROR, e.what());
            //throw e;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        try
        {
            boost::thread_group threads;
            boost::exception_ptr error;

            for (int i = 1; i < calcThreads.size(); i++)
                threads.create_thread(boost::bind(&Calculator::calculate, calcThreads[i], endTime, shared_from_this(), boost::ref(error)));

            calcThreads[0]->calculate(endTime, shared_from_this(), boost::ref(error));

            threads.join_all();
        }
        catch (std::exception& e)
        {
            UBLOG(logERROR, e.what());
            //throw e;
            exit(EXIT_FAILURE);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void CalculationManager::initCalcThreads()
{
    UBLOG(logDEBUG1, "CalculationManager::initCalcThreads() - started");

    SynchronizerPtr sync(new Synchronizer(numOfThreads));
    for (int i = 0; i < numOfThreads; i++)
        calcThreads.push_back(this->calculatorFactory->makeCalculator(grid, sync, i == 0 ? true : false, type));

    this->addBlocksToCalcThreads();
}

//////////////////////////////////////////////////////////////////////////
bool CalculationManager::balance()
{
    //if(loadBalancer->balance())
    //   return true;
    //else
    //{
    //   grid->deleteConnectors();

    //   for(CalculatorPtr calc,  calcThreads)
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
    for (std::shared_ptr<Calculator> c : calcThreads)
        c->deleteBlocks();

    addBlocksToCalcThreads();

    for (std::shared_ptr<Calculator> c : calcThreads)
        c->initConnectors();
}
//////////////////////////////////////////////////////////////////////////
void CalculationManager::addBlocksToCalcThreads()
{
    const int gridRank = grid->getRank();
    const int minInitLevel = this->grid->getCoarsestInitializedLevel();
    const int maxInitLevel = this->grid->getFinestInitializedLevel();

    for (int level = minInitLevel; level <= maxInitLevel; level++)
    {
        std::vector<Block3DPtr> blockVector;
        grid->getBlocks(level, gridRank, true, blockVector);
        for (Block3DPtr const block : blockVector)
            if (block)
                calcThreads[block->getPart()]->addBlock(block);
    }
}

//////////////////////////////////////////////////////////////////////////
void CalculationManager::setTimeAveragedValuesCoProcessor(TimeAveragedValuesCoProcessorPtr coProcessor)
{
    calcThreads[0]->setTimeAveragedValuesCoProcessor(coProcessor);
}

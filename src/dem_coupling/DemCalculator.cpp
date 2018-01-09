#include "DemCalculator.h"


#include "Grid3D.h"
#include "Synchronizer.h"
#include "CalculationManager.h"
#include "UbLogger.h"
#include "Communicator.h"
#include "TimeAveragedValuesCoProcessor.h"
#include "UbScheduler.h"


DemCalculator::DemCalculator(Grid3DPtr grid, SynchronizerPtr sync, bool mainThread) : Calculator(grid, sync, mainThread)
{

}

void DemCalculator::calculate(const double& endTime, CalculationManagerPtr cm, boost::exception_ptr& error)
{
    UBLOG(logINFO, "Calculator PE::calculate() - started");
    try
    {
        initConnectors();

        const int minInitLevel = minLevel;
        const int maxInitLevel = maxLevel - minLevel;
        int straightStartLevel = minInitLevel;
        const int internalIterations = 1 << (maxInitLevel - minInitLevel);
        int threshold;
        const int startStep = int(grid->getTimeStep()) + 1;
        const int anzCalcSteps = static_cast<int>(endTime);

        for (calcStep = startStep; calcStep <= anzCalcSteps + 1; calcStep++)
        {
            sync->wait();
            if (mainThread)
                grid->coProcess((double)(calcStep - 1));
            sync->wait();

            //new part from this calculator:
            exchangeBlockData(straightStartLevel, maxInitLevel);

            for (int staggeredStep = 1; staggeredStep <= internalIterations; staggeredStep++)
            {
                if (staggeredStep == internalIterations)
                    straightStartLevel = minInitLevel;
                else
                {
                    for (straightStartLevel = maxInitLevel, threshold = 1;
                        (staggeredStep & threshold) != threshold; straightStartLevel--, threshold <<= 1);
                }

                calculateBlocks(straightStartLevel, maxInitLevel);

                exchangeBlockData(straightStartLevel, maxInitLevel);

                applyPostCollisionBC(straightStartLevel, maxInitLevel);

                swapDistributions(straightStartLevel, maxInitLevel);

                sync->wait();
                if (taValuesCoProcessor && mainThread)
                    taValuesCoProcessor->calculateSubtotal(calcStep - 1);
                sync->wait();
            }

            //exchange data between blocks for visualization
            if (mainThread) visScheduler->isDue((double)(calcStep - 1));
            if ((int)visScheduler->getNextDueTime() == calcStep)
            {
                exchangeBlockData(straightStartLevel, maxInitLevel);
            }

        }
        UBLOG(logDEBUG1, "Calculator::calculate() - stoped");
    }
    catch (std::exception& e)
    {
        UBLOG(logERROR, e.what());
        UBLOG(logERROR, " step = " << calcStep);
        Communicator::getInstance()->abort(1);
        exit(EXIT_FAILURE);
    }
}

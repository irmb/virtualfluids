#include "DemCalculator.h"


#include "Grid3D.h"
//#include "Synchronizer.h"
#include "UbLogger.h"
#include "Communicator.h"
#include "TimeAveragedValuesCoProcessor.h"
#include "UbScheduler.h"


DemCalculator::DemCalculator(SPtr<Grid3D> grid, SPtr<UbScheduler> additionalGhostLayerUpdateScheduler, int numberOfTimeSteps) :
   BasicCalculator(grid, additionalGhostLayerUpdateScheduler, numberOfTimeSteps)
{

}

void DemCalculator::calculate()
{
    UBLOG(logDEBUG1, "OMPCalculator::calculate() - started");
    int calcStep = 0;
    try
    {
        int minInitLevel = minLevel;
        int maxInitLevel = maxLevel - minLevel;
        int straightStartLevel = minInitLevel;
        int internalIterations = 1 << (maxInitLevel - minInitLevel);
        int forwardStartLevel;
        int threshold;

#ifdef TIMING
        UbTimer timer;
        double time[6];
#endif

        for (calcStep = startTimeStep; calcStep <= numberOfTimeSteps + 1; calcStep++)
        {
            coProcess((double)(calcStep - 1));

            exchangeBlockData(straightStartLevel, maxInitLevel);


            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            UBLOG(logINFO, "calcStep = " << calcStep);
#endif
            //////////////////////////////////////////////////////////////////////////

            for (int staggeredStep = 1; staggeredStep <= internalIterations; staggeredStep++)
            {
                forwardStartLevel = straightStartLevel;
                if (staggeredStep == internalIterations) straightStartLevel = minInitLevel;
                else
                {
                    for (straightStartLevel = maxInitLevel, threshold = 1;
                        (staggeredStep&threshold) != threshold; straightStartLevel--, threshold <<= 1);
                }
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                timer.resetAndStart();
#endif
                //////////////////////////////////////////////////////////////////////////
                applyPreCollisionBC(straightStartLevel, maxInitLevel);

                //do collision for all blocks
                calculateBlocks(straightStartLevel, maxInitLevel, calcStep);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[0] = timer.stop();
                UBLOG(logINFO, "calculateBlocks time = " << time[0]);
#endif
                //////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////
                //exchange data between blocks
                exchangeBlockData(straightStartLevel, maxInitLevel);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[1] = timer.stop();
                UBLOG(logINFO, "exchangeBlockData time = " << time[1]);
#endif
                //////////////////////////////////////////////////////////////////////////
                applyPostCollisionBC(straightStartLevel, maxInitLevel);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[2] = timer.stop();
                UBLOG(logINFO, "applyBCs time = " << time[2]);
#endif
                //////////////////////////////////////////////////////////////////////////
                //swap distributions in kernel
                swapDistributions(straightStartLevel, maxInitLevel);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[3] = timer.stop();
                UBLOG(logINFO, "swapDistributions time = " << time[3]);
#endif
                //////////////////////////////////////////////////////////////////////////
                if (refinement)
                {
                    if (straightStartLevel < maxInitLevel)
                        exchangeBlockData(straightStartLevel, maxInitLevel);
                    //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                    time[4] = timer.stop();
                    UBLOG(logINFO, "refinement exchangeBlockData time = " << time[4]);
#endif
                    //////////////////////////////////////////////////////////////////////////
                    //now ghost nodes have actual values
                    //interpolation of interface nodes between grid levels
                    interpolation(straightStartLevel, maxInitLevel);
                    //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                    time[5] = timer.stop();
                    UBLOG(logINFO, "refinement interpolation time = " << time[5]);
#endif
                    //////////////////////////////////////////////////////////////////////////
                }
            }
            //exchange data between blocks for visualization
            if ((int)additionalGhostLayerUpdateScheduler->getNextDueTime() == calcStep)
            {
                exchangeBlockData(straightStartLevel, maxInitLevel);
            }
            //now ghost nodes have actual values
        }
        UBLOG(logDEBUG1, "OMPCalculator::calculate() - stoped");
    }
    catch (std::exception& e)
    {
        UBLOG(logERROR, e.what());
        UBLOG(logERROR, " step = " << calcStep);
        //throw;
        exit(EXIT_FAILURE);
    }
}

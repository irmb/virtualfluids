#include "PrePostBcCalculator.h"
#include <basics/utilities/UbException.h>
#include <boost/foreach.hpp>
#include "MathUtil.hpp"
#include "basics/writer/WbWriterVtkXmlASCII.h"

//#define TIMING

PrePostBcCalculator::PrePostBcCalculator()
{

}
//////////////////////////////////////////////////////////////////////////
PrePostBcCalculator::PrePostBcCalculator(Grid3DPtr grid, SynchronizerPtr sync, bool mainThread) : 
Calculator(grid, sync, mainThread)
{

}
//////////////////////////////////////////////////////////////////////////
void PrePostBcCalculator::calculate(const double& endTime, CalculationManagerPtr cm, boost::exception_ptr& error)
{
   UBLOG(logDEBUG1, "PrePostBcCalculator::calculate() - started");
   try
   {
      initConnectors();

      int anzLevel = maxLevel-minLevel+1;

      int minInitLevel       = minLevel;
      int maxInitLevel       = maxLevel-minLevel;
      int straightStartLevel = minInitLevel;
      int internalIterations = 1 << (maxInitLevel-minInitLevel);
      int forwardStartLevel;
      int threshold;
      int startStep = int(grid->getTimeStep())+1;

      //UBLOG(logINFO, "startStep="<<startStep);
      int anzCalcSteps = static_cast<int>(endTime);
#ifdef TIMING
      UbTimer timer;
      double time[6];
#endif

//////////////////////////////////////////////////////////////////////////
//      UBLOG(logINFO, "Number of connectors = " <<this->localConns[0].size());
//////////////////////////////////////////////////////////////////////////

      for(calcStep=startStep; calcStep<=anzCalcSteps+1; calcStep++)
      {

         //exchange data between blocks for visualization
         //sync->wait();
         ////if(visScheduler->isDue((double)(calcStep-1)))
         ////{
         //   //exchangeBlockData(minInitLevel, maxInitLevel, true);
         ////}

         ////wait for write dump files
         //sync->wait();
         //write dump 
         if (mainThread) grid->coProcess((double)(calcStep-1));
         sync->wait();


//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
         //UBLOG(logINFO, "calcStep = " <<calcStep);
#endif
//////////////////////////////////////////////////////////////////////////
         
         for(int staggeredStep=1; staggeredStep<=internalIterations; staggeredStep++)
         {
            forwardStartLevel = straightStartLevel;
            if(staggeredStep == internalIterations) straightStartLevel = minInitLevel;
            else
            {
               for(straightStartLevel=maxInitLevel,threshold=1;
                  (staggeredStep&threshold)!=threshold; straightStartLevel--,threshold<<=1);
            }
//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            timer.resetAndStart();
#endif
//////////////////////////////////////////////////////////////////////////

            //applyPreCollisionBC(straightStartLevel, maxInitLevel);


            calculateBlocks(straightStartLevel, maxInitLevel);
            ////calculateBlocks(minInitLevel, maxInitLevel, staggeredStep);
//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[0] = timer.stop();
            //UBLOG(logINFO, "calculateBlocks time = " <<time);
#endif
//////////////////////////////////////////////////////////////////////////

            //exchange data between blocks
            //Sleep(10000);
            exchangeBlockData(straightStartLevel, maxInitLevel);
//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[1] = timer.stop();
            //UBLOG(logINFO, "exchangeBlockData time = " <<time);
#endif
//////////////////////////////////////////////////////////////////////////
            //applyBCs(straightStartLevel, maxInitLevel);
            applyPostCollisionBC(straightStartLevel, maxInitLevel);
            
//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[2] = timer.stop();
            //UBLOG(logINFO, "applyBCs time = " <<time);
#endif
//////////////////////////////////////////////////////////////////////////

            //swap distributions in kernel
            swapDistributions(straightStartLevel, maxInitLevel);

            //pre-collision boundary conditions
            exchangeBlockData(straightStartLevel, maxInitLevel);
            applyPreCollisionBC(straightStartLevel, maxInitLevel);

//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[3] = timer.stop();
            //UBLOG(logINFO, "swapDistributions time = " <<time);
#endif
//////////////////////////////////////////////////////////////////////////

            if (refinement)
            {
         //      //exchange data between blocks for grid refinement
			      ////exchangeInterfaceBlockData(straightStartLevel, maxInitLevel, true);
         //DOES NOT NEED 
               //      if(straightStartLevel<maxInitLevel)
         //         exchangeBlockData(straightStartLevel, maxInitLevel, true);
         //         //exchangeInterfaceBlockData(straightStartLevel, maxInitLevel, true);
//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
               time[4] = timer.stop();
               UBLOG(logINFO, "refinement exchangeBlockData time = " <<time);
#endif
//////////////////////////////////////////////////////////////////////////
               //now ghost nodes have actual values
			      //interpolation of interface nodes between grid levels
               interpolation(straightStartLevel, maxInitLevel);
//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
               time[5] = timer.stop();
               UBLOG(logINFO, "refinement interpolation time = " <<time);
#endif
//////////////////////////////////////////////////////////////////////////
            }
            
         }
         //exchange data between blocks for visualization
        if(mainThread) visScheduler->isDue((double)(calcStep-1));
        if((int)visScheduler->getNextDueTime() == calcStep)
        {
            exchangeBlockData(straightStartLevel, maxInitLevel);
        }
         //now ghost nodes have actual values

         //dynamic load balancing
         //sync->wait();
         //if (mainThread && !loadBalancingComp)
         //{
         //   loadBalancingComp = cm->balance();
         //}
      }
      error = boost::exception_ptr();
      UBLOG(logDEBUG1, "PrePostBcCalculator::calculate() - stoped");
   }
   catch( std::exception& e )
   {
      //error = boost::current_exception();
      UBLOG(logERROR, e.what());
      UBLOG(logERROR, " step = "<<calcStep);
      boost::dynamic_pointer_cast<MPICommunicator>(Communicator::getInstance())->~MPICommunicator();
      exit(EXIT_FAILURE);
   }
}
//////////////////////////////////////////////////////////////////////////
void PrePostBcCalculator::applyPreCollisionBC(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      BOOST_FOREACH(Block3DPtr block, blocks[level])
      {
         block->getKernel()->getBCProcessor()->applyPreCollisionBC();
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void PrePostBcCalculator::applyPostCollisionBC(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      BOOST_FOREACH(Block3DPtr block, blocks[level])
      {
         block->getKernel()->getBCProcessor()->applyPostCollisionBC();
      }
   }
}


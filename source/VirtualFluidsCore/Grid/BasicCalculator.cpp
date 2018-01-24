#include "BasicCalculator.h"

#include "Block3D.h"
#include "BCProcessor.h"
#include "LBMKernel.h"
#include "Block3DConnector.h"
#include "UbScheduler.h"
#include "UbLogger.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#define OMP_SCHEDULE dynamic

//#define TIMING
//#include "UbTiming.h"

BasicCalculator::BasicCalculator(SPtr<Grid3D> grid, SPtr<UbScheduler> additionalGhostLayerUpdateScheduler, int numberOfTimeSteps) : 
   Calculator(grid, additionalGhostLayerUpdateScheduler, numberOfTimeSteps)
{

}
//////////////////////////////////////////////////////////////////////////
BasicCalculator::~BasicCalculator()
{

}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::calculate()
{
   UBLOG(logDEBUG1, "OMPCalculator::calculate() - started");
   int calcStep = 0;
   try
   {
      int minInitLevel = minLevel;
      int maxInitLevel = maxLevel-minLevel;
      int straightStartLevel = minInitLevel;
      int internalIterations = 1<<(maxInitLevel-minInitLevel);
      int forwardStartLevel;
      int threshold;

#ifdef TIMING
      UbTimer timer;
      double time[6];
#endif

      for (calcStep = startTimeStep; calcStep<=numberOfTimeSteps+1; calcStep++)
      {
         coProcess((double)(calcStep-1));

         //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
         UBLOG(logINFO, "calcStep = "<<calcStep);
#endif
         //////////////////////////////////////////////////////////////////////////

         for (int staggeredStep = 1; staggeredStep<=internalIterations; staggeredStep++)
         {
            forwardStartLevel = straightStartLevel;
            if (staggeredStep==internalIterations) straightStartLevel = minInitLevel;
            else
            {
               for (straightStartLevel = maxInitLevel, threshold = 1;
                  (staggeredStep&threshold)!=threshold; straightStartLevel--, threshold <<= 1);
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
            UBLOG(logINFO, "calculateBlocks time = "<<time[0]);
#endif
            //////////////////////////////////////////////////////////////////////////
                        //////////////////////////////////////////////////////////////////////////
                        //exchange data between blocks
            exchangeBlockData(straightStartLevel, maxInitLevel);
            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[1] = timer.stop();
            UBLOG(logINFO, "exchangeBlockData time = "<<time[1]);
#endif
            //////////////////////////////////////////////////////////////////////////
            applyPostCollisionBC(straightStartLevel, maxInitLevel);
            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[2] = timer.stop();
            UBLOG(logINFO, "applyBCs time = "<<time[2]);
#endif
            //////////////////////////////////////////////////////////////////////////
            //swap distributions in kernel
            swapDistributions(straightStartLevel, maxInitLevel);
            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[3] = timer.stop();
            UBLOG(logINFO, "swapDistributions time = "<<time[3]);
#endif
            //////////////////////////////////////////////////////////////////////////
            if (refinement)
            {
               if (straightStartLevel<maxInitLevel)
                  exchangeBlockData(straightStartLevel, maxInitLevel);
               //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
               time[4] = timer.stop();
               UBLOG(logINFO, "refinement exchangeBlockData time = "<<time[4]);
#endif
               //////////////////////////////////////////////////////////////////////////
               //now ghost nodes have actual values
               //interpolation of interface nodes between grid levels
               interpolation(straightStartLevel, maxInitLevel);
               //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
               time[5] = timer.stop();
               UBLOG(logINFO, "refinement interpolation time = "<<time[5]);
#endif
               //////////////////////////////////////////////////////////////////////////
            }
         }
         //exchange data between blocks for visualization
         if ((int)additionalGhostLayerUpdateScheduler->getNextDueTime()==calcStep)
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
      UBLOG(logERROR, " step = "<<calcStep);
      //throw;
      exit(EXIT_FAILURE);
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::calculateBlocks(int startLevel, int maxInitLevel, int calcStep)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
   {
      SPtr<Block3D> blockTemp;
      try
      {
         //startLevel bis maxInitLevel
         for (int level = startLevel; level<=maxInitLevel; level++)
         {
            //timer.resetAndStart();
            //call LBM kernel
            int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp for schedule(OMP_SCHEDULE)
#endif
            for (int i =0; i<size; i++)
            {
               blockTemp = blocks[level][i];
               blockTemp->getKernel()->calculate();
            }
            //timer.stop();
            //UBLOG(logINFO, "level = " << level << " blocks = " << blocks[level].size() << " collision time = " << timer.getTotalTime());
         }
      }
      catch (std::exception& e)
      {
         UBLOG(logERROR, e.what());
         UBLOG(logERROR, blockTemp->toString()<<" step = "<<calcStep);
         //throw;
         exit(EXIT_FAILURE);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::exchangeBlockData(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      //connectorsPrepareLocal(localConns[level]);
      connectorsSendLocal(localConns[level]);
      //connectorsReceiveLocal(localConns[level]);

      connectorsPrepareRemote(remoteConns[level]);
      connectorsSendRemote(remoteConns[level]);
      connectorsReceiveRemote(remoteConns[level]);
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::swapDistributions(int startLevel, int maxInitLevel)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
   {
      //startLevel bis maxInitLevel
      for (int level = startLevel; level<=maxInitLevel; level++)
      {
         int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp for schedule(OMP_SCHEDULE)
#endif
         for (int i =0; i<size; i++)
         {
            blocks[level][i]->getKernel()->swapDistributions();
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::connectorsPrepareLocal(std::vector< SPtr<Block3DConnector> >& connectors)
{
   int size = (int)connectors.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
   for (int i =0; i<size; i++)
   {
      connectors[i]->prepareForReceive();
      connectors[i]->prepareForSend();
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::connectorsSendLocal(std::vector< SPtr<Block3DConnector> >& connectors)
{
   int size = (int)connectors.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
   for (int i =0; i<size; i++)
   {
      connectors[i]->fillSendVectors();
      connectors[i]->sendVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::connectorsReceiveLocal(std::vector< SPtr<Block3DConnector> >& connectors)
{
   int size = (int)connectors.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
   for (int i =0; i<size; i++)
   {
      connectors[i]->receiveVectors();
      connectors[i]->distributeReceiveVectors();
   }
}
void BasicCalculator::connectorsPrepareRemote(std::vector< SPtr<Block3DConnector> >& connectors)
{
   int size = (int)connectors.size();
   for (int i =0; i<size; i++)
   {
      connectors[i]->prepareForReceive();
      connectors[i]->prepareForSend();
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::connectorsSendRemote(std::vector< SPtr<Block3DConnector> >& connectors)
{
   int size = (int)connectors.size();
   for (int i =0; i<size; i++)
   {
      connectors[i]->fillSendVectors();
      connectors[i]->sendVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::connectorsReceiveRemote(std::vector< SPtr<Block3DConnector> >& connectors)
{
   int size = (int)connectors.size();
   for (int i =0; i<size; i++)
   {
      connectors[i]->receiveVectors();
      connectors[i]->distributeReceiveVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::interpolation(int startLevel, int maxInitLevel)
{
   for (int level = startLevel; level<maxInitLevel; level++)
   {
      connectorsPrepareLocal(localInterConns[level]);
      connectorsPrepareRemote(remoteInterConns[level]);
   }

   for (int level = startLevel; level<maxInitLevel; level++)
   {
      connectorsSendLocal(localInterConns[level]);
      connectorsSendRemote(remoteInterConns[level]);
   }

   for (int level = startLevel; level<maxInitLevel; level++)
   {
      connectorsReceiveLocal(localInterConns[level]);
      connectorsReceiveRemote(remoteInterConns[level]);
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::applyPreCollisionBC(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
      for (int i =0; i<size; i++)
      {
         blocks[level][i]->getKernel()->getBCProcessor()->applyPreCollisionBC();
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::applyPostCollisionBC(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
      for (int i =0; i<size; i++)
      {
         blocks[level][i]->getKernel()->getBCProcessor()->applyPostCollisionBC();
      }
   }
}


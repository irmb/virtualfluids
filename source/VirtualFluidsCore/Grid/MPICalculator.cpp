#include "MPICalculator.h"

#include "Grid3D.h"
#include "Block3D.h"
#include "BCProcessor.h"
#include "LBMKernel.h"
#include "Block3DConnector.h"

#include "MathUtil.hpp"
#include "UbScheduler.h"
#include "UbTiming.h"
#include <UbException.h>


//#define TIMING

MPICalculator::MPICalculator()
{

}
//////////////////////////////////////////////////////////////////////////
MPICalculator::~MPICalculator()
{

}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::calculate()
{
   UBLOG(logDEBUG1, "MPICalculator::calculate() - started");
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

      for (calcStep = startTimeStep; calcStep<=lastTimeStep+1; calcStep++)
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
         if ((int)visScheduler->getNextDueTime()==calcStep)
         {
            exchangeBlockData(straightStartLevel, maxInitLevel);
         }
         //now ghost nodes have actual values
      }
      UBLOG(logDEBUG1, "MPICalculator::calculate() - stoped");
   }
   catch (std::exception& e)
   {
      UBLOG(logERROR, e.what());
      UBLOG(logERROR, " step = "<<calcStep);
   }
}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::calculateBlocks(int startLevel, int maxInitLevel, int calcStep)
{
   Block3DPtr blockTemp;
   try
   {
      //startLevel bis maxInitLevel
      for (int level = startLevel; level<=maxInitLevel; level++)
      {
         //timer.resetAndStart();
         //call LBM kernel
         for (Block3DPtr block : blocks[level])
         {
            blockTemp = block;
            block->getKernel()->calculate();
         }
         //timer.stop();
         //UBLOG(logINFO, "level = " << level << " blocks = " << blocks[level].size() << " collision time = " << timer.getTotalTime());
      }
   }
   catch (std::exception& e)
   {
      UBLOG(logERROR, e.what());
      UBLOG(logERROR, blockTemp->toString()<<" step = "<<calcStep);
      exit(EXIT_FAILURE);
   }
}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::exchangeBlockData(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      connectorsPrepare(localConns[level]);
      connectorsPrepare(remoteConns[level]);

      connectorsSend(localConns[level]);
      connectorsSend(remoteConns[level]);

      connectorsReceive(localConns[level]);
      connectorsReceive(remoteConns[level]);
   }
}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::swapDistributions(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      for (Block3DPtr block : blocks[level])
      {
         block->getKernel()->swapDistributions();
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::connectorsPrepare(std::vector< Block3DConnectorPtr >& connectors)
{
   for (Block3DConnectorPtr c : connectors)
   {
      c->prepareForReceive();
      c->prepareForSend();
   }
}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::connectorsSend(std::vector< Block3DConnectorPtr >& connectors)
{
   for (Block3DConnectorPtr c : connectors)
   {
      c->fillSendVectors();
      c->sendVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::connectorsReceive(std::vector< Block3DConnectorPtr >& connectors)
{
   for (Block3DConnectorPtr c : connectors)
   {
      c->receiveVectors();
      c->distributeReceiveVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::interpolation(int startLevel, int maxInitLevel)
{


   for (int level = startLevel; level<maxInitLevel; level++)
   {
      connectorsPrepare(localInterConns[level]);
      connectorsPrepare(remoteInterConns[level]);
   }

   for (int level = startLevel; level<maxInitLevel; level++)
   {
      connectorsSend(localInterConns[level]);
      connectorsSend(remoteInterConns[level]);
   }

   for (int level = startLevel; level<maxInitLevel; level++)
   {
      connectorsReceive(localInterConns[level]);
      connectorsReceive(remoteInterConns[level]);
   }
}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::applyPreCollisionBC(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      for (Block3DPtr block : blocks[level])
      {
         block->getKernel()->getBCProcessor()->applyPreCollisionBC();
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void MPICalculator::applyPostCollisionBC(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      for (Block3DPtr block : blocks[level])
      {
         block->getKernel()->getBCProcessor()->applyPostCollisionBC();
      }
   }
}


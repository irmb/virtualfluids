#include "OMPCalculator.h"
#include <basics/utilities/UbException.h>

#include "MathUtil.hpp"
#include "basics/writer/WbWriterVtkXmlASCII.h"

#include "BCProcessor.h"
#include "LBMKernel.h"

#ifdef _OPENMP
   #include <omp.h>
#endif
//#define TIMING

#include "Block3DConnector.h"

OMPCalculator::OMPCalculator()
{

}
//////////////////////////////////////////////////////////////////////////
OMPCalculator::OMPCalculator(Grid3DPtr grid) :
Calculator(grid, SynchronizerPtr(), true)
{

}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::calculate(const double& endTime, CalculationManagerPtr cm)
{
   UBLOG(logDEBUG1, "MPICalculator::calculate() - started");
   try
   {      
      initConnectors();

      int anzLevel = maxLevel-minLevel+1;

      int minInitLevel = minLevel;
      int maxInitLevel = maxLevel-minLevel;
      int straightStartLevel = minInitLevel;
      int internalIterations = 1<<(maxInitLevel-minInitLevel);
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

      for (calcStep = startStep; calcStep<=anzCalcSteps+1; calcStep++)
      {
         grid->coProcess((double)(calcStep-1));

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

            calculateBlocks(straightStartLevel, maxInitLevel);
//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[0] = timer.stop();
            UBLOG(logINFO, "calculateBlocks time = "<<time[0]);
#endif
            //////////////////////////////////////////////////////////////////////////
            //exchange data between blocks
            exchangeBlockData(straightStartLevel, maxInitLevel);
            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[1] = timer.stop();
            UBLOG(logINFO, "exchangeBlockData time = "<<time[1]);
#endif
            //////////////////////////////////////////////////////////////////////////
                        //applyBCs(straightStartLevel, maxInitLevel);
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
         if (mainThread) visScheduler->isDue((double)(calcStep-1));
         if ((int)visScheduler->getNextDueTime()==calcStep)
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
      UBLOG(logDEBUG1, "MPICalculator::calculate() - stoped");
   }
   catch (std::exception& e)
   {
      UBLOG(logERROR, e.what());
      UBLOG(logERROR, " step = "<<calcStep);
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::calculateBlocks(int startLevel, int maxInitLevel)
{
#ifdef _OPENMP
   #pragma omp parallel
#endif
   {
   Block3DPtr blockTemp;
   try
   {
      //startLevel bis maxInitLevel
      for (int level = startLevel; level<=maxInitLevel; level++)
      {
         //timer.resetAndStart();
         //call LBM kernel
         int size = (int)blocks[level].size();
#ifdef _OPENMP
   #pragma omp for schedule(dynamic)
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
      //error = boost::current_exception();
      UBLOG(logERROR, e.what());
      UBLOG(logERROR, blockTemp->toString()<<" step = "<<calcStep);
      exit(EXIT_FAILURE);
   }
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::exchangeBlockData(int startLevel, int maxInitLevel)
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
void OMPCalculator::swapDistributions(int startLevel, int maxInitLevel)
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
   #pragma omp for schedule(dynamic)
#endif
      for (int i =0; i<size; i++)
      {
         blocks[level][i]->getKernel()->swapDistributions();
      }
   }
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::connectorsPrepareLocal(std::vector< Block3DConnectorPtr >& connectors)
{
   int size = (int)connectors.size();
#ifdef _OPENMP
   #pragma omp parallel for schedule(dynamic)
#endif
   for (int i =0; i<size; i++)
   {
      connectors[i]->prepareForReceive();
      connectors[i]->prepareForSend();
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::connectorsSendLocal(std::vector< Block3DConnectorPtr >& connectors)
{
   int size = (int)connectors.size();
#ifdef _OPENMP
   #pragma omp parallel for schedule(dynamic)
#endif
   for (int i =0; i<size; i++)
   {
      connectors[i]->fillSendVectors();
      connectors[i]->sendVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::connectorsReceiveLocal(std::vector< Block3DConnectorPtr >& connectors)
{
   int size = (int)connectors.size();
#ifdef _OPENMP
   #pragma omp parallel for schedule(dynamic)
#endif
   for (int i =0; i<size; i++)
   {
      connectors[i]->receiveVectors();
      connectors[i]->distributeReceiveVectors();
   }
}
void OMPCalculator::connectorsPrepareRemote(std::vector< Block3DConnectorPtr >& connectors)
{
   int size = (int)connectors.size();
   for (int i =0; i<size; i++)
   {
      connectors[i]->prepareForReceive();
      connectors[i]->prepareForSend();
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::connectorsSendRemote(std::vector< Block3DConnectorPtr >& connectors)
{
   int size = (int)connectors.size();
   for (int i =0; i<size; i++)
   {
      connectors[i]->fillSendVectors();
      connectors[i]->sendVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::connectorsReceiveRemote(std::vector< Block3DConnectorPtr >& connectors)
{
   int size = (int)connectors.size();
   for (int i =0; i<size; i++)
   {
      connectors[i]->receiveVectors();
      connectors[i]->distributeReceiveVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::interpolation(int startLevel, int maxInitLevel)
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
void OMPCalculator::applyPreCollisionBC(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      int size = (int)blocks[level].size();
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (int i =0; i<size; i++)
      {
         blocks[level][i]->getKernel()->getBCProcessor()->applyPreCollisionBC();
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::applyPostCollisionBC(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      int size = (int)blocks[level].size();
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (int i =0; i<size; i++)
      {
         blocks[level][i]->getKernel()->getBCProcessor()->applyPostCollisionBC();
      }
   }
}


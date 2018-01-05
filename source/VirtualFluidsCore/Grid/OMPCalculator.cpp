#include "OMPCalculator.h"
#include <basics/utilities/UbException.h>

#include "MathUtil.hpp"
#include "basics/writer/WbWriterVtkXmlASCII.h"

#include "BCProcessor.h"
#include "LBMKernel.h"
#include <omp.h>
//#define TIMING
//#define PRECOLLISIONBC

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
      omp_set_num_threads(4);
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

         //exchange data between blocks for visualization
         //sync->wait();
         ////if(visScheduler->isDue((double)(calcStep-1)))
         ////{
         //   //exchangeBlockData(minInitLevel, maxInitLevel, true);
         ////}

         ////wait for write dump files
         //write dump 
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
            ////calculateBlocks(minInitLevel, maxInitLevel, staggeredStep);
//////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[0] = timer.stop();
            UBLOG(logINFO, "calculateBlocks time = "<<time[0]);
#endif
            //////////////////////////////////////////////////////////////////////////

                        //exchange data between blocks
                        //Sleep(10000);
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

#ifdef PRECOLLISIONBC
            exchangeBlockData(straightStartLevel, maxInitLevel);
            applyPreCollisionBC(straightStartLevel, maxInitLevel);
#endif


            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[3] = timer.stop();
            UBLOG(logINFO, "swapDistributions time = "<<time[3]);
#endif
            //////////////////////////////////////////////////////////////////////////

            if (refinement)
            {
               //      //exchange data between blocks for grid refinement
                     ////exchangeInterfaceBlockData(straightStartLevel, maxInitLevel, true);
               //DOES NOT NEED 
               if (straightStartLevel<maxInitLevel)
                  exchangeBlockData(straightStartLevel, maxInitLevel);
               //         //exchangeInterfaceBlockData(straightStartLevel, maxInitLevel, true);
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

            if (taValuesCoProcessor)
            {
               taValuesCoProcessor->calculateSubtotal(calcStep-1);
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
      //error = boost::current_exception();
      UBLOG(logERROR, e.what());
      UBLOG(logERROR, " step = "<<calcStep);
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::calculateBlocks(int startLevel, int maxInitLevel)
{
   #pragma omp parallel 
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
         #pragma omp for schedule(dynamic)
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
void OMPCalculator::calculateBlocks(int minInitLevel, int maxInitLevel, int staggeredStep)
{
   int p, maxi, maxir, maxidp, start, end;
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      p = 1<<(maxInitLevel-level);
      maxi = maxir = static_cast<int>(blocks[level].size());
      maxidp = maxi/p;
      if (p>maxi && maxi!=0) {
         maxidp = 1;
         maxi = p;
      }
      start = (staggeredStep-1)*maxidp;
      if (start>=maxi)
         start = 0;
      end = start+maxidp;
      if ((end+p)>=maxi)
         end = maxi;
      for (int i = start; i<end; i++)
      {
         if (i<maxir)
            blocks[level][i]->getKernel()->calculate();
      }
   }
}
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
void OMPCalculator::exchangeBlockData(int startLevel, int maxInitLevel)
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
void OMPCalculator::exchangeInterfaceBlockData(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      connectorsPrepare(localInterfaceBlockConns[level]);
      connectorsPrepare(remoteInterfaceBlockConns[level]);

      connectorsSend(localInterfaceBlockConns[level]);
      connectorsSend(remoteInterfaceBlockConns[level]);

      connectorsReceive(localInterfaceBlockConns[level]);
      connectorsReceive(remoteInterfaceBlockConns[level]);
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::swapDistributions(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      for(Block3DPtr block : blocks[level])
      {
         block->getKernel()->swapDistributions();
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::connectorsPrepare(std::vector< Block3DConnectorPtr >& connectors)
{
   int size = (int)connectors.size();
   #pragma omp parallel for schedule(dynamic)
   for (int i =0; i<size; i++)
   {
      connectors[i]->prepareForReceive();
      connectors[i]->prepareForSend();
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::connectorsSend(std::vector< Block3DConnectorPtr >& connectors)
{
   int size = (int)connectors.size();
   #pragma omp parallel for schedule(dynamic)
   for (int i =0; i<size; i++)
   {
      connectors[i]->fillSendVectors();
      connectors[i]->sendVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::connectorsReceive(std::vector< Block3DConnectorPtr >& connectors)
{
   int size = (int)connectors.size();
   #pragma omp parallel for schedule(dynamic)
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
//////////////////////////////////////////////////////////////////////////
void OMPCalculator::applyPreCollisionBC(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for (int level = startLevel; level<=maxInitLevel; level++)
   {
      int size = (int)blocks[level].size();
      #pragma omp for
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
      #pragma omp parallel for
      for (int i =0; i<size; i++)
      {
         blocks[level][i]->getKernel()->getBCProcessor()->applyPostCollisionBC();
      }
   }
}


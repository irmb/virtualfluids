#if defined VF_FETOL

#if defined(_WIN32) || defined(_WIN64)
   #include <winsock2.h>
#endif

#include "FETOLCalculator.h"
#include <basics/utilities/UbException.h>
#include <D3Q27OffsetInterpolationProcessor.h>
#include <FETOLSetConnectorsBlockVisitor.h>

//problem with Windows, by Unix should be uncomment 
#include <JM.h>
using namespace fetol;

FETOLCalculator::FETOLCalculator()
{

}
//////////////////////////////////////////////////////////////////////////
FETOLCalculator::FETOLCalculator(Grid3DPtr grid, SynchronizerPtr sync, bool mainThread) : 
   Calculator(grid, sync, mainThread)

{
   remoteMPIConns.resize(maxLevel+1);
   remoteBondConns.resize(maxLevel+1);
}
//////////////////////////////////////////////////////////////////////////
void FETOLCalculator::calculate(const double& endTime, CalculationManagerPtr cm, boost::exception_ptr& error)
{
   UBLOG(logDEBUG1, "FETOLCalculator::calculate() - started");
   try
   {
      int anzLevel = maxLevel-minLevel+1;

      int minInitLevel       = minLevel;
      int maxInitLevel       = maxLevel-minLevel;
      int straightStartLevel = minInitLevel;
      int internalIterations = 1 << (maxInitLevel-minInitLevel);
      int forwardStartLevel;
      int threshold;
      int startStep = int(grid->getTimeStep())+1;

      if (startStep > 1)
      {
         vector<int> cs;
         cs.push_back(startStep);
         Communicator::getInstance()->broadcastInts(cs);
         startStep = cs[0];
      }

      //UBLOG(logINFO, "startStep="<<startStep);
      int anzCalcSteps = static_cast<int>(endTime);
#ifdef TIMING
      UbTimer timer;
      double time[6];
#endif

      //////////////////////////////////////////////////////////////////////////
      //      UBLOG(logINFO, "Number of connectors = " <<this->localConns[0].size());
      //////////////////////////////////////////////////////////////////////////

      

      for(calcStep=startStep; calcStep<=anzCalcSteps; calcStep++)
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
         if (mainThread) 
         {
            try
            {
               UBLOG(logDEBUG1,"JM::getApplicationState(): "<<JM::getApplicationState());
               grid->doPostProcess((double)(calcStep-1));
            }
            catch (...)
            {
               ifRestart(straightStartLevel, maxInitLevel, false);
            }
         }
            sync->wait();


         //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
         UBLOG(logINFO, "calcStep = " <<calcStep);
#endif
         UBLOG(logDEBUG1, "FETOLCalculator::calculate() calcStep = " <<calcStep);
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
            calculateBlocks(straightStartLevel, maxInitLevel);
            //calculateBlocks(minInitLevel, maxInitLevel, staggeredStep);
            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[0] = timer.stop();
            //UBLOG(logINFO, "calculateBlocks time = " <<time);
#endif
            //////////////////////////////////////////////////////////////////////////

            //exchange data between blocks
            //Sleep(10000);
            try
            {
               UBLOG(logDEBUG1,"JM::getApplicationState(): "<<JM::getApplicationState());
               exchangeFETOLBlockData(straightStartLevel, maxInitLevel, false);
            }
            catch (...)
            {
               ifRestart(straightStartLevel, maxInitLevel, false);
            }

            UBLOG(logDEBUG1,"JM::getApplicationState(): "<<JM::getApplicationState());

            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[1] = timer.stop();
            //UBLOG(logINFO, "exchangeBlockData time = " <<time);
#endif
            //////////////////////////////////////////////////////////////////////////

            applyBCs(straightStartLevel, maxInitLevel);
            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[2] = timer.stop();
            //UBLOG(logINFO, "applyBCs time = " <<time);
#endif
            //////////////////////////////////////////////////////////////////////////

            //swap distributions in kernel
            swapDistributions(straightStartLevel, maxInitLevel);
            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            time[3] = timer.stop();
            //UBLOG(logINFO, "swapDistributions time = " <<time);
#endif
            //////////////////////////////////////////////////////////////////////////

//            if (refinement)
//            {
//               //exchange data between blocks for grid refinement
//               //exchangeInterfaceBlockData(straightStartLevel, maxInitLevel, true);
//               if(straightStartLevel<maxInitLevel)
//                  exchangeBlockData(straightStartLevel, maxInitLevel, true);
//               //exchangeInterfaceBlockData(straightStartLevel, maxInitLevel, true);
//               //////////////////////////////////////////////////////////////////////////
//#ifdef TIMING
//               time[4] = timer.stop();
//               UBLOG(logINFO, "refinement exchangeBlockData time = " <<time);
//#endif
//               //////////////////////////////////////////////////////////////////////////
//               //now ghost nodes have actual values
//               //interpolation of interface nodes between grid levels
//               interpolation(straightStartLevel, maxInitLevel);
//               //////////////////////////////////////////////////////////////////////////
//#ifdef TIMING
//               time[5] = timer.stop();
//               UBLOG(logINFO, "refinement interpolation time = " <<time);
//#endif
//               //////////////////////////////////////////////////////////////////////////
//            }

         }
         //exchange data between blocks for visualization
         if(mainThread) visScheduler->isDue((double)(calcStep-1));
         if((int)visScheduler->getNextDueTime() == calcStep)
         {
            try
            {
               UBLOG(logDEBUG1,"JM::getApplicationState(): "<<JM::getApplicationState());
               exchangeFETOLBlockData(straightStartLevel, maxInitLevel, true);
            }
            catch (...)
            {
               ifRestart(straightStartLevel, maxInitLevel, false);
            }
         }
         //now ghost nodes have actual values

      }
      error = boost::exception_ptr();
      UBLOG(logDEBUG1, "FETOLCalculator::calculate() - stoped");
   }
   catch( std::exception& e )
   {
      //error = boost::current_exception();
      UBLOG(logERROR, e.what());
      UBLOG(logERROR, " step = "<<calcStep);
      exit(EXIT_FAILURE);
   }
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void FETOLCalculator::initRemoteConnectors()
{
   Calculator::initRemoteConnectors();
   initFETOLConnectors();
}
//////////////////////////////////////////////////////////////////////////
void FETOLCalculator::initFETOLConnectors()
{
   UBLOG(logDEBUG5, "FETOLCalculator::initFETOLConnectors():start");
   for(int l = minLevel; l <= maxLevel; l++)
   {
      for(Block3DConnectorPtr c : remoteConns[l])
      {
         if(c->getTransmitterType() == Block3DConnector::MPI)
         {
            remoteMPIConns[l].push_back(c);
         }
         else if(c->getTransmitterType() == Block3DConnector::BOND)
         {
            remoteBondConns[l].push_back(c);
         }
         else
         {
            UB_THROW( UbException(UB_EXARGS,"Transmitter type isn't exist!"));
         }
      }
   }
   UBLOG(logDEBUG5, "FETOLCalculator::initFETOLConnectors():remoteMPIConns.size = "<<remoteMPIConns[0].size());
   UBLOG(logDEBUG5, "FETOLCalculator::initFETOLConnectors():remoteBondConns.size = "<<remoteBondConns[0].size());
   UBLOG(logDEBUG5, "FETOLCalculator::initFETOLConnectors():end");
}
//////////////////////////////////////////////////////////////////////////
void FETOLCalculator::exchangeFETOLBlockData(int startLevel, int maxInitLevel, bool invStep)
{
   UBLOG(logDEBUG5, "FETOLCalculator::exchangeFETOLBlockData():start");
   sync->wait();
   UBLOG(logDEBUG5, "FETOLCalculator::exchangeFETOLBlockData():localConns:start");
   //startLevel bis maxInitLevel
   for(int level=startLevel; level<=maxInitLevel; level++)
   {
      connectorsPrepare(localConns[level]);

      connectorsSend(localConns[level], invStep);

      connectorsReceive(localConns[level], invStep);
   }
   UBLOG(logDEBUG5, "FETOLCalculator::exchangeFETOLBlockData():localConns:end");
   sync->wait();

   UBLOG(logDEBUG5, "FETOLCalculator::exchangeFETOLBlockData():remoteMPIConns:start");
   for(int level=startLevel; level<=maxInitLevel; level++)
   {
      connectorsPrepare(remoteMPIConns[level]);

      connectorsSend(remoteMPIConns[level], invStep);

      connectorsReceive(remoteMPIConns[level], invStep);
   }
   UBLOG(logDEBUG5, "FETOLCalculator::exchangeFETOLBlockData():remoteMPIConns:end");

   UBLOG(logDEBUG5, "FETOLCalculator::exchangeFETOLBlockData():remoteBondConns:start");
   for(int level=startLevel; level<=maxInitLevel; level++)
   {
      connectorsPrepare(remoteBondConns[level]);

      connectorsSend(remoteBondConns[level], invStep);

      connectorsReceive(remoteBondConns[level], invStep);
   }
   UBLOG(logDEBUG5, "FETOLCalculator::exchangeFETOLBlockData():remoteBondConns:end");
   UBLOG(logDEBUG5, "FETOLCalculator::exchangeFETOLBlockData():end");
}
//////////////////////////////////////////////////////////////////////////
void FETOLCalculator::ifRestart(int startLevel, int maxInitLevel, bool invStep)
{
   UBLOG(logINFO, "before paused");
   while (JM::getApplicationState() != paused)
   {
      ; //NOP
   }

   UBLOG(logINFO,"JM::getApplicationState(): "<<JM::getApplicationState());

   while (JM::getApplicationState() == paused) 
   {
      ; //NOP
   }
   UBLOG(logINFO, "after paused");
   deleteConnectors(remoteConns);
   deleteConnectors(remoteMPIConns);
   deleteConnectors(remoteBondConns);
   grid->deleteConnectors();
   D3Q27InterpolationProcessorPtr iProcessor(new D3Q27OffsetInterpolationProcessor());
   CommunicatorPtr comm = FETOLCommunicator::getInstance();
   double nueLB=1;
   FETOLSetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nueLB, iProcessor);
   grid->accept( setConnsVisitor );
   initRemoteConnectors();
   vector<int> cs;
   cs.push_back(calcStep);
   Communicator::getInstance()->broadcastInts(cs);
   exchangeFETOLBlockData(startLevel, maxInitLevel, false);
   //UBLOG(logINFO,"new exchange JM::getApplicationState(): "<<JM::getApplicationState());
}

//////////////////////////////////////////////////////////////////////////
//void FETOLCalculator::connectorsPrepare(std::vector< Block3DConnectorPtr >& connectors)
//{
//   for(Block3DConnectorPtr c : connectors)
//   {
//      c->prepareForReceive();
//      c->prepareForSend();
//   }
//}
////////////////////////////////////////////////////////////////////////////
//void FETOLCalculator::connectorsSend(std::vector< Block3DConnectorPtr >& connectors, bool invStep)
//{
//   for(Block3DConnectorPtr c : connectors)
//   {
//      c->setInvStep(invStep);
//      c->fillSendVectors();
//      c->sendVectors();
//   }
//}
////////////////////////////////////////////////////////////////////////////
//void FETOLCalculator::connectorsReceive(std::vector< Block3DConnectorPtr >& connectors, bool invStep)
//{
//   for(Block3DConnectorPtr c : connectors)
//   {
//      c->receiveVectors();
//      c->setInvStep(invStep);
//      c->distributeReceiveVectors();
//   }
//}
////////////////////////////////////////////////////////////////////////////
//void FETOLCalculator::bondConnectorsPrepare(std::vector< Block3DConnectorPtr >& connectors, bool invStep)
//{
//   for(Block3DConnectorPtr c : connectors)
//   {
//      c->prepareForReceive();
//      c->prepareForSend();
//   }
//}
////////////////////////////////////////////////////////////////////////////
//void FETOLCalculator::bondConnectorsSend(std::vector< Block3DConnectorPtr >& connectors, bool invStep)
//{
//   for(Block3DConnectorPtr c : connectors)
//   {
//      c->setInvStep(invStep);
//      c->fillSendVectors();
//      c->sendVectors();
//   }
//}
////////////////////////////////////////////////////////////////////////////
//void FETOLCalculator::bondConnectorsReceive(std::vector< Block3DConnectorPtr >& connectors, bool invStep)
//{
//   for(Block3DConnectorPtr c : connectors)
//   {
//      c->receiveVectors();
//      c->setInvStep(invStep);
//      c->distributeReceiveVectors();
//   }
//}

//////////////////////////////////////////////////////////////////////////

#endif

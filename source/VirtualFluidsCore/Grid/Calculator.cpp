#include "Calculator.h"
#include <basics/utilities/UbException.h>
#include <boost/foreach.hpp>
#include "SimulationParameters.h"
#include "MathUtil.hpp"
#include "basics/writer/WbWriterVtkXmlASCII.h"

//#define TIMING

Calculator::Calculator()
{

}
//////////////////////////////////////////////////////////////////////////
Calculator::Calculator(Grid3DPtr grid, SynchronizerPtr sync, bool mainThread) : 
                       grid(grid),
                       sync(sync),
                       mainThread(mainThread),
                       refinement(false)
{
   minLevel = grid->getCoarsestInitializedLevel();
   maxLevel = grid->getFinestInitializedLevel();
   if(maxLevel > 0)
      refinement = true;
   else
      refinement = false;
   blocks.resize(maxLevel+1);
   localConns.resize(maxLevel+1);
   remoteConns.resize(maxLevel+1);
   localInterfaceBlockConns.resize(maxLevel+1);
   remoteInterfaceBlockConns.resize(maxLevel+1);
   localInterConns.resize(maxLevel);
   remoteInterConns.resize(maxLevel);
   loadBalancingComp = false;
}
//////////////////////////////////////////////////////////////////////////
void Calculator::calculate(const double& endTime, CalculationManagerPtr cm, boost::exception_ptr& error)
{
   UBLOG(logDEBUG1, "Calculator::calculate() - started");
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
         sync->wait();
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
            exchangeBlockData(straightStartLevel, maxInitLevel, false);
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

            //exchangeBlockData(straightStartLevel, maxInitLevel, true);
            //applyPreCollisionBC(straightStartLevel, maxInitLevel);
            

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
                     if(straightStartLevel<maxInitLevel)
                        exchangeBlockData(straightStartLevel, maxInitLevel, true);
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

            if (taValuesCoProcessor && mainThread)
            {
               taValuesCoProcessor->process(calcStep-1);
            }
            sync->wait();
            
         }
         //exchange data between blocks for visualization
        if(mainThread) visScheduler->isDue((double)(calcStep-1));
        if((int)visScheduler->getNextDueTime() == calcStep)
        {
            exchangeBlockData(straightStartLevel, maxInitLevel, true);
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
      UBLOG(logDEBUG1, "Calculator::calculate() - stoped");
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
void Calculator::calculateBlocks(int startLevel, int maxInitLevel)
{
   Block3DPtr blockTemp;
   try
   {
      //startLevel bis maxInitLevel
      for(int level=startLevel; level<=maxInitLevel; level++)
      {
         //timer.resetAndStart();
         //call LBM kernel
         BOOST_FOREACH(Block3DPtr block, blocks[level])
         {
            blockTemp = block;
            block->getKernel()->calculate();
         }
         //timer.stop();
         //UBLOG(logINFO, "level = " << level << " blocks = " << blocks[level].size() << " collision time = " << timer.getTotalTime());
      }
   }
   catch( std::exception& e )
   {      
      //error = boost::current_exception();
      UBLOG(logERROR, e.what());
      UBLOG(logERROR, blockTemp->toString()<<" step = "<<calcStep);
      exit(EXIT_FAILURE);
   }
}
//////////////////////////////////////////////////////////////////////////
void Calculator::calculateBlocks(int minInitLevel, int maxInitLevel, int staggeredStep)
{
   int p, maxi, maxir, maxidp, start, end;
   for(int level=minInitLevel; level<=maxInitLevel; level++)
   {
      p = 1<<(maxInitLevel-level);
      maxi = maxir = static_cast<int>(blocks[level].size());
      maxidp = maxi/p;
      if(p > maxi && maxi != 0){
         maxidp = 1;
         maxi = p;
      }
      start = (staggeredStep-1)*maxidp;
      if(start >= maxi)
         start = 0;
      end = start + maxidp;
      if((end + p) >= maxi)
         end = maxi;
      for (int i = start; i < end; i++)
      {
         if(i < maxir)
            blocks[level][i]->getKernel()->calculate();
      }
   }
 }
//////////////////////////////////////////////////////////////////////////
void Calculator::addBlock(Block3DPtr block)
{
   blocks[block->getLevel()].push_back(block);
}
//////////////////////////////////////////////////////////////////////////
void Calculator::initConnectors()
{
   UBLOG(logDEBUG1, "Calculator::initLocalConnectors() - start");

   for (int l = minLevel; l <= maxLevel; l++)
   {
      BOOST_FOREACH(Block3DPtr block, blocks[l])
      {     
         block->pushBackLocalSameLevelConnectors(localConns[l]);

         if(block->hasInterpolationFlag())
            block->pushBackLocalSameLevelConnectors(localInterfaceBlockConns[l]);
         if (l != maxLevel)
            block->pushBackLocalInterpolationConnectorsCF(localInterConns[l]);
      }
      if (l != maxLevel)
      {
         BOOST_FOREACH(Block3DPtr block, blocks[l+1])
         {     
            block->pushBackLocalInterpolationConnectorsFC(localInterConns[l]);
         }
      }
      UBLOG(logDEBUG5, "Calculator::initConnectors()-initConnectors(localConns["<<l<<"])");
      initConnectors(localConns[l]);

      if (l != maxLevel)
      {
         UBLOG(logDEBUG5, "Calculator::initConnectors()-initConnectors(localInterConns["<<l<<"])");
         initConnectors(localInterConns[l]);
      }
   }
   
   if (mainThread)
      initRemoteConnectors();

   UBLOG(logDEBUG1, "Calculator::initLocalConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void Calculator::initRemoteConnectors()
{
   std::vector< std::vector< Block3DConnectorPtr > > remoteInterConnsCF;
   std::vector< std::vector< Block3DConnectorPtr > > remoteInterConnsFC;
   remoteInterConnsCF.resize(maxLevel+1);
   remoteInterConnsFC.resize(maxLevel+1);

   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();
   int gridRank = grid->getRank();

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      std::vector<Block3DPtr> blockVector;
      //grid->getBlocks(level, gridRank, true, blockVector);
      grid->getBlocks(level, blockVector);
      BOOST_FOREACH(Block3DPtr block, blockVector)
      {
         int l = block->getLevel();
         block->pushBackRemoteSameLevelConnectors(remoteConns[l]);

         //if(block->isInterface())
         //   block->pushBackRemoteSameLevelConnectors(remoteInterfaceBlockConns[l]);
         block->pushBackRemoteInterpolationConnectorsCF(remoteInterConnsCF[l]);
         block->pushBackRemoteInterpolationConnectorsFC(remoteInterConnsFC[l]);
      }
   }

   for (int l = minLevel; l <= maxLevel; l++)
   {
      UBLOG(logDEBUG5, "Calculator::initRemoteConnectors()-initConnectors(remoteConns["<<l<<"])");
      initConnectors(remoteConns[l]);
      if (l != maxLevel)
      {
		 for(int i = 0; i < remoteInterConnsCF[l].size(); i++)
			remoteInterConns[l].push_back(remoteInterConnsCF[l][i]);
		 for(int i = 0; i < remoteInterConnsFC[l+1].size(); i++)
	      remoteInterConns[l].push_back(remoteInterConnsFC[l+1][i]);
       //UBLOG(logDEBUG5, "Calculator::initRemoteConnectors()-initConnectors(remoteInterConns["<<l<<"])");
       //initConnectors(remoteInterConns[l]);
      }
   }
   //////////////////////////////////////////////////////////////////////////
   //UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - start");
   for (int l = minLevel; l <= maxLevel; l++)
   {
      if (l != maxLevel)
      {
         UBLOG(logDEBUG5, "Calculator::initRemoteConnectors()-initConnectors(remoteInterConns["<<l<<"])");
         BOOST_FOREACH(Block3DConnectorPtr c, remoteInterConns[l] ) c->init();
      }
   }
   //UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - end");
   //////////////////////////////////////////////////////////////////////////
   //sendTransmitterDataSize
   //UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - start");
   for (int l = minLevel; l <= maxLevel; l++)
   {
      if (l != maxLevel)
      {
         UBLOG(logDEBUG5, "Calculator::initRemoteConnectors()-sendTransmitterDataSize(remoteInterConns["<<l<<"])");
         BOOST_FOREACH(Block3DConnectorPtr c, remoteInterConns[l] ) c->sendTransmitterDataSize();
      }
   }
   //UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - end");
   //////////////////////////////////////////////////////////////////////////
   //receiveTransmitterDataSize
   //wenn er hier bei verteilten berechnungen stopped, dann ist vermutlich auf einer seite ein nicht aktiver block!!!
   //UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - start");
   for (int l = minLevel; l <= maxLevel; l++)
   {
      if (l != maxLevel)
      {
         UBLOG(logDEBUG5, "Calculator::initRemoteConnectors()-receiveTransmitterDataSize(remoteInterConns["<<l<<"])");
         BOOST_FOREACH(Block3DConnectorPtr c, remoteInterConns[l] ) c->receiveTransmitterDataSize();
      }
   }
   //UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - end");
   //////////////////////////////////////////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////
void Calculator::initConnectors(std::vector<Block3DConnectorPtr>& connectors)
{
   UBLOG(logDEBUG1, "Calculator::initConnectors() - start");

   //initialization
   //////////////////////////////////////////////////////////////////////////
   //initialize connectors
   UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - start");
   BOOST_FOREACH(Block3DConnectorPtr c, connectors ) c->init();
   UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - end");
   //////////////////////////////////////////////////////////////////////////
   //sendTransmitterDataSize
   UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - start");
   BOOST_FOREACH(Block3DConnectorPtr c, connectors ) c->sendTransmitterDataSize();
   UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - end");
   //////////////////////////////////////////////////////////////////////////
   //receiveTransmitterDataSize
   //wenn er hier bei verteilten berechnungen stopped, dann ist vermutlich auf einer seite ein nicht aktiver block!!!
   UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - start");
   BOOST_FOREACH(Block3DConnectorPtr c, connectors ) c->receiveTransmitterDataSize();
   UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - end");

   UBLOG(logDEBUG1, "Calculator::initConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void Calculator::exchangeBlockData(int startLevel, int maxInitLevel, bool invStep)
{
   sync->wait();
   //startLevel bis maxInitLevel
   for(int level=startLevel; level<=maxInitLevel; level++)
   {
      connectorsPrepare(localConns[level]);
      connectorsPrepare(remoteConns[level]);

      connectorsSend(localConns[level], invStep);
      connectorsSend(remoteConns[level], invStep);

      connectorsReceive(localConns[level], invStep);
      connectorsReceive(remoteConns[level], invStep);
   }
   sync->wait();
}
//////////////////////////////////////////////////////////////////////////
void Calculator::exchangeInterfaceBlockData(int startLevel, int maxInitLevel, bool invStep)
{
   sync->wait();
   //startLevel bis maxInitLevel
   for(int level=startLevel; level<=maxInitLevel; level++)
   {
      connectorsPrepare(localInterfaceBlockConns[level]);
      connectorsPrepare(remoteInterfaceBlockConns[level]);

      connectorsSend(localInterfaceBlockConns[level], invStep);
      connectorsSend(remoteInterfaceBlockConns[level], invStep);

      connectorsReceive(localInterfaceBlockConns[level], invStep);
      connectorsReceive(remoteInterfaceBlockConns[level], invStep);
   }
   sync->wait();
}
//////////////////////////////////////////////////////////////////////////
void Calculator::swapDistributions(int startLevel, int maxInitLevel)
{
   //startLevel bis maxInitLevel
   for(int level=startLevel; level<=maxInitLevel; level++)
   {
      BOOST_FOREACH(Block3DPtr block, blocks[level])
      {
         block->getKernel()->swapDistributions();
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void Calculator::connectorsPrepare(std::vector< Block3DConnectorPtr >& connectors)
{
   BOOST_FOREACH(Block3DConnectorPtr c, connectors)
   {
      c->prepareForReceive();
      c->prepareForSend();
   }
}
//////////////////////////////////////////////////////////////////////////
void Calculator::connectorsSend(std::vector< Block3DConnectorPtr >& connectors, bool invStep)
{
   BOOST_FOREACH(Block3DConnectorPtr c, connectors)
   {
	   c->setInvStep(invStep);
      c->fillSendVectors();
      c->sendVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void Calculator::connectorsReceive(std::vector< Block3DConnectorPtr >& connectors, bool invStep)
{
   BOOST_FOREACH(Block3DConnectorPtr c, connectors)
   {
	   c->setInvStep(invStep);
      c->receiveVectors();
      c->distributeReceiveVectors();
   }
}
//////////////////////////////////////////////////////////////////////////
void Calculator::connectorsSetInvStep(std::vector< Block3DConnectorPtr >& connectors, bool invStep)
{
   BOOST_FOREACH(Block3DConnectorPtr c, connectors)
   {
      c->setInvStep(invStep);
   }
}
//////////////////////////////////////////////////////////////////////////
void Calculator::interpolation(int startLevel, int maxInitLevel)
{
   sync->wait();

   for(int level=startLevel; level<maxInitLevel; level++)
   {
      connectorsPrepare(localInterConns[level]);
      connectorsPrepare(remoteInterConns[level]);
   }

   sync->wait();

   for(int level=startLevel; level<maxInitLevel; level++)
   {
      connectorsSend(localInterConns[level], true);
      connectorsSend(remoteInterConns[level], true);
   }

   sync->wait();

   for(int level=startLevel; level<maxInitLevel; level++)
   {
      connectorsReceive(localInterConns[level], true);
      connectorsReceive(remoteInterConns[level], true);
   }

   sync->wait();

}
//////////////////////////////////////////////////////////////////////////
void Calculator::setVisScheduler(UbSchedulerPtr s)
{
   visScheduler = s;
}
//////////////////////////////////////////////////////////////////////////
//double Calculator::getCallculationTime()
//{
//   return timer.getTotalTime();
//}
//////////////////////////////////////////////////////////////////////////
std::vector< std::vector< Block3DPtr > > Calculator::getBlocks()
{
   return blocks;
}
//////////////////////////////////////////////////////////////////////////
void Calculator::deleteBlocks()
{
   BOOST_FOREACH(std::vector< Block3DPtr > &bs, blocks)
      bs.resize(0);
}
//////////////////////////////////////////////////////////////////////////
void Calculator::deleteConnectors()
{
   deleteConnectors(localConns);
   deleteConnectors(remoteConns);

   deleteConnectors(localInterfaceBlockConns);
   deleteConnectors(remoteInterfaceBlockConns);

   deleteConnectors(localInterConns);
   deleteConnectors(remoteInterConns);
}
//////////////////////////////////////////////////////////////////////////
void Calculator::deleteConnectors(std::vector< std::vector< Block3DConnectorPtr > >& conns)
{
   BOOST_FOREACH(std::vector< Block3DConnectorPtr > &c, conns)
      c.resize(0);
}
//////////////////////////////////////////////////////////////////////////
void Calculator::setTimeAveragedValuesCoProcessor(TimeAveragedValuesCoProcessorPtr coProcessor)
{
   taValuesCoProcessor = coProcessor;
}

//////////////////////////////////////////////////////////////////////////
//void Calculator::applyBCs( int startLevel, int maxInitLevel )
//{
//   //startLevel bis maxInitLevel
//   for(int level=startLevel; level<=maxInitLevel; level++)
//   {
//      //call LBM kernel
//      BOOST_FOREACH(Block3DPtr block, blocks[level])
//      {
//         block->getKernel()->getBCProcessor()->applyBC();
//      }
//   }
//}
//////////////////////////////////////////////////////////////////////////
void Calculator::applyPreCollisionBC(int startLevel, int maxInitLevel)
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
void Calculator::applyPostCollisionBC(int startLevel, int maxInitLevel)
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


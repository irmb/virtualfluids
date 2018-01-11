#include "Calculator.h"

#include "Grid3D.h"
#include "Block3DConnector.h"
#include "Block3D.h"
#include "UbScheduler.h"
#include "CoProcessor.h"

#include <basics/utilities/UbException.h>

Calculator::Calculator()
{

}
//////////////////////////////////////////////////////////////////////////
Calculator::~Calculator()
{

}
//////////////////////////////////////////////////////////////////////////
void Calculator::setGrid(std::shared_ptr<Grid3D> grid)
{
   this->grid = grid;
   startTimeStep = int(grid->getTimeStep())+1;
   minLevel = grid->getCoarsestInitializedLevel();
   maxLevel = grid->getFinestInitializedLevel();
   if (maxLevel > 0)
      refinement = true;
   else
      refinement = false;
   blocks.resize(maxLevel+1);
   localConns.resize(maxLevel+1);
   remoteConns.resize(maxLevel+1);
   localInterConns.resize(maxLevel);
   remoteInterConns.resize(maxLevel);

   int gridRank = grid->getRank();

   for (int level = minLevel; level <= maxLevel; level++)
   {
      std::vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, true, blockVector);
      for (Block3DPtr const block : blockVector)
         if (block)
            blocks[block->getLevel()].push_back(block);
   }

   initLocalConnectors();
   initRemoteConnectors();
}
//////////////////////////////////////////////////////////////////////////
void Calculator::setLastTimeStep(int t)
{
   this->lastTimeStep = t;
}
//////////////////////////////////////////////////////////////////////////
void Calculator::setVisScheduler(std::shared_ptr<UbScheduler> s)
{
   visScheduler = s;
}
//////////////////////////////////////////////////////////////////////////
void Calculator::addCoProcessor(std::shared_ptr<CoProcessor> coProcessor)
{
   coProcessors.push_back(coProcessor);
}
//////////////////////////////////////////////////////////////////////////
void Calculator::coProcess(double step)
{
   for (CoProcessorPtr cp : coProcessors)
   {
      cp->process(step);
   }
}
//////////////////////////////////////////////////////////////////////////
void Calculator::initLocalConnectors()
{
   UBLOG(logDEBUG1, "Calculator::initLocalConnectors() - start");

   for (int l = minLevel; l <= maxLevel; l++)
   {
      for(Block3DPtr block : blocks[l])
      {     
         block->pushBackLocalSameLevelConnectors(localConns[l]);

         if (l != maxLevel)
            block->pushBackLocalInterpolationConnectorsCF(localInterConns[l]);
      }
      if (l != maxLevel)
      {
         for(Block3DPtr block : blocks[l+1])
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
   
   UBLOG(logDEBUG1, "Calculator::initLocalConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void Calculator::initRemoteConnectors()
{
   std::vector< std::vector< Block3DConnectorPtr > > remoteInterConnsCF;
   std::vector< std::vector< Block3DConnectorPtr > > remoteInterConnsFC;
   remoteInterConnsCF.resize(maxLevel+1);
   remoteInterConnsFC.resize(maxLevel+1);

   for(int l = minLevel; l<=maxLevel;l++)
   {
      std::vector<Block3DPtr> blockVector;
      //grid->getBlocks(level, gridRank, true, blockVector);
      grid->getBlocks(l, blockVector);
      for(Block3DPtr block : blockVector)
      {
         int l = block->getLevel();
         block->pushBackRemoteSameLevelConnectors(remoteConns[l]);

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
      }
   }
   //////////////////////////////////////////////////////////////////////////
   //UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - start");
   for (int l = minLevel; l <= maxLevel; l++)
   {
      if (l != maxLevel)
      {
         UBLOG(logDEBUG5, "Calculator::initRemoteConnectors()-initConnectors(remoteInterConns["<<l<<"])");
         for(Block3DConnectorPtr c : remoteInterConns[l] ) c->init();
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
         for(Block3DConnectorPtr c : remoteInterConns[l] ) c->sendTransmitterDataSize();
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
         for(Block3DConnectorPtr c : remoteInterConns[l] ) c->receiveTransmitterDataSize();
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
   for(Block3DConnectorPtr c : connectors ) c->init();
   UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - end");
   //////////////////////////////////////////////////////////////////////////
   //sendTransmitterDataSize
   UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - start");
   for(Block3DConnectorPtr c : connectors ) c->sendTransmitterDataSize();
   UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - end");
   //////////////////////////////////////////////////////////////////////////
   //receiveTransmitterDataSize
   //wenn er hier bei verteilten berechnungen stopped, dann ist vermutlich auf einer seite ein nicht aktiver block!!!
   UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - start");
   for(Block3DConnectorPtr c : connectors ) c->receiveTransmitterDataSize();
   UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - end");

   UBLOG(logDEBUG1, "Calculator::initConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void Calculator::deleteBlocks()
{
   for(std::vector< Block3DPtr > &bs : blocks)
      bs.resize(0);
}
//////////////////////////////////////////////////////////////////////////
void Calculator::deleteConnectors()
{
   deleteConnectors(localConns);
   deleteConnectors(remoteConns);

   deleteConnectors(localInterConns);
   deleteConnectors(remoteInterConns);
}
//////////////////////////////////////////////////////////////////////////
void Calculator::deleteConnectors(std::vector< std::vector< Block3DConnectorPtr > >& conns)
{
   for(std::vector< Block3DConnectorPtr > &c : conns)
      c.resize(0);
}
//////////////////////////////////////////////////////////////////////////


#ifndef ConnectorBlockVisitor_H
#define ConnectorBlockVisitor_H

#include "Block3DVisitor.h"
#include "D3Q27System.h"
#include "Communicator.h"
#include "InterpolationProcessor.h"
#include "CreateTransmittersHelper.h"
#include "ConnectorFactory.h"

class ConnectorBlockVisitor : public Block3DVisitor
{
public:
   ConnectorBlockVisitor(CommunicatorPtr comm, LBMReal nu, InterpolationProcessorPtr iProcessor, ConnectorFactoryPtr cFactory);
   virtual ~ConnectorBlockVisitor();
   void visit(Grid3DPtr grid, Block3DPtr block);
   //////////////////////////////////////////////////////////////////////////
protected:
   void setSameLevelConnectors(Grid3DPtr grid, Block3DPtr block);
   void setRemoteConnectors(Block3DPtr sblock, Block3DPtr tblock, int dir);
   void setInterpolationConnectors(Grid3DPtr grid, Block3DPtr block);
   void setInterpolationConnectors(Block3DPtr fBlockSW, Block3DPtr fBlockSE, Block3DPtr fBlockNW, Block3DPtr fBlockNE, Block3DPtr cBlock, int dir);
   void createTransmitters(Block3DPtr cBlock, Block3DPtr fBlock, int dir,
      CreateTransmittersHelper::IBlock ib,
      CreateTransmittersHelper::TransmitterPtr& senderCF,
      CreateTransmittersHelper::TransmitterPtr& receiverCF,
      CreateTransmittersHelper::TransmitterPtr& senderFC,
      CreateTransmittersHelper::TransmitterPtr& receiverFC);
   CommunicatorPtr comm;
   int gridRank;
   LBMReal nu;
   InterpolationProcessorPtr iProcessor;
   ConnectorFactoryPtr cFactory;
};

#endif //ConnectorBlockVisitor_H


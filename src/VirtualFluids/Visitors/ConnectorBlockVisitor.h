#ifndef ConnectorBlockVisitor_H
#define ConnectorBlockVisitor_H

#include <memory>

#include "Block3DVisitor.h"
#include "D3Q27System.h"
#include "CreateTransmittersHelper.h"

class Grid3D;
class Block3D;
class InterpolationProcessor;
class ConnectorFactory;

class ConnectorBlockVisitor : public Block3DVisitor
{
public:
   ConnectorBlockVisitor(CommunicatorPtr comm, LBMReal nu, std::shared_ptr<InterpolationProcessor> iProcessor, std::shared_ptr<ConnectorFactory> cFactory);
   virtual ~ConnectorBlockVisitor();
      void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;
   //////////////////////////////////////////////////////////////////////////
protected:
   void setSameLevelConnectors(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block);
   void setRemoteConnectors(std::shared_ptr<Block3D> sblock, std::shared_ptr<Block3D> tblock, int dir);
   void setInterpolationConnectors(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block);
   void setInterpolationConnectors(std::shared_ptr<Block3D> fBlockSW, std::shared_ptr<Block3D> fBlockSE, std::shared_ptr<Block3D> fBlockNW, std::shared_ptr<Block3D> fBlockNE, std::shared_ptr<Block3D> cBlock, int dir);
   void createTransmitters(std::shared_ptr<Block3D> cBlock, std::shared_ptr<Block3D> fBlock, int dir,
      CreateTransmittersHelper::IBlock ib,
      CreateTransmittersHelper::TransmitterPtr& senderCF,
      CreateTransmittersHelper::TransmitterPtr& receiverCF,
      CreateTransmittersHelper::TransmitterPtr& senderFC,
      CreateTransmittersHelper::TransmitterPtr& receiverFC);
   CommunicatorPtr comm;
   int gridRank;
   LBMReal nu;
   std::shared_ptr<InterpolationProcessor> iProcessor;
   std::shared_ptr<ConnectorFactory> cFactory;
};

#endif //ConnectorBlockVisitor_H


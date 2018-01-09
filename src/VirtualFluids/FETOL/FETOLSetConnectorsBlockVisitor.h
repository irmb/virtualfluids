#if defined VF_MPI && defined VF_FETOL

#ifndef D3Q27BondSetConnectorsBlockVisitor_H
#define D3Q27BondSetConnectorsBlockVisitor_H

#include "Block3DVisitor.h"
#include "D3Q27System.h"
#include "FETOLCommunicator.h"
#include "InterpolationProcessor.h"
#include "CreateTransmittersHelper.h"

class FETOLSetConnectorsBlockVisitor : public Block3DVisitor
{
public:
   FETOLSetConnectorsBlockVisitor(CommunicatorPtr comm, bool fullConnector, int dirs, LBMReal nue, D3Q27InterpolationProcessorPtr iProcessor);
   virtual ~FETOLSetConnectorsBlockVisitor();
      void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;
   //////////////////////////////////////////////////////////////////////////
protected:
   void setSameLevelConnectors(Grid3DPtr grid, Block3DPtr block);
   void setRemoteConnectors(Block3DPtr sblock, Block3DPtr tblock, int dir);
   void setBundleConnectors(Block3DPtr sblock, Block3DPtr tblock, int dir);
   //void setInterpolationConnectors(Grid3DPtr grid, Block3DPtr block);
   //void setInterpolationConnectors(Grid3DPtr grid, Block3DPtr fBlockSW, Block3DPtr fBlockSE, Block3DPtr fBlockNW, Block3DPtr fBlockNE, Block3DPtr cBlock, int dir);
   //void createTransmitters(Block3DPtr cBlock, Block3DPtr fBlock, int dir, 
   //                        D3Q27CreateTransmittersHelper::TransmitterPtr& senderCF, 
   //                        D3Q27CreateTransmittersHelper::TransmitterPtr& receiverCF, 
   //                        D3Q27CreateTransmittersHelper::TransmitterPtr& senderFC, 
   //                        D3Q27CreateTransmittersHelper::TransmitterPtr& receiverFC);
   BondCommunicatorPtr comm;
   bool fullConnector;
   int dirs;
   int gridBundle;
   int gridRank;
   LBMReal nue;
   D3Q27InterpolationProcessorPtr iProcessor;
};

#endif //D3Q27SETCONNECTORSVISITOR_H

#endif

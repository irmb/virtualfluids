#ifndef D3Q27SETCONNECTORSBLOCKVISITOR_H
#define D3Q27SETCONNECTORSBLOCKVISITOR_H

#include "Block3DVisitor.h"
#include "D3Q27System.h"
#include "Communicator.h"
#include "D3Q27InterpolationProcessor.h"
#include "CreateTransmittersHelper.h"

class D3Q27SetConnectorsBlockVisitor : public Block3DVisitor
{
public:
	D3Q27SetConnectorsBlockVisitor(CommunicatorPtr comm, bool fullConnector, int dirs, LBMReal nue, D3Q27InterpolationProcessorPtr iProcessor);
	virtual ~D3Q27SetConnectorsBlockVisitor();
	void visit(Grid3DPtr grid, Block3DPtr block);
	//////////////////////////////////////////////////////////////////////////
protected:
	void setSameLevelConnectors(Grid3DPtr grid, Block3DPtr block);
	void setRemoteConnectors(Block3DPtr sblock, Block3DPtr tblock, int dir, bool fullConnector);
	void setInterpolationConnectors(Grid3DPtr grid, Block3DPtr block);
	void setInterpolationConnectors(Block3DPtr fBlockSW, Block3DPtr fBlockSE, Block3DPtr fBlockNW, Block3DPtr fBlockNE, Block3DPtr cBlock, int dir);
	void createTransmitters(Block3DPtr cBlock, Block3DPtr fBlock, int dir,
      CreateTransmittersHelper::IBlock ib,
		CreateTransmittersHelper::TransmitterPtr& senderCF, 
		CreateTransmittersHelper::TransmitterPtr& receiverCF, 
		CreateTransmittersHelper::TransmitterPtr& senderFC, 
		CreateTransmittersHelper::TransmitterPtr& receiverFC);
	CommunicatorPtr comm;
	bool fullConnector;
	int dirs;
	int gridRank;
	LBMReal nue;
	D3Q27InterpolationProcessorPtr iProcessor;
};

#endif //D3Q27SETCONNECTORSVISITOR_H

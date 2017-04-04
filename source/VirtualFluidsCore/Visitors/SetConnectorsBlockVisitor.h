#ifndef SETCONNECTORSBLOCKVISITOR_H
#define SETCONNECTORSBLOCKVISITOR_H

#include "Block3DVisitor.h"
#include "D3Q27System.h"
#include "Communicator.h"
#include "InterpolationProcessor.h"
#include "CreateTransmittersHelper.h"

class SetConnectorsBlockVisitor : public Block3DVisitor
{
public:
	SetConnectorsBlockVisitor(CommunicatorPtr comm, bool fullConnector, int dirs, LBMReal nue, InterpolationProcessorPtr iProcessor);
	virtual ~SetConnectorsBlockVisitor();
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
	InterpolationProcessorPtr iProcessor;
};

#endif //D3Q27SETCONNECTORSVISITOR_H

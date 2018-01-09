#ifndef SETCONNECTORSBLOCKVISITOR_H
#define SETCONNECTORSBLOCKVISITOR_H

#include <memory>

#include "Block3DVisitor.h"
#include "D3Q27System.h"

#include "CreateTransmittersHelper.h"

class Grid3D;
class Block3D;
class Communicator;
class InterpolationProcessor;

class SetConnectorsBlockVisitor : public Block3DVisitor
{
public:
	SetConnectorsBlockVisitor(std::shared_ptr<Communicator> comm, bool fullConnector, int dirs, LBMReal nue, std::shared_ptr<InterpolationProcessor> iProcessor);
	virtual ~SetConnectorsBlockVisitor();
	void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;
	//////////////////////////////////////////////////////////////////////////
protected:
	void setSameLevelConnectors(std::shared_ptr<Grid3D> grid, Block3DPtr block);
	void setRemoteConnectors(Block3DPtr sblock, Block3DPtr tblock, int dir, bool fullConnector);
	void setInterpolationConnectors(std::shared_ptr<Grid3D> grid, Block3DPtr block);
	void setInterpolationConnectors(Block3DPtr fBlockSW, Block3DPtr fBlockSE, Block3DPtr fBlockNW, Block3DPtr fBlockNE, Block3DPtr cBlock, int dir);
	void createTransmitters(Block3DPtr cBlock, Block3DPtr fBlock, int dir,
      CreateTransmittersHelper::IBlock ib,
		CreateTransmittersHelper::TransmitterPtr& senderCF, 
		CreateTransmittersHelper::TransmitterPtr& receiverCF, 
		CreateTransmittersHelper::TransmitterPtr& senderFC, 
		CreateTransmittersHelper::TransmitterPtr& receiverFC);
    std::shared_ptr<Communicator> comm;
	bool fullConnector;
	int dirs;
	int gridRank;
	LBMReal nue;
    std::shared_ptr<InterpolationProcessor> iProcessor;
};

#endif //D3Q27SETCONNECTORSVISITOR_H

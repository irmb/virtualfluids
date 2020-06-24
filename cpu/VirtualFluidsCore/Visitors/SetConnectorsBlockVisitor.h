#ifndef SETCONNECTORSBLOCKVISITOR_H
#define SETCONNECTORSBLOCKVISITOR_H

#include <PointerDefinitions.h>

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
	SetConnectorsBlockVisitor(SPtr<Communicator> comm, bool fullConnector, int dirs, LBMReal nue, SPtr<InterpolationProcessor> iProcessor);
	virtual ~SetConnectorsBlockVisitor();
	void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;
	//////////////////////////////////////////////////////////////////////////
protected:
	void setSameLevelConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block);
	void setRemoteConnectors(SPtr<Block3D> sblock, SPtr<Block3D> tblock, int dir, bool fullConnector);
	void setInterpolationConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block);
	void setInterpolationConnectors(SPtr<Block3D> fBlockSW, SPtr<Block3D> fBlockSE, SPtr<Block3D> fBlockNW, SPtr<Block3D> fBlockNE, SPtr<Block3D> cBlock, int dir);
	void createTransmitters(SPtr<Block3D> cBlock, SPtr<Block3D> fBlock, int dir,
      CreateTransmittersHelper::IBlock ib,
		CreateTransmittersHelper::TransmitterPtr& senderCF, 
		CreateTransmittersHelper::TransmitterPtr& receiverCF, 
		CreateTransmittersHelper::TransmitterPtr& senderFC, 
		CreateTransmittersHelper::TransmitterPtr& receiverFC);
    SPtr<Communicator> comm;
	bool fullConnector;
	int dirs;
	int gridRank;
	LBMReal nue;
    SPtr<InterpolationProcessor> iProcessor;
};

#endif //D3Q27SETCONNECTORSVISITOR_H

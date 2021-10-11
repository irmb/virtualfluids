#ifndef GridReaderGenerator_H
#define GridReaderGenerator_H

#include "../GridProvider.h"

#include <vector>
#include <string>
#include <memory>

#include "LBM/LB.h"

class Parameter;
class GridBuilder;
namespace vf
{
	namespace gpu
	{
		class Communicator;
	}
} // namespace vf

class GridGenerator
	: public GridProvider
{
private:
	std::vector<std::string> channelDirections;
	std::vector<std::string> channelBoundaryConditions;

	std::shared_ptr<GridBuilder> builder;

public:
    VIRTUALFLUIDS_GPU_EXPORT GridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager);
	VIRTUALFLUIDS_GPU_EXPORT virtual ~GridGenerator();

	void allocArrays_CoordNeighborGeo() override;
    void allocArrays_BoundaryValues() override;
    void initalValuesDomainDecompostion();

	// communication after coarse to fine
    void initCommunicationArraysForCommAfterFinetoCoarseX(const uint &level, int j, int direction);
    void initCommunicationArraysForCommAfterFinetoCoarseY(const uint &level, int j, int direction);
    void initCommunicationArraysForCommAfterFinetoCoarseZ(const uint &level, int j, int direction);
    void reorderSendIndicesForCommAfterFtoCX(int direction, int level, int j,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoCY(int direction, int level, int j,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoCZ(int direction, int level, int j,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoC(int *sendIndices, int &numberOfSendNeighborsAfterFtoC, int direction,
                                            int level, int j, std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    bool isSparseIndexInICellFCC(uint sizeOfICellFCC, int sparseIndexSend, int level);
    void aggregateNodesInICellCFC(int level, std::vector<uint> &nodesCFC);
    void addUniqueIndexToCommunicationVectors(std::vector<int> &sendIndicesAfterFtoC, int &sparseIndexSend,
                                              std::vector<unsigned int> &sendIndicesForCommAfterFtoCPositions,
                                              uint &posInSendIndices) const;
    void findIfSparseIndexIsInSendIndicesAndAddToCommVectors(int sparseIndex, int *sendIndices, uint numberOfSendIndices,
                                                         std::vector<int> &sendIndicesAfterFtoC,
                                                         std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    void findIndicesNotInCommAfterFtoC(const uint &numberOfSendIndices, int *sendIndices,
                                           std::vector<int> &sendIndicesAfterFtoC, std::vector<int> &sendIndicesOther);
    void reorderRecvIndicesForCommAfterFtoCX(int direction, int level, int j,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCY(int direction, int level, int j,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCZ(int direction, int level, int j,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoC(int *recvIndices, int &numberOfRecvNeighborsAfterFtoC, int direction,
                                            int level, int j, std::vector<uint> &sendIndicesForCommAfterFtoCPositions);

	void allocArrays_BoundaryQs() override;
    void allocArrays_OffsetScale() override;
    void allocArrays_fluidNodeIndices() override;
    void allocArrays_fluidNodeIndicesBorder() override;

	virtual void setDimensions() override;
	virtual void setBoundingBox() override;

	virtual void initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int> > > periodV, std::vector<std::vector<unsigned int> > periodIndex, std::string way) override;
	
private:
	void setPressureValues(int channelSide) const;
	void setPressRhoBC(int sizePerLevel, int level, int channelSide) const;

	void setVelocityValues(int channelSide) const;
	void setVelocity(int level, int sizePerLevel, int channelSide) const;

	void setOutflowValues(int channelSide) const;
	void setOutflow(int level, int sizePerLevel, int channelSide) const;

	void setPressQs(int channelSide) const;
	void setVelocityQs(int channelSide) const;
	void setOutflowQs(int channelSide) const;
	void setNoSlipQs(int channelSide) const;
	void setGeoQs() const;
	void modifyQElement(int channelSide, unsigned int level) const;

	void initalQStruct(QforBoundaryConditions& Q,int channelSide, unsigned int level) const;
	void printQSize(std::string bc,int channelSide, unsigned int level) const;
	void setSizeNoSlip(int channelSide, unsigned int level) const;
	void setSizeGeoQs(unsigned int level) const;
	void setQ27Size(QforBoundaryConditions &Q, real* QQ, unsigned int sizeQ) const;
	bool hasQs(int channelSide, unsigned int level) const;

public:
    void initalGridInformations() override;


private:
    std::string verifyNeighborIndices(int level) const;
    std::string verifyNeighborIndex(int level, int index, int &invalidNodes, int &stopperNodes, int &wrongNeighbors) const;
    std::string checkNeighbor(int level, real x, real y, real z, int index, int& numberOfWrongNeihgbors, int neighborIndex, real neighborX, real neighborY, real neighborZ, std::string direction) const;

};

#endif

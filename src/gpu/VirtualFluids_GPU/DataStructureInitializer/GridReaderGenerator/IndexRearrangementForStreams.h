#ifndef IndexRearrangementForStreams_H
#define IndexRearrangementForStreams_H

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
}

class IndexRearrangementForStreams
{
private:
	std::shared_ptr<GridBuilder> builder;
    std::shared_ptr<Parameter> para;

public:
    IndexRearrangementForStreams(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder);
    
    // communication after coarse to fine
    void initCommunicationArraysForCommAfterFinetoCoarseX(const uint &level, int j, int direction);
    void initCommunicationArraysForCommAfterFinetoCoarseY(const uint &level, int j, int direction);
    void initCommunicationArraysForCommAfterFinetoCoarseZ(const uint &level, int j, int direction);

    // split interpolation cells
    void splitCoarseToFineIntoBorderAndBulk(const uint &level);
    void splitFineToCoarseIntoBorderAndBulk(const uint &level);


private:
    // communication after coarse to fine
    void copyProcessNeighborToCommAfterFtoCX(const uint &level, int j);
    void copyProcessNeighborToCommAfterFtoCY(const uint &level, int j);
    void copyProcessNeighborToCommAfterFtoCZ(const uint &level, int j);

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
    void findIndicesNotInCommAfterFtoC(const uint &numberOfSendOrRecvIndices, int *sendOrReceiveIndices,
                                       std::vector<int> &sendOrReceiveIndicesAfterFtoC,
                                       std::vector<int> &sendOrIndicesOther);

    void reorderRecvIndicesForCommAfterFtoCX(int direction, int level, int j,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCY(int direction, int level, int j,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCZ(int direction, int level, int j,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoC(int *recvIndices, int &numberOfRecvNeighborsAfterFtoC, int direction,
                                            int level, int j, std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    
    // split interpolation cells
    void getGridInterfaceIndicesBorderBulkCF(int level);
    void getGridInterfaceIndicesBorderBulkFC(int level);
};

#endif

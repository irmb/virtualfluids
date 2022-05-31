//! \file IndexRearrangementForStreams.h
//! \ingroup GPU
//! \author Anna Wellmann
//! \ref master thesis of Anna Wellmann

#ifndef IndexRearrangementForStreams_H
#define IndexRearrangementForStreams_H

#include <gpu/VirtualFluids_GPU/DataStructureInitializer/GridProvider.h>

#include <memory>
#include <string>
#include <vector>

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

class IndexRearrangementForStreams
{
private:
    std::shared_ptr<GridBuilder> builder;
    std::shared_ptr<Parameter> para;

public:
    //! \brief construct IndexRearrangementForStreams object
    IndexRearrangementForStreams(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder);

    //////////////////////////////////////////////////////////////////////////
    // communication after coarse to fine
    //////////////////////////////////////////////////////////////////////////
    void initCommunicationArraysForCommAfterFinetoCoarseX(const uint &level, int j, int direction);
    void initCommunicationArraysForCommAfterFinetoCoarseY(const uint &level, int j, int direction);
    void initCommunicationArraysForCommAfterFinetoCoarseZ(const uint &level, int j, int direction);

public:
    //////////////////////////////////////////////////////////////////////////
    // split interpolation cells
    //////////////////////////////////////////////////////////////////////////

    //! \brief split the interpolation cells from coarse to fine into border an bulk
    //! \details For communication hiding, the interpolation cells from the coarse to the fine grid need to be split
    //! into two groups:
    //!
    //! - cells which are at the border between two gpus --> "border"
    //!
    //! - the other cells which are not directly related to the communication betweeen the two gpus --> "bulk"
    //! \ref see master thesis of Anna Wellmann (p. 62-68)
    void splitCoarseToFineIntoBorderAndBulk(const uint &level);

    //! \brief split the interpolation cells from fine to coarse into border an bulk
    //! \details For communication hiding, the interpolation cells from the fine to the coarse grid need to be split
    //! into two groups:
    //!
    //! - cells which are at the border between two gpus --> "border"
    //!
    //! - the other cells which are not directly related to the communication betweeen the two gpus --> "bulk"
    //!
    //! \ref see master thesis of Anna Wellmann (p. 62-68)
    void splitFineToCoarseIntoBorderAndBulk(const uint &level);

private:
    //////////////////////////////////////////////////////////////////////////
    // communication after coarse to fine
    //////////////////////////////////////////////////////////////////////////
    void copyProcessNeighborToCommAfterFtoCX(const uint &level, int indexOfProcessNeighbor);
    void copyProcessNeighborToCommAfterFtoCY(const uint &level, int indexOfProcessNeighbor);
    void copyProcessNeighborToCommAfterFtoCZ(const uint &level, int indexOfProcessNeighbor);

    void reorderSendIndicesForCommAfterFtoCX(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoCY(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoCZ(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoC(int *sendIndices, int &numberOfSendNeighborsAfterFtoC, int direction,
                                            int level, std::vector<uint> &sendIndicesForCommAfterFtoCPositions);

    bool isSparseIndexInICellFCC(uint sizeOfICellFCC, int sparseIndexSend, int level);
    void aggregateNodesInICellCFC(int level, std::vector<uint> &nodesCFC);
    void addUniqueIndexToCommunicationVectors(std::vector<int> &sendIndicesAfterFtoC, int &sparseIndexSend,
                                              std::vector<unsigned int> &sendIndicesForCommAfterFtoCPositions,
                                              uint &posInSendIndices) const;
    void
    findIfSparseIndexIsInSendIndicesAndAddToCommVectors(int sparseIndex, int *sendIndices, uint numberOfSendIndices,
                                                        std::vector<int> &sendIndicesAfterFtoC,
                                                        std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    void findIndicesNotInCommAfterFtoC(const uint &numberOfSendOrRecvIndices, int *sendOrReceiveIndices,
                                       std::vector<int> &sendOrReceiveIndicesAfterFtoC,
                                       std::vector<int> &sendOrIndicesOther);

    void reorderRecvIndicesForCommAfterFtoCX(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCY(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCZ(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoC(int *recvIndices, int &numberOfRecvNeighborsAfterFtoC, int direction,
                                            int level, std::vector<uint> &sendIndicesForCommAfterFtoCPositions);

private:
    //////////////////////////////////////////////////////////////////////////
    // split interpolation cells
    //////////////////////////////////////////////////////////////////////////

    //! \brief This function reorders the arrays of CFC/CFF indices and sets the pointers and sizes of the new
    //! subarrays: \details The coarse cells for interpolation from coarse to fine (iCellCFC) are divided into two
    //! subgroups: border and bulk. The fine cells (iCellCFF) are reordered accordingly. The offset cells (xOffCF,
    //! yOffCF, zOffCF) must be reordered in the same way.
    void getGridInterfaceIndicesBorderBulkCF(int level);

    //! \brief This function reorders the arrays of FCC/FCF indices and return pointers and sizes of the new subarrays:
    //! \details The coarse cells for interpolation from fine to coarse (iCellFCC) are divided into two subgroups:
    //! border and bulk. The fine cells (iCellFCF) are reordered accordingly.
    void getGridInterfaceIndicesBorderBulkFC(int level);
};

#endif

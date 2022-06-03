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
namespace vf::gpu
{
class Communicator;
}

class IndexRearrangementForStreams
{
public:
    //! \brief construct IndexRearrangementForStreams object
    IndexRearrangementForStreams(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder, vf::gpu::Communicator& communicator);

    //////////////////////////////////////////////////////////////////////////
    // communication after coarse to fine
    //////////////////////////////////////////////////////////////////////////

    //! \brief initialize the arrays for the communication after the interpolation from fine to coarse in x direction
    //! \details Only the nodes involved in the interpolation need to be exchanged. Therefore in this method all nodes,
    //! which are part of the interpolation as well as the communication, are identified.
    //!
    //! \ref see master thesis of Anna
    //! Wellmann (p. 59-62: "Reduzieren der auszutauschenden Knoten")
    void initCommunicationArraysForCommAfterFinetoCoarseX(const uint &level, int j, int direction);
    //! \brief initialize the arrays for the communication after the interpolation from fine to coarse in y direction
    //! \details --> see x direction
    void initCommunicationArraysForCommAfterFinetoCoarseY(const uint &level, int j, int direction);
    //! \brief initialize the arrays for the communication after the interpolation from fine to coarse in z direction
    //! \details --> see x direction
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
    //! - the other cells which are not directly related to the communication between the two gpus --> "bulk"
    //!
    //! \ref see master thesis of Anna Wellmann (p. 62-68: "Überdeckung der reduzierten Kommunikation")
    void splitCoarseToFineIntoBorderAndBulk(const uint &level);

    //! \brief split the interpolation cells from fine to coarse into border an bulk
    //! \details For communication hiding, the interpolation cells from the fine to the coarse grid need to be split
    //! into two groups:
    //!
    //! - cells which are at the border between two gpus --> "border"
    //!
    //! - the other cells which are not directly related to the communication between the two gpus --> "bulk"
    //!
    //! \ref see master thesis of Anna Wellmann (p. 62-68: "Überdeckung der reduzierten Kommunikation")
    void splitFineToCoarseIntoBorderAndBulk(const uint &level);

private:
    //////////////////////////////////////////////////////////////////////////
    // communication after coarse to fine
    //////////////////////////////////////////////////////////////////////////

    //! \brief inits pointers for reduced communication after interpolation fine to coarse by copying them from "normal"
    //! communication
    void copyProcessNeighborToCommAfterFtoCX(const uint &level, int indexOfProcessNeighbor);
    void copyProcessNeighborToCommAfterFtoCY(const uint &level, int indexOfProcessNeighbor);
    void copyProcessNeighborToCommAfterFtoCZ(const uint &level, int indexOfProcessNeighbor);

    void reorderSendIndicesForCommAfterFtoCX(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoCY(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoCZ(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);

    //! \brief the send indices are reordered for the communication after the interpolation from fine to coarse
    //! \details The indices of nodes which are part of the interpolation are moved to the front of vector with the send
    //! indices. 
    //! \pre para->getParH(level)->intCF needs to be inititalized 
    //! \param sendIndices is the pointer to the vector with the send indices, which will be reordered in this function
    //! \param numberOfSendNodesAfterFtoC will be set in this method 
    //! \param sendIndicesForCommAfterFtoCPositions stores each sendIndex's positions before reordering
    void reorderSendIndicesForCommAfterFtoC(int *sendIndices, int &numberOfSendNodesAfterFtoC, int direction,
                                            int level, std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    //! \brief check if a sparse index occurs in the ICellFCC
    bool isSparseIndexInICellFCC(uint sizeOfICellFCC, int sparseIndexSend, int level);
    //! \brief aggregate all nodes in the coarse cells for the interpolation in coarse to fine
    //! \details For the coarse cells in the interpolation from coarse to fine only one node is stored. This methods
    //! looks for the other nodes of each cell and puts them into vector. Duplicate nodes are only stored once.
    void aggregateNodesInICellCFC(int level, std::vector<uint> &nodesCFC);
    //! \brief add index to sendIndicesAfterFtoC and sendIndicesForCommAfterFtoCPositions, but omit indices which are already in sendIndicesAfterFtoC
    void addUniqueIndexToCommunicationVectors(std::vector<int> &sendIndicesAfterFtoC, int &sparseIndexSend,
                                              std::vector<unsigned int> &sendIndicesForCommAfterFtoCPositions,
                                              uint &posInSendIndices) const;
    //! \brief find if a sparse index is a send index. If true, call addUniqueIndexToCommunicationVectors()
    void
    findIfSparseIndexIsInSendIndicesAndAddToCommVectors(int sparseIndex, int *sendIndices, uint numberOfSendIndices,
                                                        std::vector<int> &sendIndicesAfterFtoC,
                                                        std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    //! \brief find all indices which are not part of the communication after the interpolation from fine to coarse
    void findIndicesNotInCommAfterFtoC(const uint &numberOfSendOrRecvIndices, int *sendOrReceiveIndices,
                                       std::vector<int> &sendOrReceiveIndicesAfterFtoC,
                                       std::vector<int> &sendOrIndicesOther);

    void reorderRecvIndicesForCommAfterFtoCX(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCY(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCZ(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
                                             
    //! \brief reorder the receive indices in the same way that the send indices were reordered.
    //! \details When the send indices are reordered, the receive indices need to be reordered accordingly.
    //! \pre sendIndicesForCommAfterFtoCPositions should not be empty
    //! \param recvIndices is the pointer to the vector with the receive indices, which will be reordered in this function
    //! \param numberOfRecvNodesAfterFtoC will be set in this function
    //! \param sendIndicesForCommAfterFtoCPositions stores each sendIndex's positions before reordering and is used to reorder the receive indices in the same way
    void reorderRecvIndicesForCommAfterFtoC(int *recvIndices, int &numberOfRecvNodesAfterFtoC, int direction,
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


private:
    std::shared_ptr<GridBuilder> builder;
    std::shared_ptr<Parameter> para;
    vf::gpu::Communicator& communicator;

    // used for tests
    friend class IndexRearrangementForStreamsTest_reorderSendIndices;
};

#endif

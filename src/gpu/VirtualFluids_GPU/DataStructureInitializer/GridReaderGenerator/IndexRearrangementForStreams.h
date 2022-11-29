//! \file IndexRearrangementForStreams.h
//! \ingroup GPU
//! \author Anna Wellmann
//! \details See [master thesis of Anna Wellmann]

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
    //! \brief Construct IndexRearrangementForStreams object
    IndexRearrangementForStreams(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder, vf::gpu::Communicator& communicator);

    //////////////////////////////////////////////////////////////////////////
    // communication after fine to coarse
    //////////////////////////////////////////////////////////////////////////

    //! \brief Initialize the arrays for the communication after the interpolation from fine to coarse in x direction
    //! \details Only the nodes involved in the interpolation need to be exchanged. Therefore in this method all nodes,
    //! which are part of the interpolation as well as the communication, are identified.
    //!See [master thesis of Anna Wellmann (p. 59-62: "Reduzieren der auszutauschenden Knoten")]
    void initCommunicationArraysForCommAfterFinetoCoarseX(uint level, int j, int direction);
    //! \brief Initialize the arrays for the communication after the interpolation from fine to coarse in y direction
    //! \details --> see x direction
    void initCommunicationArraysForCommAfterFinetoCoarseY(uint level, int j, int direction);
    //! \brief Initialize the arrays for the communication after the interpolation from fine to coarse in z direction
    //! \details --> see x direction
    void initCommunicationArraysForCommAfterFinetoCoarseZ(uint level, int j, int direction);

protected:
    //////////////////////////////////////////////////////////////////////////
    // communication after fine to coarse
    //////////////////////////////////////////////////////////////////////////

    //! \brief Initializes the send indices for the communication after the interpolation from fine to coarse
    std::vector<uint> initSendIndicesForCommAfterFToCX(uint level, int indexOfProcessNeighbor, int direction);
    std::vector<uint> initSendIndicesForCommAfterFToCY(uint level, int indexOfProcessNeighbor, int direction);
    std::vector<uint> initSendIndicesForCommAfterFToCZ(uint level, int indexOfProcessNeighbor, int direction);

    //! \brief send sendIndicesForCommAfterFtoCPositions to receiving process and receive
    //! recvIndicesForCommAfterFtoCPositions from neighboring process
    std::vector<uint> exchangeIndicesForCommAfterFtoCX(uint level, int indexOfProcessNeighbor, int direction,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    std::vector<uint> exchangeIndicesForCommAfterFtoCY(uint level, int indexOfProcessNeighbor, int direction,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    std::vector<uint> exchangeIndicesForCommAfterFtoCZ(uint level, int indexOfProcessNeighbor, int direction,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions);

    //! \brief Initializes the send indices for the communication after the interpolation from fine to coarse
    void initRecvIndicesForCommAfterFToCX(uint level, int indexOfProcessNeighbor, int direction, std::vector<uint>& recvIndicesForCommAfterFtoCPositions);
    void initRecvIndicesForCommAfterFToCY(uint level, int indexOfProcessNeighbor, int direction, std::vector<uint>& recvIndicesForCommAfterFtoCPositions);
    void initRecvIndicesForCommAfterFToCZ(uint level, int indexOfProcessNeighbor, int direction, std::vector<uint>& recvIndicesForCommAfterFtoCPositions);

    //! \brief Initializes pointers for reduced communication after the interpolation from fine to coarse by copying
    //! them from "normal" communication
    void copyProcessNeighborToCommAfterFtoCX(uint level, int indexOfProcessNeighbor);
    void copyProcessNeighborToCommAfterFtoCY(uint level, int indexOfProcessNeighbor);
    void copyProcessNeighborToCommAfterFtoCZ(uint level, int indexOfProcessNeighbor);

    //! \brief --> see reorderSendIndicesForCommAfterFtoC
    void reorderSendIndicesForCommAfterFtoCX(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoCY(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderSendIndicesForCommAfterFtoCZ(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);

    //! \brief The send indices are reordered for the communication after the interpolation from fine to coarse
    //! \details The indices of nodes which are part of the interpolation are moved to the front of vector with the send
    //! indices. 
    //! \pre para->getParH(level)->intCF needs to be inititalized 
    //! \param sendIndices is the pointer to the vector with the send indices, which will be reordered in this function
    //! \param numberOfSendNodesAfterFtoC will be set in this method 
    //! \param sendIndicesForCommAfterFtoCPositions stores each sendIndex's positions before reordering
    void reorderSendIndicesForCommAfterFtoC(int *sendIndices, int &numberOfSendNodesAfterFtoC, int direction,
                                            int level, std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    //! \brief Check if a sparse index occurs in the ICellFCC
    bool isSparseIndexInICellFCC(uint sizeOfICellFCC, int sparseIndexSend, int level);
    //! \brief Aggregate all nodes in the coarse cells for the interpolation in coarse to fine
    //! \details For the coarse cells in the interpolation from coarse to fine only one node is stored. This methods
    //! looks for the other nodes of each cell and puts them into vector. Duplicate nodes are only stored once.
    void aggregateNodesInICellCFC(int level, std::vector<uint> &nodesCFC);
    //! \brief Add index to sendIndicesAfterFtoC and sendIndicesForCommAfterFtoCPositions, but omit indices which are already in sendIndicesAfterFtoC
    void addUniqueIndexToCommunicationVectors(std::vector<int> &sendIndicesAfterFtoC, int &sparseIndexSend,
                                              std::vector<unsigned int> &sendIndicesForCommAfterFtoCPositions,
                                              uint &posInSendIndices) const;
    //! \brief Find if a sparse index is a send index. If true, call addUniqueIndexToCommunicationVectors()
    void
    findIfSparseIndexIsInSendIndicesAndAddToCommVectors(int sparseIndex, int *sendIndices, uint numberOfSendIndices,
                                                        std::vector<int> &sendIndicesAfterFtoC,
                                                        std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    //! \brief Find all indices which are not part of the communication after the interpolation from fine to coarse
    void findIndicesNotInCommAfterFtoC(const uint &numberOfSendOrRecvIndices, int *sendOrReceiveIndices,
                                       std::vector<int> &sendOrReceiveIndicesAfterFtoC,
                                       std::vector<int> &sendOrIndicesOther);

    //! \brief --> see reorderRecvIndicesForCommAfterFtoC
    void reorderRecvIndicesForCommAfterFtoCX(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCY(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
    void reorderRecvIndicesForCommAfterFtoCZ(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions);
                                             
    //! \brief Reorder the receive indices in the same way that the send indices were reordered.
    //! \details When the send indices are reordered, the receive indices need to be reordered accordingly.
    //! \pre sendIndicesForCommAfterFtoCPositions should not be empty
    //! \param recvIndices is the pointer to the vector with the receive indices, which will be reordered in this function
    //! \param numberOfRecvNodesAfterFtoC will be set in this function
    //! \param sendIndicesForCommAfterFtoCPositions stores each sendIndex's positions before reordering and is used to reorder the receive indices in the same way
    void reorderRecvIndicesForCommAfterFtoC(int *recvIndices, int &numberOfRecvNodesAfterFtoC, int direction,
                                            int level, std::vector<uint> &sendIndicesForCommAfterFtoCPositions);

private:
    std::shared_ptr<GridBuilder> builder;
    std::shared_ptr<Parameter> para;
    vf::gpu::Communicator& communicator;

    // used for tests
    friend class IndexRearrangementForStreamsTest_reorderSendIndices;
    friend class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX;
    friend class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY;
    friend class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ;
};

#endif

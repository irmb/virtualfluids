//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Anna Wellmann
//! \details See [master thesis of Anna Wellmann]
//=======================================================================================
#ifndef IndexRearrangementForStreams_H
#define IndexRearrangementForStreams_H

#include <memory>
#include <vector>

#include <basics/DataTypes.h>

class Parameter;
class GridBuilder;
namespace vf::parallel
{
class Communicator;
}

class IndexRearrangementForStreams
{
public:
    //! \brief Construct IndexRearrangementForStreams object
    IndexRearrangementForStreams(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder, vf::parallel::Communicator& communicator);

    virtual ~IndexRearrangementForStreams() = default;

    //////////////////////////////////////////////////////////////////////////
    // communication after fine to coarse
    //////////////////////////////////////////////////////////////////////////

    //! \brief Initialize the arrays for the communication after the interpolation from fine to coarse in x direction
    //! \details Only the nodes involved in the interpolation need to be exchanged. Therefore in this method all nodes,
    //! which are part of the interpolation as well as the communication, are identified.
    //!See [master thesis of Anna Wellmann (p. 59-62: "Reduzieren der auszutauschenden Knoten")]
    virtual void initCommunicationArraysForCommAfterFinetoCoarseX(uint level, int j, int direction) const;
    //! \brief Initialize the arrays for the communication after the interpolation from fine to coarse in y direction
    //! \details --> see x direction
    virtual void initCommunicationArraysForCommAfterFinetoCoarseY(uint level, int j, int direction) const;
    //! \brief Initialize the arrays for the communication after the interpolation from fine to coarse in z direction
    //! \details --> see x direction
    virtual void initCommunicationArraysForCommAfterFinetoCoarseZ(uint level, int j, int direction) const;

protected:
    //////////////////////////////////////////////////////////////////////////
    // communication after fine to coarse
    //////////////////////////////////////////////////////////////////////////

    //! \brief Initializes the send indices for the communication after the interpolation from fine to coarse
    std::vector<uint> initSendIndicesForCommAfterFToCX(uint level, int indexOfProcessNeighbor, int direction) const;
    std::vector<uint> initSendIndicesForCommAfterFToCY(uint level, int indexOfProcessNeighbor, int direction) const;
    std::vector<uint> initSendIndicesForCommAfterFToCZ(uint level, int indexOfProcessNeighbor, int direction) const;

    //! \brief send sendIndicesForCommAfterFtoCPositions to receiving process and receive
    //! recvIndicesForCommAfterFtoCPositions from neighboring process
    std::vector<uint> exchangeIndicesForCommAfterFtoCX(uint level, int indexOfProcessNeighbor,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    std::vector<uint> exchangeIndicesForCommAfterFtoCY(uint level, int indexOfProcessNeighbor,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    std::vector<uint> exchangeIndicesForCommAfterFtoCZ(uint level, int indexOfProcessNeighbor,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;

    //! \brief Initializes the send indices for the communication after the interpolation from fine to coarse
    void initRecvIndicesForCommAfterFToCX(uint level, int indexOfProcessNeighbor, int direction,
                                          std::vector<uint> &recvIndicesForCommAfterFtoCPositions) const;
    void initRecvIndicesForCommAfterFToCY(uint level, int indexOfProcessNeighbor, int direction,
                                          std::vector<uint> &recvIndicesForCommAfterFtoCPositions) const;
    void initRecvIndicesForCommAfterFToCZ(uint level, int indexOfProcessNeighbor, int direction,
                                          std::vector<uint> &recvIndicesForCommAfterFtoCPositions) const;

    //! \brief Initializes pointers for reduced communication after the interpolation from fine to coarse by copying
    //! them from "normal" communication
    void copyProcessNeighborToCommAfterFtoCX(uint level, int indexOfProcessNeighbor) const;
    void copyProcessNeighborToCommAfterFtoCY(uint level, int indexOfProcessNeighbor) const;
    void copyProcessNeighborToCommAfterFtoCZ(uint level, int indexOfProcessNeighbor) const;

    //! \brief --> see reorderSendIndicesForCommAfterFtoC
    void reorderSendIndicesForCommAfterFtoCX(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    void reorderSendIndicesForCommAfterFtoCY(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    void reorderSendIndicesForCommAfterFtoCZ(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;

    //! \brief The send indices are reordered for the communication after the interpolation from fine to coarse
    //! \details The indices of nodes which are part of the interpolation are moved to the front of vector with the send
    //! indices.
    //! \pre para->getParH(level)->intCF needs to be inititalized
    //! \param sendIndices is the pointer to the vector with the send indices, which will be reordered in this function
    //! \param numberOfSendNodesAfterFtoC will be set in this method
    //! \param sendIndicesForCommAfterFtoCPositions stores each sendIndex's positions before reordering
    void reorderSendIndicesForCommAfterFtoC(int *sendIndices, int &numberOfSendNodesAfterFtoC, int direction,
                                            int level, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    //! \brief Check if a sparse index occurs in the coarse nodes for the interpolation from fine to coarse
    bool isSparseIndexInCoarseIndexForFtoC(uint numberOfCoarseNodesForFtoC, int sparseIndexSend, int level) const;
    //! \brief Aggregate all nodes in the coarse cells for the interpolation in coarse to fine
    //! \details For the coarse cells in the interpolation from coarse to fine only one node is stored. This methods
    //! looks for the other nodes of each cell and puts them into vector. Duplicate nodes are only stored once.
    void aggregateCoarseNodesForCtoF(int level, std::vector<uint> &nodesCFC) const;
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
    static void findIndicesNotInCommAfterFtoC(const uint &numberOfSendOrRecvIndices, int *sendOrReceiveIndices,
                                              std::vector<int> &sendOrReceiveIndicesAfterFtoC,
                                              std::vector<int> &sendOrIndicesOther);

    //! \brief --> see reorderRecvIndicesForCommAfterFtoC
    void reorderRecvIndicesForCommAfterFtoCX(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    void reorderRecvIndicesForCommAfterFtoCY(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    void reorderRecvIndicesForCommAfterFtoCZ(int direction, int level, int indexOfProcessNeighbor,
                                             std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;

    //! \brief Reorder the receive indices in the same way that the send indices were reordered.
    //! \details When the send indices are reordered, the receive indices need to be reordered accordingly.
    //! \pre sendIndicesForCommAfterFtoCPositions should not be empty
    //! \param recvIndices is the pointer to the vector with the receive indices, which will be reordered in this
    //! function \param numberOfRecvNodesAfterFtoC will be set in this function \param
    //! sendIndicesForCommAfterFtoCPositions stores each sendIndex's positions before reordering and is used to reorder
    //! the receive indices in the same way
    void reorderRecvIndicesForCommAfterFtoC(int *recvIndices, int &numberOfRecvNodesAfterFtoC, int direction, int level,
                                            std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;

private:
    std::shared_ptr<GridBuilder> builder;
    std::shared_ptr<Parameter> para;
    vf::parallel::Communicator &communicator;

    // used for tests
    friend class IndexRearrangementForStreamsTest_reorderSendIndices;
    friend class IndexRearrangementForStreamsTest_reorderRecvIndicesX;
    friend class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX;
    friend class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY;
    friend class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ;
};

#endif

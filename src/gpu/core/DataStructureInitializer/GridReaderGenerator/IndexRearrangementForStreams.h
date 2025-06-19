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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_DataStructureInitializer DataStructureInitializer
//! \ingroup gpu_core core
//! \{
//! \author Anna Wellmann
//=======================================================================================
#ifndef IndexRearrangementForStreams_H
#define IndexRearrangementForStreams_H

#include "Calculation/Calculation.h"
#include <memory>
#include <vector>

#include <basics/DataTypes.h>

class Parameter;
class GridBuilder;
namespace vf::parallel
{
class Communicator;
}

//! \brief class that is used to rearrange the arrays of node indices for communication between gpus. The rearrangement is
//! needed for communication hiding with cuda streams
//! \details This class changes the order of the node indices that are needed for communication between gpus. The indices are
//! reordered so that they can be split into two groups: nodes that are part if the interpolation between grid levels, and
//! nodes that are not. These groups are needed for communication hiding. For details see [master thesis of Anna Wellmann]
class IndexRearrangementForStreams
{
public:
    IndexRearrangementForStreams(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder, vf::parallel::Communicator& communicator);

    virtual ~IndexRearrangementForStreams() = default;

    //////////////////////////////////////////////////////////////////////////
    // communication after fine to coarse
    //////////////////////////////////////////////////////////////////////////

    //! \brief Initialize the arrays for the communication after the interpolation from fine to coarse
    //! \details Only the nodes involved in the interpolation need to be exchanged. Therefore in this method all nodes,
    //! which are part of the interpolation as well as the communication, are identified.
    //!See [master thesis of Anna Wellmann (p. 59-62: "Reduzieren der auszutauschenden Knoten")]
    virtual void initCommunicationArraysForCommAfterFinetoCoarse(ProcessNeighbor27& sendNeighborHost, ProcessNeighbor27& sendNeighborDevice, ProcessNeighbor27& sendNeighborAfterFtoCHost,  ProcessNeighbor27& sendNeighborAfterFtoCDevice, ProcessNeighbor27& recvNeighborHost,ProcessNeighbor27& recvNeighborDevice, ProcessNeighbor27& recvNeighborAfterFtoCHost, ProcessNeighbor27& recvNeighborAfterFtoCDevice, int level, int direction) const;


protected:
    //////////////////////////////////////////////////////////////////////////
    // communication after fine to coarse
    //////////////////////////////////////////////////////////////////////////

    //! \brief send sendIndicesForCommAfterFtoCPositions to receiving process and receive
    //! recvIndicesForCommAfterFtoCPositions from neighboring process
    std::vector<uint> exchangeIndicesForCommAfterFtoC(ProcessNeighbor27& sendNeighbor, ProcessNeighbor27& recvNeighbor,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;

    //! \brief Initializes pointers for reduced communication after the interpolation from fine to coarse by copying
    //! them from "normal" communication
    void copyProcessNeighborToCommAfterFtoC(ProcessNeighbor27& neighbor, ProcessNeighbor27& neighborAfterFtoC) const;
    void setNumberOfNodes( ProcessNeighbor27& neighborAfterFtoCHost,  ProcessNeighbor27& neighborAfterFtoCDevice, uint numberOfNodes) const;


    //! \brief The send indices are reordered for the communication after the interpolation from fine to coarse
    //! \details The indices of nodes which are part of the interpolation are moved to the front of vector with the send
    //! indices.
    //! \pre para->getParH(level)->intCF needs to be initialized
    //! \param sendIndicesForCommAfterFtoCPositions stores each sendIndex's positions before reordering
    uint reorderSendIndicesForCommAfterFtoC(ProcessNeighbor27& sendNeighbor, int direction,
                                            int level, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    //! \brief Check if a sparse index occurs in the coarse nodes for the interpolation from fine to coarse
    bool isSparseIndexInCoarseIndexForFtoC(uint numberOfCoarseNodesForFtoC, uint sparseIndexSend, int level) const;
    //! \brief Aggregate all nodes in the coarse cells for the interpolation in coarse to fine
    //! \details For the coarse cells in the interpolation from coarse to fine only one node is stored. This methods
    //! looks for the other nodes of each cell and puts them into vector. Duplicate nodes are only stored once.
    void aggregateCoarseNodesForCtoF(int level, std::vector<uint> &nodesCFC) const;
    //! \brief Add index to sendIndicesAfterFtoC and sendIndicesForCommAfterFtoCPositions, but omit indices which are already in sendIndicesAfterFtoC
    void addUniqueIndexToCommunicationVectors(std::vector<uint> &sendIndicesAfterFtoC, uint &sparseIndexSend,
                                              std::vector<uint> &sendIndicesForCommAfterFtoCPositions,
                                              uint &posInSendIndices) const;
    //! \brief Find if a sparse index is a send index. If true, call addUniqueIndexToCommunicationVectors()
    void
    findIfSparseIndexIsInSendIndicesAndAddToCommVectors(uint sparseIndex, const uint *sendIndices, uint numberOfSendIndices,
                                                        std::vector<uint> &sendIndicesAfterFtoC,
                                                        std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;
    //! \brief Find all indices which are not part of the communication after the interpolation from fine to coarse
    static void findIndicesNotInCommAfterFtoC(const uint &numberOfSendOrRecvIndices, const uint *sendOrReceiveIndices,
                                              std::vector<uint> &sendOrReceiveIndicesAfterFtoC,
                                              std::vector<uint> &sendOrIndicesOther);

    //! \brief Reorder the receive indices in the same way that the send indices were reordered.
    //! \details When the send indices are reordered, the receive indices need to be reordered accordingly.
    //! \pre sendIndicesForCommAfterFtoCPositions should not be empty
    //! \param recvIndices is the pointer to the vector with the receive indices, which will be reordered in this
    //! function \param numberOfRecvNodesAfterFtoC will be set in this function \param
    //! sendIndicesForCommAfterFtoCPositions stores each sendIndex's positions before reordering and is used to reorder
    //! the receive indices in the same way
    uint reorderRecvIndicesForCommAfterFtoC(ProcessNeighbor27& recvNeighbor, int direction, int level,
                                            const std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;

private:
    std::shared_ptr<GridBuilder> builder;
    std::shared_ptr<Parameter> para;
    vf::parallel::Communicator &communicator;

    // used for tests
    friend class IndexRearrangementForStreamsTest_reorderSendIndices;
    friend class IndexRearrangementForStreamsTest_reorderRecvIndices;
    friend class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoC;
};

#endif

//! \}

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
#include <array>

#include <basics/DataTypes.h>

namespace vf::parallel
{
class Communicator;
}

namespace vf::gpu {

class Parameter;
class GridBuilder;


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
    //! See [master thesis of Anna Wellmann (p. 59-62: "Reduzieren der auszutauschenden Knoten")]
    //! \return array of Process Neighbors after FtoC in same order as input parameters
    virtual std::array<ProcessNeighbor27, 4> initCommunicationArraysForCommAfterFinetoCoarse(
        const ProcessNeighbor27& sendNeighborHost, const ProcessNeighbor27& sendNeighborDevice,
        const ProcessNeighbor27& recvNeighborHost, const ProcessNeighbor27& recvNeighborDevice, int level,
        int direction) const;

protected:
    //////////////////////////////////////////////////////////////////////////
    // communication after fine to coarse
    //////////////////////////////////////////////////////////////////////////

    //! \brief send sendIndicesForCommAfterFtoCPositions to receiving process and receive
    //! recvIndicesForCommAfterFtoCPositions from neighboring process
    //! \return recvIndicesForCommAfterFtoCPositions
    std::vector<uint> exchangeIndicesForCommAfterFtoC(const ProcessNeighbor27& sendNeighbor, const ProcessNeighbor27& recvNeighbor,
                                                      const std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const;

    //! \brief Initializes pointers for reduced communication after the interpolation from fine to coarse by copying
    //! them from "normal" communication
    static ProcessNeighbor27 makeProcessNeighborToCommAfterFtoC(const ProcessNeighbor27& neighbor, uint numberOfNodes);


    //! \brief The send indices are reordered for the communication after the interpolation from fine to coarse
    //! \details The indices of nodes which are part of the interpolation are moved to the front of vector with the send
    //! indices.
    //! \pre para->getParH(level)->intCF needs to be initialized
    //! \return each sendindex's positions before reordering
    std::vector<uint> reorderSendIndicesForCommAfterFtoC(uint* indices, int direction, int level) const;
    //! \brief Aggregate all nodes in the coarse cells for the interpolation in coarse to fine
    //! \details For the coarse cells in the interpolation from coarse to fine only one node is stored. This methods
    //! looks for the other nodes of each cell and puts them into vector. Duplicate nodes are only stored once.
    std::vector<uint> aggregateCoarseNodesForCtoF(int level) const;

    //! \brief Find all indices which are not part of the communication after the interpolation from fine to coarse
    static std::vector<uint> findIndicesNotInCommAfterFtoC(uint numberOfIndices, const uint *indices, const std::vector<uint> &indicesAfterFtoC);

    //! \brief Reorder the receive indices in the same way that the send indices were reordered.
    //! \details When the send indices are reordered, the receive indices need to be reordered accordingly.
    //! \pre sendIndicesForCommAfterFtoCPositions should not be empty
    //! \param recvIndices is the pointer to the vector with the receive indices, which will be reordered in this
    //! function \param recvIndicesForCommAfterFtoCPositions stores the positions of the neighbors send indices before reordering and is used to reorder
    //! the receive indices in the same way
    void reorderRecvIndicesForCommAfterFtoC(uint* recvIndices, int direction, int level,
                                            const std::vector<uint> &recvIndicesForCommAfterFtoCPositions) const;

private:
    std::shared_ptr<GridBuilder> builder;
    std::shared_ptr<Parameter> para;
    vf::parallel::Communicator &communicator;

    // used for tests
    friend class IndexRearrangementForStreamsTest_reorderSendIndices;
    friend class IndexRearrangementForStreamsTest_reorderRecvIndices;
    friend class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoC;
};

}

#endif

//! \}

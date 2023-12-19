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
//! \addtogroup cpu_Visitors Visitors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================
#include "RenumberGridVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
//#include <mpi.h>

RenumberGridVisitor::RenumberGridVisitor(std::shared_ptr<vf::parallel::Communicator> com) : comm(com) {}

//////////////////////////////////////////////////////////////////////////
void RenumberGridVisitor::visit(SPtr<Grid3D> grid)
{
    int counter = 0;

    // UBLOG(logDEBUG5, "RenumberGridVisitor::visit() - start");
    std::vector<SPtr<Block3D>> blocks;
    //   int gridRank = grid->getRank();
    int size;
    // MPI_Comm_size(MPI_COMM_WORLD, &size);
    size = comm->getNumberOfProcesses();

    int minInitLevel = grid->getCoarsestInitializedLevel();
    int maxInitLevel = grid->getFinestInitializedLevel();

    Grid3D::BlockIDMap &blockIdMap = grid->getBlockIDs();
    blockIdMap.clear();

    for (int rank = 0; rank < size; rank++) {
        for (int level = minInitLevel; level <= maxInitLevel; level++) {
            std::vector<SPtr<Block3D>> blockVector;
            grid->getBlocks(level, blockVector);
            for (SPtr<Block3D> block : blockVector) {
                if (block->getRank() == rank) {
                    block->setGlobalID(counter);
                    blockIdMap.insert(std::make_pair(counter, block));
                    counter++;
                }
            }
        }
    }

    // UBLOG(logDEBUG5, "RenumberGridVisitor::visit() - end");
}

//! \}

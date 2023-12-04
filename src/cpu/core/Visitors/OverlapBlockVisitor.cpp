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
//! \file OverlapBlockVisitor.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher, Soeren Freudiger, Sebastian Geller
//=======================================================================================
#include "OverlapBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"

OverlapBlockVisitor::OverlapBlockVisitor(int levelDepth /*shut be maxGridLevel*/, bool includeNotActiveBlocks)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), levelDepth(levelDepth), includeNotActiveBlocks(includeNotActiveBlocks)
{
}
//////////////////////////////////////////////////////////////////////////
void OverlapBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    int ix1, ix2, ix3, level;
    ix1   = block->getX1();
    ix2   = block->getX2();
    ix3   = block->getX3();
    level = block->getLevel();

    int nix1, nix2, nix3, nlev;
    int neighix1, neighix2, neighix3, neighlev;
    std::vector<SPtr<Block3D>> neighbors;
    grid->getAllNeighbors(ix1, ix2, ix3, level, this->levelDepth, neighbors);
    //   bool hasAdded = false;
    for (size_t i = 0; i < neighbors.size(); i++) {
        if ((neighbors[i]->isActive() || includeNotActiveBlocks) && neighbors[i]->getLevel() > level) {
            neighix1 = neighbors[i]->getX1();
            neighix2 = neighbors[i]->getX2();
            neighix3 = neighbors[i]->getX3();
            neighlev = neighbors[i]->getLevel();
            nix1     = neighix1 >> 1;
            nix2     = neighix2 >> 1;
            nix3     = neighix3 >> 1;
            nlev     = neighlev - 1;

            if (nlev != level) {
                throw UbException(UB_EXARGS, "OverlapBlockVisitor::adaptBlock - leveldifferenz passt nicht, block: " +
                                                 block->toString());
            }

            SPtr<Block3D> newBlock = grid->getBlock(nix1, nix2, nix3, nlev);
            if (!newBlock) {
                newBlock = SPtr<Block3D>(new Block3D(nix1, nix2, nix3, nlev));
                grid->addBlock(newBlock);
                // hasAdded=true;
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
std::string OverlapBlockVisitor::getSpecificDescription()
{
    std::string str("Overlap:");
    return str;
}

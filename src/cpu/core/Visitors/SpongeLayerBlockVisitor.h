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
//! \file SpongeLayerBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef SpongeLayerBlockVisitor_h__
#define SpongeLayerBlockVisitor_h__

#include "Block3DVisitor.h"
#include "D3Q27System.h"

class Grid3D;
class Block3D;
class GbCuboid3D;
class LBMKernel;

//! \brief Set sponge layer for all blocks inside boundingBox
//! \details This visitor sets viscosity gradient inside bounding box.
//! \author K. Kutscher
class SpongeLayerBlockVisitor : public Block3DVisitor
{
public:
    SpongeLayerBlockVisitor(SPtr<GbCuboid3D> boundingBox, SPtr<LBMKernel> kernel, real nue, int dir);
    ~SpongeLayerBlockVisitor() override;

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    SPtr<GbCuboid3D> boundingBox;
    SPtr<LBMKernel> kernel;
    real nue;
    int dir;
};

#endif // SetSpongeLayerBlockVisitor_h__

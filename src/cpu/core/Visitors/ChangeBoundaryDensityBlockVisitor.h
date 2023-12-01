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
//! \file ChangeBoundaryDensityBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef ChangeBoundaryDensityBlockVisitor_h__
#define ChangeBoundaryDensityBlockVisitor_h__

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "lbm/constants/D3Q27.h"

class Block3D;
class Grid3D;
class BoundaryConditions;

class ChangeBoundaryDensityBlockVisitor : public Block3DVisitor
{
public:
    ChangeBoundaryDensityBlockVisitor(real oldBoundaryDensity, real newBoundaryDensity);
    ~ChangeBoundaryDensityBlockVisitor() override;

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    real oldBoundaryDensity;
    real newBoundaryDensity;
    SPtr<BoundaryConditions> bcPtr;
};
#endif // ChangeBoundaryDensityBlockVisitor_h__

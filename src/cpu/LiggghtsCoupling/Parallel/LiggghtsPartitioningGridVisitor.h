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
//! \file LiggghtsPartitioningGridVisitor.h
//! \ingroup LiggghtsCoupling
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef LiggghtsPartitioningGridVisitor_h
#define LiggghtsPartitioningGridVisitor_h

#include <lammps.h>
#include <vector>
#include "basics/PointerDefinitions.h"
#include "cpu/core/Visitors/Grid3DVisitor.h"

class LiggghtsCouplingWrapper;
class Grid3D;

class LiggghtsPartitioningGridVisitor : public Grid3DVisitor
{
public:
    LiggghtsPartitioningGridVisitor(int nx, int ny, int nz, LAMMPS_NS::LAMMPS *lmp);

    ~LiggghtsPartitioningGridVisitor() override;

    void visit(SPtr<Grid3D> grid) override;

private:
    int nx, ny, nz;
    LAMMPS_NS::LAMMPS &lmp;
    int npx{ 0 }, npy{ 0 }, npz{ 0 };
    std::vector<int> xVal, yVal, zVal;
};
#endif
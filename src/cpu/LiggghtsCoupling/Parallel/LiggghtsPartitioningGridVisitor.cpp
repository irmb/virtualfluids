#include "LiggghtsPartitioningGridVisitor.h"
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
//! \file LiggghtsPartitioningGridVisitor.cpp
//! \ingroup LiggghtsCoupling
//! \author Konstantin Kutscher
//=======================================================================================

#include <comm.h>
#include "LiggghtsPartitioningGridVisitor.h"
#include "cpu/core/Simulation/Grid3D.h"
#include "cpu/core/Simulation/Block3D.h"

LiggghtsPartitioningGridVisitor::LiggghtsPartitioningGridVisitor(int nx, int ny, int nz, LAMMPS_NS::LAMMPS *lmp) : nx(nx), ny(ny), nz(nz), lmp(*lmp)
{
 
}

LiggghtsPartitioningGridVisitor::~LiggghtsPartitioningGridVisitor()
{

}

void LiggghtsPartitioningGridVisitor::visit(SPtr<Grid3D> grid)
{
    npx = lmp.comm->procgrid[0];
    npy = lmp.comm->procgrid[1];
    npz = lmp.comm->procgrid[2];

    for (int i = 0; i <= npx; i++)
        xVal.push_back(round(lmp.comm->xsplit[i] * (double)nx));
    for (int i = 0; i <= npy; i++)
        yVal.push_back(round(lmp.comm->ysplit[i] * (double)ny));
    for (int i = 0; i <= npz; i++)
        zVal.push_back(round(lmp.comm->zsplit[i] * (double)nz));

    UbTupleInt3 blockNX = grid->getBlockNX();

    for (int iX = 0; iX < xVal.size() - 1; ++iX) {
        for (int iY = 0; iY < yVal.size() - 1; ++iY) {
            for (int iZ = 0; iZ < zVal.size() - 1; ++iZ) {

                int rank = (int)lmp.comm->grid2proc[iX][iY][iZ];
                int blockX1 = xVal[iX] / val<1>(blockNX);
                int blockX2 = yVal[iY] / val<2>(blockNX);
                int blockX3 = zVal[iZ] / val<3>(blockNX);
                SPtr<Block3D> block = grid->getBlock(blockX1, blockX2, blockX3, 0);
                block->setRank(rank);
            }
        }
    }

    xVal.clear();
    yVal.clear();
    zVal.clear();
}

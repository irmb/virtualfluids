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
//! \file Distribution.cpp
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "Distribution.h"

#include <stdio.h>

#include "grid/distributions/D3Q27.h"
#include "lbm/constants/D3Q27.h"
#include "grid/Grid.h"
using namespace vf::lbm::dir;

Distribution DistributionHelper::getDistribution27() 
{
    Distribution d27;
    d27.name = "D3Q27";
    d27.dir_start = STARTDIR;
    d27.dir_end = ENDDIR;

    d27.dirs = new int[(ENDDIR + 1) * DIMENSION];

    d27.directions = new Direction[ENDDIR + 1];
    d27.directions[E] = Direction(DIR_27_E_X, DIR_27_E_Y, DIR_27_E_Z);
    d27.directions[W] = Direction(DIR_27_W_X, DIR_27_W_Y, DIR_27_W_Z);
    d27.directions[N] = Direction(DIR_27_N_X, DIR_27_N_Y, DIR_27_N_Z);
    d27.directions[S] = Direction(DIR_27_S_X, DIR_27_S_Y, DIR_27_S_Z);
    d27.directions[T] = Direction(DIR_27_T_X, DIR_27_T_Y, DIR_27_T_Z);
    d27.directions[B] = Direction(DIR_27_B_X, DIR_27_B_Y, DIR_27_B_Z);

    d27.directions[NE] = Direction(DIR_27_NE_X, DIR_27_NE_Y, DIR_27_NE_Z);
    d27.directions[SW] = Direction(DIR_27_SW_X, DIR_27_SW_Y, DIR_27_SW_Z);
    d27.directions[SE] = Direction(DIR_27_SE_X, DIR_27_SE_Y, DIR_27_SE_Z);
    d27.directions[NW] = Direction(DIR_27_NW_X, DIR_27_NW_Y, DIR_27_NW_Z);

    d27.directions[TE] = Direction(DIR_27_TE_X, DIR_27_TE_Y, DIR_27_TE_Z);
    d27.directions[BW] = Direction(DIR_27_BW_X, DIR_27_BW_Y, DIR_27_BW_Z);
    d27.directions[BE] = Direction(DIR_27_BE_X, DIR_27_BE_Y, DIR_27_BE_Z);
    d27.directions[TW] = Direction(DIR_27_TW_X, DIR_27_TW_Y, DIR_27_TW_Z);

    d27.directions[TN] = Direction(DIR_27_TN_X, DIR_27_TN_Y, DIR_27_TN_Z);
    d27.directions[BS] = Direction(DIR_27_BS_X, DIR_27_BS_Y, DIR_27_BS_Z);
    d27.directions[BN] = Direction(DIR_27_BN_X, DIR_27_BN_Y, DIR_27_BN_Z);
    d27.directions[TS] = Direction(DIR_27_TS_X, DIR_27_TS_Y, DIR_27_TS_Z);

    d27.directions[REST] = Direction(DIR_27_REST_X, DIR_27_REST_Y, DIR_27_REST_Z);

    d27.directions[TNE] = Direction(DIR_27_TNE_X, DIR_27_TNE_Y, DIR_27_TNE_Z);
    d27.directions[TNW] = Direction(DIR_27_TNW_X, DIR_27_TNW_Y, DIR_27_TNW_Z);
    d27.directions[TSE] = Direction(DIR_27_TSE_X, DIR_27_TSE_Y, DIR_27_TSE_Z);
    d27.directions[TSW] = Direction(DIR_27_TSW_X, DIR_27_TSW_Y, DIR_27_TSW_Z);

    d27.directions[BNE] = Direction(DIR_27_BNE_X, DIR_27_BNE_Y, DIR_27_BNE_Z);
    d27.directions[BNW] = Direction(DIR_27_BNW_X, DIR_27_BNW_Y, DIR_27_BNW_Z);
    d27.directions[BSE] = Direction(DIR_27_BSE_X, DIR_27_BSE_Y, DIR_27_BSE_Z);
    d27.directions[BSW] = Direction(DIR_27_BSW_X, DIR_27_BSW_Y, DIR_27_BSW_Z);


    d27.dirs[E * 3    ] = DIR_27_E_X;
    d27.dirs[E * 3 + 1] = DIR_27_E_Y;
    d27.dirs[E * 3 + 2] = DIR_27_E_Z;

    d27.dirs[W * 3    ] = DIR_27_W_X;
    d27.dirs[W * 3 + 1] = DIR_27_W_Y;
    d27.dirs[W * 3 + 2] = DIR_27_W_Z;
    
    d27.dirs[N * 3    ] = DIR_27_N_X;
    d27.dirs[N * 3 + 1] = DIR_27_N_Y;
    d27.dirs[N * 3 + 2] = DIR_27_N_Z;

    d27.dirs[S * 3    ] = DIR_27_S_X;
    d27.dirs[S * 3 + 1] = DIR_27_S_Y;
    d27.dirs[S * 3 + 2] = DIR_27_S_Z;
    
    d27.dirs[T * 3    ] = DIR_27_T_X;
    d27.dirs[T * 3 + 1] = DIR_27_T_Y;
    d27.dirs[T * 3 + 2] = DIR_27_T_Z;
    
    d27.dirs[B * 3    ] = DIR_27_B_X;
    d27.dirs[B * 3 + 1] = DIR_27_B_Y;
    d27.dirs[B * 3 + 2] = DIR_27_B_Z;

    d27.dirs[NE * 3    ] = DIR_27_NE_X;
    d27.dirs[NE * 3 + 1] = DIR_27_NE_Y;
    d27.dirs[NE * 3 + 2] = DIR_27_NE_Z;
    
    d27.dirs[SW * 3    ] = DIR_27_SW_X;
    d27.dirs[SW * 3 + 1] = DIR_27_SW_Y;
    d27.dirs[SW * 3 + 2] = DIR_27_SW_Z;

    d27.dirs[SE * 3    ] = DIR_27_SE_X;
    d27.dirs[SE * 3 + 1] = DIR_27_SE_Y;
    d27.dirs[SE * 3 + 2] = DIR_27_SE_Z;

    d27.dirs[NW * 3    ] = DIR_27_NW_X;
    d27.dirs[NW * 3 + 1] = DIR_27_NW_Y;
    d27.dirs[NW * 3 + 2] = DIR_27_NW_Z;

    d27.dirs[TE * 3    ] = DIR_27_TE_X;
    d27.dirs[TE * 3 + 1] = DIR_27_TE_Y;
    d27.dirs[TE * 3 + 2] = DIR_27_TE_Z;

    d27.dirs[BW * 3    ] = DIR_27_BW_X;
    d27.dirs[BW * 3 + 1] = DIR_27_BW_Y;
    d27.dirs[BW * 3 + 2] = DIR_27_BW_Z;
                              
    d27.dirs[BE * 3    ] = DIR_27_BE_X;
    d27.dirs[BE * 3 + 1] = DIR_27_BE_Y;
    d27.dirs[BE * 3 + 2] = DIR_27_BE_Z;
                              
    d27.dirs[TW * 3    ] = DIR_27_TW_X;
    d27.dirs[TW * 3 + 1] = DIR_27_TW_Y;
    d27.dirs[TW * 3 + 2] = DIR_27_TW_Z;
                              
    d27.dirs[TN * 3    ] = DIR_27_TN_X;
    d27.dirs[TN * 3 + 1] = DIR_27_TN_Y;
    d27.dirs[TN * 3 + 2] = DIR_27_TN_Z;
                              
    d27.dirs[BS * 3    ] = DIR_27_BS_X;
    d27.dirs[BS * 3 + 1] = DIR_27_BS_Y;
    d27.dirs[BS * 3 + 2] = DIR_27_BS_Z;
                              
    d27.dirs[BN * 3    ] = DIR_27_BN_X;
    d27.dirs[BN * 3 + 1] = DIR_27_BN_Y;
    d27.dirs[BN * 3 + 2] = DIR_27_BN_Z;

    d27.dirs[TS * 3    ] = DIR_27_TS_X;
    d27.dirs[TS * 3 + 1] = DIR_27_TS_Y;
    d27.dirs[TS * 3 + 2] = DIR_27_TS_Z;

    d27.dirs[REST * 3    ] = DIR_27_REST_X;   //
    d27.dirs[REST * 3 + 1] = DIR_27_REST_Y;   //  ZERO ELEMENT
    d27.dirs[REST * 3 + 2] = DIR_27_REST_Z;   //

    d27.dirs[TNE * 3    ] = DIR_27_TNE_X;
    d27.dirs[TNE * 3 + 1] = DIR_27_TNE_Y;
    d27.dirs[TNE * 3 + 2] = DIR_27_TNE_Z;

    d27.dirs[BNE * 3    ] = DIR_27_BNE_X;
    d27.dirs[BNE * 3 + 1] = DIR_27_BNE_Y;
    d27.dirs[BNE * 3 + 2] = DIR_27_BNE_Z;

    d27.dirs[TSE * 3    ] = DIR_27_TSE_X;
    d27.dirs[TSE * 3 + 1] = DIR_27_TSE_Y;
    d27.dirs[TSE * 3 + 2] = DIR_27_TSE_Z;

    d27.dirs[BSE * 3    ] = DIR_27_BSE_X;
    d27.dirs[BSE * 3 + 1] = DIR_27_BSE_Y;
    d27.dirs[BSE * 3 + 2] = DIR_27_BSE_Z;

    d27.dirs[TNW * 3    ] = DIR_27_TNW_X;
    d27.dirs[TNW * 3 + 1] = DIR_27_TNW_Y;
    d27.dirs[TNW * 3 + 2] = DIR_27_TNW_Z;

    d27.dirs[BNW * 3    ] = DIR_27_BNW_X;
    d27.dirs[BNW * 3 + 1] = DIR_27_BNW_Y;
    d27.dirs[BNW * 3 + 2] = DIR_27_BNW_Z;

    d27.dirs[TSW * 3    ] = DIR_27_TSW_X;
    d27.dirs[TSW * 3 + 1] = DIR_27_TSW_Y;
    d27.dirs[TSW * 3 + 2] = DIR_27_TSW_Z;

    d27.dirs[BSW * 3    ] = DIR_27_BSW_X;
    d27.dirs[BSW * 3 + 1] = DIR_27_BSW_Y;
    d27.dirs[BSW * 3 + 2] = DIR_27_BSW_Z;

    return d27;
}

Distribution DistributionHelper::getDistribution(std::string name)
{
    if (name == "D3Q27")
        return getDistribution27();
    
    printf("wrong Distribution name... D3Q27 is chosen!\n");
    return getDistribution27();

}

std::vector<std::vector<real> > DistributionHelper::getQsWithoutRowsWithOnlyZeroValues(const Grid &grid, const Distribution &d)
{
    return getVectorWithoutRowsWithOnlyZeroValues(getAllQsOnFluidNodes(grid, d));
}

std::vector<std::vector<real> > DistributionHelper::getAllQsOnFluidNodes(const Grid &grid, const Distribution &d)
{
    std::vector<std::vector<real> > qs(grid.getSize(), std::vector<real>(d.dir_end + 1, 0));
    for (unsigned int node = 0; node < grid.getSize(); node++) {
        qs[node][0] = (real)node;
        for (int i = d.dir_start; i < d.dir_end; i++) {
            real qVal = (d.f)[i * grid.getSize() + node];
            if (qVal == 0.0f)
                continue;
            int index_fluidNode = getNeighborNodeIndexInGivenDirection(d, grid, node, i);
            qs[index_fluidNode][i + 1] = qVal;
        }
    }
    return qs;
}

int DistributionHelper::getNeighborNodeIndexInGivenDirection(const Distribution &d, const Grid &grid, const int node, const int dir_index)
{
    Vertex dir = Vertex((real)d.dirs[dir_index * DIMENSION + 0], (real)d.dirs[dir_index * DIMENSION + 1], (real)d.dirs[dir_index * DIMENSION + 2]);
    real x, y, z;
    grid.transIndexToCoords(node, x, y, z);
    Vertex solid_node = Vertex(x,y,z);
    Vertex fluid_node = Vertex(solid_node.x - dir.x, solid_node.y - dir.y, solid_node.z - dir.z);
    return grid.transCoordToIndex(fluid_node.x, fluid_node.y, fluid_node.z);
}

std::vector<std::vector<real> > DistributionHelper::getVectorWithoutRowsWithOnlyZeroValues(std::vector<std::vector<real> > qs)
{
    std::vector<std::vector<real> > qs_ausgeduennt;
    bool hasQs = false;
    for (std::size_t node = 0; node < qs.size(); node++) {
        for (std::size_t dir = 0; dir < qs[node].size() - 1; dir++) {
            if (qs[node][dir + 1] != 0)
                hasQs = true;
        }
        if (hasQs) {
            std::vector<real> qRow(qs[node].begin(), qs[node].end());
            qs_ausgeduennt.push_back(qRow);
            hasQs = false;
        }
    }
    return qs_ausgeduennt;
}

void DistributionHelper::printQs(std::vector<std::vector<real> > qs, int decimalPlaces)
{
    for (std::size_t node = 0; node < qs.size(); node++) {
        printf("index %zu: ", node);
        for (std::size_t dir = 1; dir < qs[node].size(); dir++) {
            printf("%d ", (int)qs[node][0]);
            printf("%.*f ", decimalPlaces, qs[node][dir]);
        }
        printf("\n");
    }
}

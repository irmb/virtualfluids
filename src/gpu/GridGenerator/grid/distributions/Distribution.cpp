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
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//=======================================================================================
#include "Distribution.h"

#include <stdio.h>

#include "grid/distributions/D3Q27.h"
#include "grid/Grid.h"
#include "lbm/constants/D3Q27.h"
using namespace vf::lbm::dir;

Distribution DistributionHelper::getDistribution27() 
{
    Distribution d27;
    d27.name = "D3Q27";
    d27.dir_start = STARTDIR;
    d27.dir_end = ENDDIR;

    d27.dirs.resize((ENDDIR + 1) * DIMENSION);

    d27.directions.resize(ENDDIR + 1);
    d27.directions[dP00] = Direction(DIR_27_E_X, DIR_27_E_Y, DIR_27_E_Z);
    d27.directions[dM00] = Direction(DIR_27_W_X, DIR_27_W_Y, DIR_27_W_Z);
    d27.directions[d0P0] = Direction(DIR_27_N_X, DIR_27_N_Y, DIR_27_N_Z);
    d27.directions[d0M0] = Direction(DIR_27_S_X, DIR_27_S_Y, DIR_27_S_Z);
    d27.directions[d00P] = Direction(DIR_27_T_X, DIR_27_T_Y, DIR_27_T_Z);
    d27.directions[d00M] = Direction(DIR_27_B_X, DIR_27_B_Y, DIR_27_B_Z);

    d27.directions[dPP0] = Direction(DIR_27_NE_X, DIR_27_NE_Y, DIR_27_NE_Z);
    d27.directions[dMM0] = Direction(DIR_27_SW_X, DIR_27_SW_Y, DIR_27_SW_Z);
    d27.directions[dPM0] = Direction(DIR_27_SE_X, DIR_27_SE_Y, DIR_27_SE_Z);
    d27.directions[dMP0] = Direction(DIR_27_NW_X, DIR_27_NW_Y, DIR_27_NW_Z);

    d27.directions[dP0P] = Direction(DIR_27_TE_X, DIR_27_TE_Y, DIR_27_TE_Z);
    d27.directions[dM0M] = Direction(DIR_27_BW_X, DIR_27_BW_Y, DIR_27_BW_Z);
    d27.directions[dP0M] = Direction(DIR_27_BE_X, DIR_27_BE_Y, DIR_27_BE_Z);
    d27.directions[dM0P] = Direction(DIR_27_TW_X, DIR_27_TW_Y, DIR_27_TW_Z);

    d27.directions[d0PP] = Direction(DIR_27_TN_X, DIR_27_TN_Y, DIR_27_TN_Z);
    d27.directions[d0MM] = Direction(DIR_27_BS_X, DIR_27_BS_Y, DIR_27_BS_Z);
    d27.directions[d0PM] = Direction(DIR_27_BN_X, DIR_27_BN_Y, DIR_27_BN_Z);
    d27.directions[d0MP] = Direction(DIR_27_TS_X, DIR_27_TS_Y, DIR_27_TS_Z);

    d27.directions[d000] = Direction(DIR_27_REST_X, DIR_27_REST_Y, DIR_27_REST_Z);

    d27.directions[dPPP] = Direction(DIR_27_TNE_X, DIR_27_TNE_Y, DIR_27_TNE_Z);
    d27.directions[dMPP] = Direction(DIR_27_TNW_X, DIR_27_TNW_Y, DIR_27_TNW_Z);
    d27.directions[dPMP] = Direction(DIR_27_TSE_X, DIR_27_TSE_Y, DIR_27_TSE_Z);
    d27.directions[dMMP] = Direction(DIR_27_TSW_X, DIR_27_TSW_Y, DIR_27_TSW_Z);

    d27.directions[dPPM] = Direction(DIR_27_BNE_X, DIR_27_BNE_Y, DIR_27_BNE_Z);
    d27.directions[dMPM]= Direction(DIR_27_BNW_X, DIR_27_BNW_Y, DIR_27_BNW_Z);
    d27.directions[dPMM]= Direction(DIR_27_BSE_X, DIR_27_BSE_Y, DIR_27_BSE_Z);
    d27.directions[dMMM] = Direction(DIR_27_BSW_X, DIR_27_BSW_Y, DIR_27_BSW_Z);


    d27.dirs[dP00 * 3    ] = DIR_27_E_X;
    d27.dirs[dP00 * 3 + 1] = DIR_27_E_Y;
    d27.dirs[dP00 * 3 + 2] = DIR_27_E_Z;

    d27.dirs[dM00 * 3    ] = DIR_27_W_X;
    d27.dirs[dM00 * 3 + 1] = DIR_27_W_Y;
    d27.dirs[dM00 * 3 + 2] = DIR_27_W_Z;
    
    d27.dirs[d0P0 * 3    ] = DIR_27_N_X;
    d27.dirs[d0P0 * 3 + 1] = DIR_27_N_Y;
    d27.dirs[d0P0 * 3 + 2] = DIR_27_N_Z;

    d27.dirs[d0M0 * 3    ] = DIR_27_S_X;
    d27.dirs[d0M0 * 3 + 1] = DIR_27_S_Y;
    d27.dirs[d0M0 * 3 + 2] = DIR_27_S_Z;
    
    d27.dirs[d00P * 3    ] = DIR_27_T_X;
    d27.dirs[d00P * 3 + 1] = DIR_27_T_Y;
    d27.dirs[d00P * 3 + 2] = DIR_27_T_Z;
    
    d27.dirs[d00M * 3    ] = DIR_27_B_X;
    d27.dirs[d00M * 3 + 1] = DIR_27_B_Y;
    d27.dirs[d00M * 3 + 2] = DIR_27_B_Z;

    d27.dirs[dPP0 * 3    ] = DIR_27_NE_X;
    d27.dirs[dPP0 * 3 + 1] = DIR_27_NE_Y;
    d27.dirs[dPP0 * 3 + 2] = DIR_27_NE_Z;
    
    d27.dirs[dMM0 * 3    ] = DIR_27_SW_X;
    d27.dirs[dMM0 * 3 + 1] = DIR_27_SW_Y;
    d27.dirs[dMM0 * 3 + 2] = DIR_27_SW_Z;

    d27.dirs[dPM0 * 3    ] = DIR_27_SE_X;
    d27.dirs[dPM0 * 3 + 1] = DIR_27_SE_Y;
    d27.dirs[dPM0 * 3 + 2] = DIR_27_SE_Z;

    d27.dirs[dMP0 * 3    ] = DIR_27_NW_X;
    d27.dirs[dMP0 * 3 + 1] = DIR_27_NW_Y;
    d27.dirs[dMP0 * 3 + 2] = DIR_27_NW_Z;

    d27.dirs[dP0P * 3    ] = DIR_27_TE_X;
    d27.dirs[dP0P * 3 + 1] = DIR_27_TE_Y;
    d27.dirs[dP0P * 3 + 2] = DIR_27_TE_Z;

    d27.dirs[dM0M * 3    ] = DIR_27_BW_X;
    d27.dirs[dM0M * 3 + 1] = DIR_27_BW_Y;
    d27.dirs[dM0M * 3 + 2] = DIR_27_BW_Z;
                              
    d27.dirs[dP0M * 3    ] = DIR_27_BE_X;
    d27.dirs[dP0M * 3 + 1] = DIR_27_BE_Y;
    d27.dirs[dP0M * 3 + 2] = DIR_27_BE_Z;
                              
    d27.dirs[dM0P * 3    ] = DIR_27_TW_X;
    d27.dirs[dM0P * 3 + 1] = DIR_27_TW_Y;
    d27.dirs[dM0P * 3 + 2] = DIR_27_TW_Z;
                              
    d27.dirs[d0PP * 3    ] = DIR_27_TN_X;
    d27.dirs[d0PP * 3 + 1] = DIR_27_TN_Y;
    d27.dirs[d0PP * 3 + 2] = DIR_27_TN_Z;
                              
    d27.dirs[d0MM * 3    ] = DIR_27_BS_X;
    d27.dirs[d0MM * 3 + 1] = DIR_27_BS_Y;
    d27.dirs[d0MM * 3 + 2] = DIR_27_BS_Z;
                              
    d27.dirs[d0PM * 3    ] = DIR_27_BN_X;
    d27.dirs[d0PM * 3 + 1] = DIR_27_BN_Y;
    d27.dirs[d0PM * 3 + 2] = DIR_27_BN_Z;

    d27.dirs[d0MP * 3    ] = DIR_27_TS_X;
    d27.dirs[d0MP * 3 + 1] = DIR_27_TS_Y;
    d27.dirs[d0MP * 3 + 2] = DIR_27_TS_Z;

    d27.dirs[d000 * 3    ] = DIR_27_REST_X;   //
    d27.dirs[d000 * 3 + 1] = DIR_27_REST_Y;   //  ZERO ELEMENT
    d27.dirs[d000 * 3 + 2] = DIR_27_REST_Z;   //

    d27.dirs[dPPP * 3    ] = DIR_27_TNE_X;
    d27.dirs[dPPP * 3 + 1] = DIR_27_TNE_Y;
    d27.dirs[dPPP * 3 + 2] = DIR_27_TNE_Z;

    d27.dirs[dPPM * 3    ] = DIR_27_BNE_X;
    d27.dirs[dPPM * 3 + 1] = DIR_27_BNE_Y;
    d27.dirs[dPPM * 3 + 2] = DIR_27_BNE_Z;

    d27.dirs[dPMP * 3    ] = DIR_27_TSE_X;
    d27.dirs[dPMP * 3 + 1] = DIR_27_TSE_Y;
    d27.dirs[dPMP * 3 + 2] = DIR_27_TSE_Z;

    d27.dirs[dPMM * 3    ] = DIR_27_BSE_X;
    d27.dirs[dPMM * 3 + 1] = DIR_27_BSE_Y;
    d27.dirs[dPMM * 3 + 2] = DIR_27_BSE_Z;

    d27.dirs[dMPP * 3    ] = DIR_27_TNW_X;
    d27.dirs[dMPP * 3 + 1] = DIR_27_TNW_Y;
    d27.dirs[dMPP * 3 + 2] = DIR_27_TNW_Z;

    d27.dirs[dMPM * 3    ] = DIR_27_BNW_X;
    d27.dirs[dMPM * 3 + 1] = DIR_27_BNW_Y;
    d27.dirs[dMPM * 3 + 2] = DIR_27_BNW_Z;

    d27.dirs[dMMP * 3    ] = DIR_27_TSW_X;
    d27.dirs[dMMP * 3 + 1] = DIR_27_TSW_Y;
    d27.dirs[dMMP * 3 + 2] = DIR_27_TSW_Z;

    d27.dirs[dMMM * 3    ] = DIR_27_BSW_X;
    d27.dirs[dMMM * 3 + 1] = DIR_27_BSW_Y;
    d27.dirs[dMMM * 3 + 2] = DIR_27_BSW_Z;

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

//! \}

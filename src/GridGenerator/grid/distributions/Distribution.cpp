#include "Distribution.h"
  
#include "D3Q7.h"
#include "D3Q13.h"
#include "D3Q19.h"
#include "D3Q27.h"
#include <stdio.h>

#include <GridGenerator/grid/Grid.cuh>

#include <GridGenerator/geometries/Vertex/Vertex.cuh>

Distribution DistributionHelper::getDistribution7() 
{
    Distribution d7;
    d7.name = "D3Q7";
    d7.dir_start = DIR_7_START;
    d7.dir_end = DIR_7_END;

    d7.dirs = new int[DIR_7_END * DIMENSION];

    int dir_num = 0;
    d7.dirs[dir_num++] = DIR_7_E_X;
    d7.dirs[dir_num++] = DIR_7_E_Y;
    d7.dirs[dir_num++] = DIR_7_E_Z;

    d7.dirs[dir_num++] = DIR_7_W_X;
    d7.dirs[dir_num++] = DIR_7_W_Y;
    d7.dirs[dir_num++] = DIR_7_W_Z;
          
    d7.dirs[dir_num++] = DIR_7_N_X;
    d7.dirs[dir_num++] = DIR_7_N_Y;
    d7.dirs[dir_num++] = DIR_7_N_Z;
        
    d7.dirs[dir_num++] = DIR_7_S_X;
    d7.dirs[dir_num++] = DIR_7_S_Y;
    d7.dirs[dir_num++] = DIR_7_S_Z;
          
    d7.dirs[dir_num++] = DIR_7_T_X;
    d7.dirs[dir_num++] = DIR_7_T_Y;
    d7.dirs[dir_num++] = DIR_7_T_Z;
           
    d7.dirs[dir_num++] = DIR_7_B_X;
    d7.dirs[dir_num++] = DIR_7_B_Y;
    d7.dirs[dir_num++] = DIR_7_B_Z;

    return d7;
}

Distribution DistributionHelper::getDistribution13() 
{
    Distribution d13;
    d13.name = "D3Q13";
    d13.dir_start = DIR_13_START;
    d13.dir_end = DIR_13_END;

    d13.dirs = new int[DIR_13_END * DIMENSION];

    int dir_num = 0;
    d13.dirs[dir_num++] = DIR_13_NE_X;
    d13.dirs[dir_num++] = DIR_13_NE_Y;
    d13.dirs[dir_num++] = DIR_13_NE_Z;

    d13.dirs[dir_num++] = DIR_13_SW_X;
    d13.dirs[dir_num++] = DIR_13_SW_Y;
    d13.dirs[dir_num++] = DIR_13_SW_Z;

    d13.dirs[dir_num++] = DIR_13_SE_X;
    d13.dirs[dir_num++] = DIR_13_SE_Y;
    d13.dirs[dir_num++] = DIR_13_SE_Z;

    d13.dirs[dir_num++] = DIR_13_NW_X;
    d13.dirs[dir_num++] = DIR_13_NW_Y;
    d13.dirs[dir_num++] = DIR_13_NW_Z;

    d13.dirs[dir_num++] = DIR_13_TE_X;
    d13.dirs[dir_num++] = DIR_13_TE_Y;
    d13.dirs[dir_num++] = DIR_13_TE_Z;

    d13.dirs[dir_num++] = DIR_13_BW_X;
    d13.dirs[dir_num++] = DIR_13_BW_Y;
    d13.dirs[dir_num++] = DIR_13_BW_Z;

    d13.dirs[dir_num++] = DIR_13_BE_X;
    d13.dirs[dir_num++] = DIR_13_BE_Y;
    d13.dirs[dir_num++] = DIR_13_BE_Z;

    d13.dirs[dir_num++] = DIR_13_TW_X;
    d13.dirs[dir_num++] = DIR_13_TW_Y;
    d13.dirs[dir_num++] = DIR_13_TW_Z;

    d13.dirs[dir_num++] = DIR_13_TN_X;
    d13.dirs[dir_num++] = DIR_13_TN_Y;
    d13.dirs[dir_num++] = DIR_13_TN_Z;

    d13.dirs[dir_num++] = DIR_13_BS_X;
    d13.dirs[dir_num++] = DIR_13_BS_Y;
    d13.dirs[dir_num++] = DIR_13_BS_Z;

    d13.dirs[dir_num++] = DIR_13_BN_X;
    d13.dirs[dir_num++] = DIR_13_BN_Y;
    d13.dirs[dir_num++] = DIR_13_BN_Z;

    d13.dirs[dir_num++] = DIR_13_TS_X;
    d13.dirs[dir_num++] = DIR_13_TS_Y;
    d13.dirs[dir_num++] = DIR_13_TS_Z;

    return d13;
}

Distribution DistributionHelper::getDistribution19() 
{
    Distribution d19;
    d19.name = "D3Q19";
    d19.dir_start = DIR_19_START;
    d19.dir_end = DIR_19_END;

    d19.dirs = new int[DIR_19_END * DIMENSION];

    int dir_num = 0;
    d19.dirs[dir_num++] = DIR_19_E_X;
    d19.dirs[dir_num++] = DIR_19_E_Y;
    d19.dirs[dir_num++] = DIR_19_E_Z;

    d19.dirs[dir_num++] = DIR_19_W_X;
    d19.dirs[dir_num++] = DIR_19_W_Y;
    d19.dirs[dir_num++] = DIR_19_W_Z;

    d19.dirs[dir_num++] = DIR_19_N_X;
    d19.dirs[dir_num++] = DIR_19_N_Y;
    d19.dirs[dir_num++] = DIR_19_N_Z;

    d19.dirs[dir_num++] = DIR_19_S_X;
    d19.dirs[dir_num++] = DIR_19_S_Y;
    d19.dirs[dir_num++] = DIR_19_S_Z;

    d19.dirs[dir_num++] = DIR_19_T_X;
    d19.dirs[dir_num++] = DIR_19_T_Y;
    d19.dirs[dir_num++] = DIR_19_T_Z;

    d19.dirs[dir_num++] = DIR_19_B_X;
    d19.dirs[dir_num++] = DIR_19_B_Y;
    d19.dirs[dir_num++] = DIR_19_B_Z;

    d19.dirs[dir_num++] = DIR_19_NE_X;
    d19.dirs[dir_num++] = DIR_19_NE_Y;
    d19.dirs[dir_num++] = DIR_19_NE_Z;

    d19.dirs[dir_num++] = DIR_19_SW_X;
    d19.dirs[dir_num++] = DIR_19_SW_Y;
    d19.dirs[dir_num++] = DIR_19_SW_Z;

    d19.dirs[dir_num++] = DIR_19_SE_X;
    d19.dirs[dir_num++] = DIR_19_SE_Y;
    d19.dirs[dir_num++] = DIR_19_SE_Z;

    d19.dirs[dir_num++] = DIR_19_NW_X;
    d19.dirs[dir_num++] = DIR_19_NW_Y;
    d19.dirs[dir_num++] = DIR_19_NW_Z;

    d19.dirs[dir_num++] = DIR_19_TE_X;
    d19.dirs[dir_num++] = DIR_19_TE_Y;
    d19.dirs[dir_num++] = DIR_19_TE_Z;

    d19.dirs[dir_num++] = DIR_19_BW_X;
    d19.dirs[dir_num++] = DIR_19_BW_Y;
    d19.dirs[dir_num++] = DIR_19_BW_Z;

    d19.dirs[dir_num++] = DIR_19_BE_X;
    d19.dirs[dir_num++] = DIR_19_BE_Y;
    d19.dirs[dir_num++] = DIR_19_BE_Z;
    
    d19.dirs[dir_num++] = DIR_19_TW_X;
    d19.dirs[dir_num++] = DIR_19_TW_Y;
    d19.dirs[dir_num++] = DIR_19_TW_Z;
    
    d19.dirs[dir_num++] = DIR_19_TN_X;
    d19.dirs[dir_num++] = DIR_19_TN_Y;
    d19.dirs[dir_num++] = DIR_19_TN_Z;
    
    d19.dirs[dir_num++] = DIR_19_BS_X;
    d19.dirs[dir_num++] = DIR_19_BS_Y;
    d19.dirs[dir_num++] = DIR_19_BS_Z;
    
    d19.dirs[dir_num++] = DIR_19_BN_X;
    d19.dirs[dir_num++] = DIR_19_BN_Y;
    d19.dirs[dir_num++] = DIR_19_BN_Z;
    
    d19.dirs[dir_num++] = DIR_19_TS_X;
    d19.dirs[dir_num++] = DIR_19_TS_Y;
    d19.dirs[dir_num++] = DIR_19_TS_Z;

    return d19;
}

Distribution DistributionHelper::getDistribution27() 
{
    Distribution d27;
    d27.name = "D3Q27";
    d27.dir_start = DIR_27_START;
    d27.dir_end = DIR_27_END;

    d27.dirs = new int[(DIR_27_END + 1) * DIMENSION];

    int dir_num = 0;
    d27.dirs[dir_num++] = DIR_19_E_X;
    d27.dirs[dir_num++] = DIR_19_E_Y;
    d27.dirs[dir_num++] = DIR_19_E_Z;

    d27.dirs[dir_num++] = DIR_19_W_X;
    d27.dirs[dir_num++] = DIR_19_W_Y;
    d27.dirs[dir_num++] = DIR_19_W_Z;

    d27.dirs[dir_num++] = DIR_19_N_X;
    d27.dirs[dir_num++] = DIR_19_N_Y;
    d27.dirs[dir_num++] = DIR_19_N_Z;

    d27.dirs[dir_num++] = DIR_19_S_X;
    d27.dirs[dir_num++] = DIR_19_S_Y;
    d27.dirs[dir_num++] = DIR_19_S_Z;

    d27.dirs[dir_num++] = DIR_19_T_X;
    d27.dirs[dir_num++] = DIR_19_T_Y;
    d27.dirs[dir_num++] = DIR_19_T_Z;

    d27.dirs[dir_num++] = DIR_19_B_X;
    d27.dirs[dir_num++] = DIR_19_B_Y;
    d27.dirs[dir_num++] = DIR_19_B_Z;

    d27.dirs[dir_num++] = DIR_19_NE_X;
    d27.dirs[dir_num++] = DIR_19_NE_Y;
    d27.dirs[dir_num++] = DIR_19_NE_Z;

    d27.dirs[dir_num++] = DIR_19_SW_X;
    d27.dirs[dir_num++] = DIR_19_SW_Y;
    d27.dirs[dir_num++] = DIR_19_SW_Z;

    d27.dirs[dir_num++] = DIR_19_SE_X;
    d27.dirs[dir_num++] = DIR_19_SE_Y;
    d27.dirs[dir_num++] = DIR_19_SE_Z;

    d27.dirs[dir_num++] = DIR_19_NW_X;
    d27.dirs[dir_num++] = DIR_19_NW_Y;
    d27.dirs[dir_num++] = DIR_19_NW_Z;

    d27.dirs[dir_num++] = DIR_19_TE_X;
    d27.dirs[dir_num++] = DIR_19_TE_Y;
    d27.dirs[dir_num++] = DIR_19_TE_Z;

    d27.dirs[dir_num++] = DIR_19_BW_X;
    d27.dirs[dir_num++] = DIR_19_BW_Y;
    d27.dirs[dir_num++] = DIR_19_BW_Z;

    d27.dirs[dir_num++] = DIR_19_BE_X;
    d27.dirs[dir_num++] = DIR_19_BE_Y;
    d27.dirs[dir_num++] = DIR_19_BE_Z;

    d27.dirs[dir_num++] = DIR_19_TW_X;
    d27.dirs[dir_num++] = DIR_19_TW_Y;
    d27.dirs[dir_num++] = DIR_19_TW_Z;

    d27.dirs[dir_num++] = DIR_19_TN_X;
    d27.dirs[dir_num++] = DIR_19_TN_Y;
    d27.dirs[dir_num++] = DIR_19_TN_Z;

    d27.dirs[dir_num++] = DIR_19_BS_X;
    d27.dirs[dir_num++] = DIR_19_BS_Y;
    d27.dirs[dir_num++] = DIR_19_BS_Z;

    d27.dirs[dir_num++] = DIR_19_BN_X;
    d27.dirs[dir_num++] = DIR_19_BN_Y;
    d27.dirs[dir_num++] = DIR_19_BN_Z;

    d27.dirs[dir_num++] = DIR_19_TS_X;
    d27.dirs[dir_num++] = DIR_19_TS_Y;
    d27.dirs[dir_num++] = DIR_19_TS_Z;


    d27.dirs[dir_num++] = 0;   //
    d27.dirs[dir_num++] = 0;   //  ZERO ELEMENT
    d27.dirs[dir_num++] = 0;   //


    d27.dirs[dir_num++] = DIR_27_TNE_X;
    d27.dirs[dir_num++] = DIR_27_TNE_Y;
    d27.dirs[dir_num++] = DIR_27_TNE_Z;

    d27.dirs[dir_num++] = DIR_27_BNE_X;
    d27.dirs[dir_num++] = DIR_27_BNE_Y;
    d27.dirs[dir_num++] = DIR_27_BNE_Z;

    d27.dirs[dir_num++] = DIR_27_TSE_X;
    d27.dirs[dir_num++] = DIR_27_TSE_Y;
    d27.dirs[dir_num++] = DIR_27_TSE_Z;

    d27.dirs[dir_num++] = DIR_27_BSE_X;
    d27.dirs[dir_num++] = DIR_27_BSE_Y;
    d27.dirs[dir_num++] = DIR_27_BSE_Z;


    d27.dirs[dir_num++] = DIR_27_TNW_X;
    d27.dirs[dir_num++] = DIR_27_TNW_Y;
    d27.dirs[dir_num++] = DIR_27_TNW_Z;

    d27.dirs[dir_num++] = DIR_27_BNW_X;
    d27.dirs[dir_num++] = DIR_27_BNW_Y;
    d27.dirs[dir_num++] = DIR_27_BNW_Z;

    d27.dirs[dir_num++] = DIR_27_TSW_X;
    d27.dirs[dir_num++] = DIR_27_TSW_Y;
    d27.dirs[dir_num++] = DIR_27_TSW_Z;

    d27.dirs[dir_num++] = DIR_27_BSW_X;
    d27.dirs[dir_num++] = DIR_27_BSW_Y;
    d27.dirs[dir_num++] = DIR_27_BSW_Z;

    return d27;
}

Distribution DistributionHelper::getDistribution(std::string name)
{
    if (name == "D3Q7")
        return getDistribution7();
    if (name == "D3Q13")
        return getDistribution13();
    if (name == "D3Q19")
        return getDistribution19();
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
    std::vector<std::vector<real> > qs(grid.size, std::vector<real>(d.dir_end + 1, 0));
    for (unsigned int node = 0; node < grid.size; node++) {
        qs[node][0] = (real)node;
        for (int i = d.dir_start; i < d.dir_end; i++) {
            real qVal = (d.f)[i * grid.size + node];
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
    unsigned int x, y, z;
    grid.transIndexToCoords(node, x, y, z);
    Vertex solid_node = Vertex((real)x, (real)y,(real)z);
    Vertex fluid_node = Vertex(solid_node.x - dir.x, solid_node.y - dir.y, solid_node.z - dir.z);
    return grid.transCoordToIndex(fluid_node);
}

std::vector<std::vector<real> > DistributionHelper::getVectorWithoutRowsWithOnlyZeroValues(std::vector<std::vector<real> > qs)
{
    std::vector<std::vector<real> > qs_ausgeduennt;
    bool hasQs = false;
    for (int node = 0; node < qs.size(); node++) {
        for (int dir = 0; dir < qs[node].size() - 1; dir++) {
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
    for (int node = 0; node < qs.size(); node++) {
        printf("index %d: ", node);
        for (int dir = 1; dir < qs[node].size(); dir++) {
            printf("%d ", (int)qs[node][0]);
            printf("%.*f ", decimalPlaces, qs[node][dir]);
        }
        printf("\n");
    }
}

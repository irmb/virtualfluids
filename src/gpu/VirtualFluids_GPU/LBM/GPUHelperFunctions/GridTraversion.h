#ifndef GRID_TRAVERSION
#define GRID_TRAVERSION

#include <basics/DataTypes.h>
#include <math.h>

__inline__ __host__ __device__ uint traverseSourceCell(real coordDestinationX, real coordDestinationY, real coordDestinationZ,
                                   uint indexOfInverseNeighborOfSourceCell, real coordInvNeighborCellX, real coordInvNeighborCellY,
                                   real coordInvNeighborCellZ, const uint *neighborX, const uint *neighborY, const uint *neighborZ, real dx)
{
    uint xTraverse = (uint)floor((coordDestinationX - coordInvNeighborCellX) / dx);
    uint yTraverse = (uint)floor((coordDestinationY - coordInvNeighborCellY) / dx);
    uint zTraverse = (uint)floor((coordDestinationZ - coordInvNeighborCellZ) / dx);

    // printf("s %.4f, d %.4f, delta %.4f, dx %.4f, xTraverse %d\n", coordInvNeighborCellX, coordDestinationX, coordDestinationX - coordInvNeighborCellX, dx, xTraverse);

    uint newIndexOfSourceCell = indexOfInverseNeighborOfSourceCell;
    for (uint ix = 1; ix <= xTraverse; ix++) {
        newIndexOfSourceCell = neighborX[newIndexOfSourceCell];
    }
    for (uint iy = 1; iy <= yTraverse; iy++) {
        newIndexOfSourceCell = neighborY[newIndexOfSourceCell];
    }
    for (uint iz = 1; iz <= zTraverse; iz++) {
        newIndexOfSourceCell = neighborZ[newIndexOfSourceCell];
    }
    return newIndexOfSourceCell;
}

#endif
#ifndef FIND_NEIGHBORS_H
#define FIND_NEIGHBORS_H

#include "Parameter/Parameter.h"
#include "lbm/constants/D3Q27.h"

using namespace vf::lbm::dir;

// Only use for fluid nodes!
inline uint getNeighborIndex(LBMSimulationParameter *parH, const uint position, const int direction)
{
    uint nodeIndex = position;

    if (mapForPointerChasing.at(direction).counterInverse != 0) {
        nodeIndex = parH->neighborInverse[nodeIndex];
    }

    for (uint x = 0; x < mapForPointerChasing.at(direction).counterX; x++) {
        nodeIndex = parH->neighborX[nodeIndex];
    }

    for (uint y = 0; y < mapForPointerChasing.at(direction).counterY; y++) {
        nodeIndex = parH->neighborY[nodeIndex];
    }

    for (uint z = 0; z < mapForPointerChasing.at(direction).counterZ; z++) {
        nodeIndex = parH->neighborZ[nodeIndex];
    }

    return nodeIndex;
}

#endif

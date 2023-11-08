#ifndef FIND_NEIGHBORS_H
#define FIND_NEIGHBORS_H

#include <map>

#include "Parameter/Parameter.h"
#include "lbm/constants/D3Q27.h"

using namespace vf::lbm::dir;

struct countersForPointerChasing
{
    uint counterInverse;
    uint counterX;
    uint counterY;
    uint counterZ;
};

const std::map<const size_t, const countersForPointerChasing> mapForPointerChasing =
{
    {d000, countersForPointerChasing{0, 0, 0, 0}},
    {dP00, countersForPointerChasing{0, 1, 0, 0}},
    {dM00, countersForPointerChasing{1, 0, 1, 1}},
    {d0P0, countersForPointerChasing{0, 0, 1, 0}},
    {d0M0, countersForPointerChasing{1, 1, 0, 1}},
    {d00P, countersForPointerChasing{0, 0, 0, 1}},
    {d00M, countersForPointerChasing{1, 1, 1, 0}},

    {dPP0, countersForPointerChasing{0, 1, 1, 0}},
    {dMM0, countersForPointerChasing{1, 0, 0, 1}},
    {dPM0, countersForPointerChasing{1, 2, 0, 1}},
    {dMP0, countersForPointerChasing{1, 0, 2, 1}},
    {dP0P, countersForPointerChasing{0, 1, 0, 1}},
    {dM0M, countersForPointerChasing{1, 0, 1, 0}},
    {dP0M, countersForPointerChasing{1, 2, 1, 0}},
    {dM0P, countersForPointerChasing{1, 0, 1, 2}},
    {d0PP, countersForPointerChasing{0, 0, 1, 1}},
    {d0MM, countersForPointerChasing{1, 1, 0, 0}},
    {d0PM, countersForPointerChasing{1, 1, 2, 0}},
    {d0MP, countersForPointerChasing{1, 1, 0, 2}},

    {dPPP, countersForPointerChasing{0, 1, 1, 1}},
    {dMPP, countersForPointerChasing{1, 0, 2, 2}},
    {dPMP, countersForPointerChasing{1, 2, 0, 2}},
    {dMMP, countersForPointerChasing{1, 0, 0, 2}},
    {dPPM, countersForPointerChasing{1, 2, 2, 0}},
    {dMPM, countersForPointerChasing{1, 0, 2, 0}},
    {dPMM, countersForPointerChasing{1, 2, 0, 0}},
    {dMMM, countersForPointerChasing{1, 0, 0, 0}}
};

// Only use for fluid nodes!
inline uint getNeighborIndex(LBMSimulationParameter* parH, const uint position, const int direction)
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

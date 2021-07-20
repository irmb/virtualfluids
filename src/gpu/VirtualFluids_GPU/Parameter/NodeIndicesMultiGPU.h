#ifndef INDICES_MULTIGPU_H
#define INDICES_MULTIGPU_H

#include <vector>
#include <memory>
#include "basics/Core/DataTypes.h"

class NodeIndicesMultiGPU
{
    const std::vector<uint> *geoFluidSizes;
    const std::vector<const uint*> *geoFluidNodeIndices; 

public:
    NodeIndicesMultiGPU(const std::vector<uint> *geoFluidSizes, const std::vector<const uint*> * geoFluidNodes);

    uint getGeoFluidSize(uint level);
    const uint *getGeoFluidNodeIndices(uint level);
};

#endif
#ifndef INDICES_MULTIGPU_H
#define INDICES_MULTIGPU_H

#include <vector>
#include <memory>
#include "basics/Core/DataTypes.h"

class NodeIndicesMultiGPU
{
    const std::vector<uint> *geoFluidSizes;
    const std::vector<const std::vector<uint>*> *geoFluidNodeIndices; 

public:
    NodeIndicesMultiGPU(const std::vector<uint> *geoFluidSizes, const std::vector<const std::vector<uint> *> * geoFluidNodes);

    uint getGeoFluidSize(uint gridNumber);
    const std::vector<uint>* getGeoFluidNodeIndices(uint gridNumber);
};

#endif
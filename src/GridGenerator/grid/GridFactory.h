#ifndef GRID_FACTORY_H
#define GRID_FACTORY_H

#include <VirtualFluidsDefinitions.h>
#include <core/PointerDefinitions.h>

#include "Grid.cuh"

class VF_PUBLIC GridFactory
{
public:
    static SPtr<Grid> makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta,
                               const std::string& device, const std::string& d3Qxx);
};


#endif

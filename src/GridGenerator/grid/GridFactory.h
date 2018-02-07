#ifndef GRID_FACTORY_H
#define GRID_FACTORY_H

#include <VirtualFluidsDefinitions.h>
#include <core/PointerDefinitions.h>

#include "grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h"
#include "grid/GridStrategy/GridGpuStrategy/GridGpuStrategy.h"
#include "distributions/Distribution.h"

enum class Device
{
    CPU, GPU
};



template <typename Grid>
class VF_PUBLIC GridFactory
{
public:
    SPtr<Grid> makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, const std::string& d3Qxx = "D3Q27")
    {
        if (!gridStrategy)
            throw "GridStrategy has to be set before make Grid!";

        Distribution distribution = DistributionHelper::getDistribution(d3Qxx);
        return Grid::makeShared(startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, distribution);
    }

    void setGridStrategy(Device device)
    {
        switch (device)
        {
        case Device::CPU:
            gridStrategy = SPtr<GridStrategy>(new GridCpuStrategy()); break;
        case Device::GPU:
            gridStrategy = SPtr<GridStrategy>(new GridGpuStrategy()); break;
        }
    }

    void setGridStrategy(SPtr<GridStrategy> gridStrategy)
    {
        this->gridStrategy = gridStrategy;
    }

private:
    SPtr<GridStrategy> gridStrategy;
};


#endif

#ifndef GRID_FACTORY_H
#define GRID_FACTORY_H

#include <VirtualFluidsDefinitions.h>
#include <core/PointerDefinitions.h>

#include "grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h"
#include "grid/GridStrategy/GridGpuStrategy/GridGpuStrategy.h"
#include "distributions/Distribution.h"
#include "GridImp.cuh"
#include "GridMocks.h"
#include "geometries/Cuboid/Cuboid.h"

enum class Device
{
    CPU, GPU
};


class VF_PUBLIC GridFactory
{
public:
    SPtr<Grid> makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, const std::string& d3Qxx = "D3Q27")
    {
        if (!gridStrategy)
            throw "GridStrategy has to be set before make Grid!";

        Distribution distribution = DistributionHelper::getDistribution(d3Qxx);
        if(this->grid == "stub")
            return GridStub::makeShared(startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, distribution);
        else if(this->grid == "spy")
            return GridSpy::makeShared(startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, distribution);


        return GridImp::makeShared(new Cuboid(startX, startY, startZ, endX, endY, endZ), delta, gridStrategy, distribution);
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

    void setGrid(const std::string& grid)
    {
        this->grid = grid;
    }

private:
    SPtr<GridStrategy> gridStrategy;
    std::string grid;
};


#endif

#include "GridFactory.h"

#include <iostream>

#include "grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h"
#include "grid/GridStrategy/GridGpuStrategy/GridGpuStrategy.h"
#include "distributions/Distribution.h"


SPtr<Grid> GridFactory::makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta,
                                 const std::string& device, const std::string& d3Qxx)
{
    SPtr<GridStrategy> gridStrategy;
    if (device == "gpu")
        gridStrategy = SPtr<GridStrategy>(new GridGpuStrategy());
    else if (device == "cpu")
        gridStrategy = SPtr<GridStrategy>(new GridCpuStrategy());
    else
    {
        std::cout << "Cannot identify device: " << device << ". cpu is chosen!\n";
        gridStrategy = SPtr<GridStrategy>(new GridGpuStrategy());
    }

    Distribution distribution = DistributionHelper::getDistribution(d3Qxx);

    return Grid::makeShared(startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, distribution);
}

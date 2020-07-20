#ifndef GRID_FACTORY_H
#define GRID_FACTORY_H

#include "global.h"

#include "geometries/Cuboid/Cuboid.h"
#include "geometries/Sphere/Sphere.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"

#include "grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h"
#include "grid/GridStrategy/GridGpuStrategy/GridGpuStrategy.h"
#include "grid/distributions/Distribution.h"
#include "grid/GridImp.h"

enum class Device
{
    CPU, GPU
};

enum class TriangularMeshDiscretizationMethod
{
    RAYCASTING, POINT_IN_OBJECT, POINT_UNDER_TRIANGLE
};

class VIRTUALFLUIDS_GPU_EXPORT GridFactory
{
public:
    static SPtr<GridFactory> make()
    {
        return SPtr<GridFactory>(new GridFactory());
    }

private:
    GridFactory()
    {
        gridStrategy = SPtr<GridStrategy>(new GridCpuStrategy());
        triangularMeshDiscretizationStrategy = new RayCastingDiscretizationStrategy();
    }

public:
    SPtr<Grid> makeGrid(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level, const std::string& d3Qxx = "D3Q27")
    {
        Distribution distribution = DistributionHelper::getDistribution(d3Qxx);

        SPtr<GridImp> grid;

        grid = GridImp::makeShared(gridShape, startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, distribution, level);

        grid->setTriangularMeshDiscretizationStrategy(this->triangularMeshDiscretizationStrategy);

        return grid;
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

    void setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod triangularMeshDiscretizationMethod)
    {
        switch (triangularMeshDiscretizationMethod)
        {
        case TriangularMeshDiscretizationMethod::POINT_UNDER_TRIANGLE:
            triangularMeshDiscretizationStrategy = new PointUnderTriangleStrategy();
            break;
        case TriangularMeshDiscretizationMethod::RAYCASTING:
            triangularMeshDiscretizationStrategy = new RayCastingDiscretizationStrategy();
            break;
        case TriangularMeshDiscretizationMethod::POINT_IN_OBJECT:
            triangularMeshDiscretizationStrategy = new PointInObjectDiscretizationStrategy();
            break;
        }
    }

private:
    TriangularMeshDiscretizationStrategy* triangularMeshDiscretizationStrategy;
    SPtr<GridStrategy> gridStrategy;
};


#endif

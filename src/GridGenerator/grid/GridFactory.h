#ifndef GRID_FACTORY_H
#define GRID_FACTORY_H

#include <VirtualFluidsDefinitions.h>
#include <core/PointerDefinitions.h>

#include "grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h"
#include "grid/GridStrategy/GridGpuStrategy/GridGpuStrategy.h"
#include "distributions/Distribution.h"
#include "GridImp.h"
#include "grid/gridmocks.h"
#include "geometries/Cuboid/Cuboid.h"
#include "geometries/Sphere/Sphere.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"

enum class Device
{
    CPU, GPU
};

enum class TriangularMeshDiscretizationMethod
{
    RAYCASTING, POINT_IN_OBJECT, POINT_UNDER_TRIANGLE
};

enum class TestDouble
{
    PRODUCTIONCLASS, STUB, SPY
};


class VF_PUBLIC GridFactory
{
public:
    static SPtr<GridFactory> make()
    {
        return SPtr<GridFactory>(new GridFactory());
    }

private:
    GridFactory()
    {
        gridType = TestDouble::PRODUCTIONCLASS;
        gridStrategy = SPtr<GridStrategy>(new GridCpuStrategy());
        triangularMeshDiscretizationStrategy = new RayCastingDiscretizationStrategy();
    }

public:
    SPtr<Grid> makeGrid(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level, const std::string& d3Qxx = "D3Q27")
    {
        Distribution distribution = DistributionHelper::getDistribution(d3Qxx);

        SPtr<GridImp> grid;

        switch (gridType)
        {
        case TestDouble::PRODUCTIONCLASS:
            grid = GridImp::makeShared(gridShape, startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, distribution, level);
            break;
        case TestDouble::STUB:
            return GridStub::makeShared(gridShape, startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, distribution, level);
        case TestDouble::SPY:
            return GridSpy::makeShared(gridShape, startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, distribution, level);
        }

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

    void setGridType(TestDouble gridType)
    {
        this->gridType = gridType;
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
    TestDouble gridType;
};


#endif

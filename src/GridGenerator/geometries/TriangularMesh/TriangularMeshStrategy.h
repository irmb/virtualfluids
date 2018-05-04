#ifndef TriangularMeshStrategy_H
#define TriangularMeshStrategy_H

#include "GridGenerator/global.h"

#include "../Triangle/Triangle.h"
#include "TriangularMesh.h"
#include "numerics/geometry3d/GbTriFaceMesh3D.h"
#include "grid/GridImp.h"
#include "utilities/logger/Logger.h"
#include "grid/NodeValues.h"


class TriangularMeshDiscretizationStrategy
{
public:
    TriangularMeshDiscretizationStrategy() {}
    virtual ~TriangularMeshDiscretizationStrategy() {}


    void discretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType)
    {  
        this->doDiscretize(triangularMesh, grid, InnerType, OuterType);
        this->removeOddBoundaryCellNodes(grid);
    }

private:
    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType) = 0;
    void removeOddBoundaryCellNodes(GridImp* grid)
    {
#pragma omp parallel for
        for (uint index = 0; index < grid->getSize(); index++)
            grid->removeOddBoundaryCellNode(index);
    }


};



class PointInObjectDiscretizationStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    PointInObjectDiscretizationStrategy() {}
    virtual ~PointInObjectDiscretizationStrategy() {}

    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType)
    {
        auto mesh = triangularMesh->getGbTriFaceMesh3D();

        for (uint i = 0; i < grid->getSize(); i++)
        {
            real x, y, z;
            grid->transIndexToCoords(i, x, y, z);

            if (mesh->isPointInGbObject3D(x, y, z))
                grid->setNodeTo(i, InnerType);
            //else
            //    grid->setNodeTo(i, OuterType);
        }

        delete mesh;
    }
};

class RayCastingDiscretizationStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    RayCastingDiscretizationStrategy() {}
    virtual ~RayCastingDiscretizationStrategy() {}

    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType)
    {
        auto mesh = triangularMesh->getGbTriFaceMesh3D();

        const real minXExact = triangularMesh->minmax.minX;
        const real minYExact = triangularMesh->minmax.minY;
        const real minZExact = triangularMesh->minmax.minZ;

        const real maxXExact = triangularMesh->minmax.maxX;
        const real maxYExact = triangularMesh->minmax.maxY;
        const real maxZExact = triangularMesh->minmax.maxZ;

        const auto min = grid->getMinimumOnNode(Vertex(minXExact, minYExact, minZExact));

        const real minX = min.x;
        const real minY = min.y;
        const real minZ = min.z;

        const auto max = grid->getMaximumOnNode(Vertex(maxXExact, maxYExact, maxZExact));

        const real maxX = max.x;
        const real maxY = max.y;
        const real maxZ = max.z;


        real x, y, z;
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (x = minX; x <= maxX; x += grid->getDelta())
                {
                    grid->setNodeTo(grid->transCoordToIndex(x, y, z), InnerType);
                }
            }
        }



        int counter = 0;

        // Test line intersection
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (x = minX; x <= maxX; x += grid->getDelta())
                {
                    counter++;
                    if (mesh->intersectLine((x - grid->getDelta()), y, z, x, y, z)) 
                        break;
                    grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                }
            }
        }

        // Test line intersection from opposite direction
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (x = maxX; x >= minX; x -= grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        counter++;
                        if (mesh->intersectLine((x + grid->getDelta()), y, z, x, y, z))
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        // Test line intersection
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (x = minX; x <= maxX; x += grid->getDelta())
            {
                for (y = minY; y <= maxY; y += grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        counter++;
                        if (mesh->intersectLine(x, (y - grid->getDelta()), z, x, y, z)) 
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        // Test line intersection from opposite direction
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (x = minX; x <= maxX; x += grid->getDelta())
            {
                for (y = maxY; y >= minY; y -= grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        counter++;
                        if (mesh->intersectLine(x, (y + grid->getDelta()), z, x, y, z))
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        // Test line intersection
        for (x = minX; x <= maxX; x += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (z = minZ; z <= maxZ; z += grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        counter++;
                        if (mesh->intersectLine(x, y, (z - grid->getDelta()), x, y, z)) 
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        // Test line intersection from opposite direction
        for (x = minX; x <= maxX; x += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (z = maxZ; z >= minZ; z -= grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        counter++;
                        if (mesh->intersectLine(x, y, (z + grid->getDelta()), x, y, z)) 
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        delete mesh;
    }
};

class PointUnderTriangleStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    PointUnderTriangleStrategy() {}
    virtual ~PointUnderTriangleStrategy() {}
    void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char innerType, char outerType) override
    {
#pragma omp parallel for
        for (long i = 0; i < triangularMesh->size; i++)
            this->meshReverse(triangularMesh->triangles[i], grid, innerType);

        this->findInsideNodes(grid, innerType);

#pragma omp parallel for
        for (uint i = 0; i < grid->getSize(); i++)
            this->setNegativeDirBorderTo(grid, i, innerType);
    }

private:
    HOSTDEVICE void meshReverse(Triangle &triangle, GridImp* grid, char innerType)
    {
        auto box = grid->getBoundingBoxOnNodes(triangle);

        const real delta = grid->getDelta();
        triangle.initalLayerThickness(delta);

        for (real x = box.minX; x <= box.maxX; x += delta)
        {
            for (real y = box.minY; y <= box.maxY; y += delta)
            {
                for (real z = box.minZ; z <= box.maxZ; z += delta)
                {
                    const uint index = grid->transCoordToIndex(x, y, z);

                    const Vertex point(x, y, z);

                    const char pointValue = triangle.isUnderFace(point);

                    if (pointValue == NEGATIVE_DIRECTION_BORDER)
                        grid->setNodeTo(index, NEGATIVE_DIRECTION_BORDER);
                    else if (pointValue == INSIDE)
                        grid->setNodeTo(index, innerType);
                }
            }
        }
    }

    HOSTDEVICE void findInsideNodes(GridImp* grid, char innerType)
    {
        bool foundInsideNode = true;
        while (foundInsideNode)
        {
            foundInsideNode = false;
            for (uint index = 0; index < grid->getSize(); index++)
                this->setInsideNode(grid, index, foundInsideNode, innerType);
        }
    }

    HOSTDEVICE void setInsideNode(GridImp* grid, const uint &index, bool &insideNodeFound, char innerType)
    {
        if (grid->isNode(index, NEGATIVE_DIRECTION_BORDER))
            return;

        if (!grid->isNode(index, innerType) && grid->nodeInNextCellIs(index, innerType))
        {
            grid->setNodeTo(index, innerType);
            insideNodeFound = true;
        }
    }

    HOSTDEVICE void setNegativeDirBorderTo(GridImp* grid, const uint &index, char innerType)
    {
        if (grid->isNode(index, NEGATIVE_DIRECTION_BORDER))
            grid->setNodeTo(index, innerType);
    }

};



#endif


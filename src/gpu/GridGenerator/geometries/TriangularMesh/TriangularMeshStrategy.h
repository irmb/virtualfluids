#ifndef TriangularMeshStrategy_H
#define TriangularMeshStrategy_H

#include "global.h"
#include "GridGenerator_export.h"

class GridImp;
class TriangularMesh;
struct Triangle;

class GRIDGENERATOR_EXPORT TriangularMeshDiscretizationStrategy
{
public:
    TriangularMeshDiscretizationStrategy() {}
    virtual ~TriangularMeshDiscretizationStrategy() {}


    void discretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType)
    {  
        this->doDiscretize(triangularMesh, grid, InnerType, OuterType);
    }

private:
    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType) = 0;
    void removeOddBoundaryCellNodes(GridImp* grid);
};



class GRIDGENERATOR_EXPORT PointInObjectDiscretizationStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    PointInObjectDiscretizationStrategy() {}
    virtual ~PointInObjectDiscretizationStrategy() {}

    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType);
};

class GRIDGENERATOR_EXPORT RayCastingDiscretizationStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    RayCastingDiscretizationStrategy() {}
    virtual ~RayCastingDiscretizationStrategy() {}

    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType);
};

class GRIDGENERATOR_EXPORT PointUnderTriangleStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    PointUnderTriangleStrategy() {}
    virtual ~PointUnderTriangleStrategy() {}
    void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char innerType, char outerType) override;

private:
    HOSTDEVICE void meshReverse(Triangle& triangle, GridImp* grid, char innerType);

    HOSTDEVICE void findInsideNodes(GridImp* grid, char innerType);

    HOSTDEVICE void setInsideNode(GridImp* grid, const uint &index, bool &insideNodeFound, char innerType);

    HOSTDEVICE void setNegativeDirBorderTo(GridImp* grid, const uint &index, char innerType);

};



#endif


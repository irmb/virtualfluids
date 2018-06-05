#ifndef TriangularMeshStrategy_H
#define TriangularMeshStrategy_H

#include "GridGenerator/global.h"

#include "VirtualFluidsDefinitions.h"

class GridImp;
class TriangularMesh;
struct Triangle;

class VF_PUBLIC TriangularMeshDiscretizationStrategy
{
public:
    TriangularMeshDiscretizationStrategy() {}
    virtual ~TriangularMeshDiscretizationStrategy() {}


    void discretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType)
    {  
        this->doDiscretize(triangularMesh, grid, InnerType, OuterType);
        //this->removeOddBoundaryCellNodes(grid);
    }

private:
    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType) = 0;
    void removeOddBoundaryCellNodes(GridImp* grid);
};



class VF_PUBLIC PointInObjectDiscretizationStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    PointInObjectDiscretizationStrategy() {}
    virtual ~PointInObjectDiscretizationStrategy() {}

    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType);
};

class VF_PUBLIC RayCastingDiscretizationStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    RayCastingDiscretizationStrategy() {}
    virtual ~RayCastingDiscretizationStrategy() {}

    virtual void RayCastingDiscretizationStrategy::doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType);
};

class VF_PUBLIC PointUnderTriangleStrategy : public TriangularMeshDiscretizationStrategy
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


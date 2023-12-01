//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file TriangularMeshStrategy.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef TriangularMeshStrategy_H
#define TriangularMeshStrategy_H

#include "global.h"

class GridImp;
class TriangularMesh;
struct Triangle;

class TriangularMeshDiscretizationStrategy
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



class PointInObjectDiscretizationStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    PointInObjectDiscretizationStrategy() {}
    virtual ~PointInObjectDiscretizationStrategy() {}

    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType);
};

class RayCastingDiscretizationStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    RayCastingDiscretizationStrategy() {}
    virtual ~RayCastingDiscretizationStrategy() {}

    virtual void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType);
};

class PointUnderTriangleStrategy : public TriangularMeshDiscretizationStrategy
{
public:
    PointUnderTriangleStrategy() {}
    virtual ~PointUnderTriangleStrategy() {}
    void doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char innerType, char outerType) override;

private:
    void meshReverse(Triangle& triangle, GridImp* grid, char innerType);

    void findInsideNodes(GridImp* grid, char innerType);

    void setInsideNode(GridImp* grid, const uint &index, bool &insideNodeFound, char innerType);

    void setNegativeDirBorderTo(GridImp* grid, const uint &index, char innerType);

};



#endif


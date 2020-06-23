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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file Grid.h
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef GRID_H
#define GRID_H

#include "Core/LbmOrGks.h"

#include "global.h"

#include "geometries/Vertex/Vertex.h"

#include "grid/Cell.h"

struct Vertex;
struct Triangle;
class GridStrategy;
class GridInterface;
class Object;
class BoundingBox;

class VF_PUBLIC Grid
{
public:
    virtual ~Grid() {}

    virtual const Object* getObject() const = 0;

    virtual real getDelta() const = 0;
    virtual uint getSparseSize() const = 0;
    virtual uint getSize() const = 0;

    virtual real getStartX() const = 0;
    virtual real getStartY() const = 0;
    virtual real getStartZ() const = 0;

    virtual real getEndX() const = 0;
    virtual real getEndY() const = 0;
    virtual real getEndZ() const = 0;

    virtual uint getNumberOfNodesX() const = 0;
    virtual uint getNumberOfNodesY() const = 0;
    virtual uint getNumberOfNodesZ() const = 0;

    virtual int getSparseIndex(uint matrixIndex) const = 0;
    virtual char getFieldEntry(uint matrixIndex) const = 0;
    virtual void setFieldEntry(uint matrixIndex, char type) = 0;

    virtual int *getNeighborsX() const = 0;
    virtual int *getNeighborsY() const = 0;
    virtual int *getNeighborsZ() const = 0;
    virtual int *getNeighborsNegative() const = 0;

    virtual real* getDistribution() const = 0;
    virtual int* getDirection() const = 0;
    virtual int getStartDirection() const = 0;
    virtual int getEndDirection() const = 0;

    virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, uint *geo) const = 0;

    virtual SPtr<GridStrategy> getGridStrategy() const = 0;
    virtual void transIndexToCoords(uint index, real &x, real &y, real &z) const = 0;
    virtual uint transCoordToIndex(const real &x, const real &y, const real &z) const = 0;

    virtual void inital(const SPtr<Grid> fineGrid, uint numberOfLayers) = 0;
    
    virtual void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) = 0;
    virtual void setPeriodicityX(bool periodicity) = 0;
    virtual void setPeriodicityY(bool periodicity) = 0;
    virtual void setPeriodicityZ(bool periodicity) = 0;

    virtual bool getPeriodicityX() = 0;
    virtual bool getPeriodicityY() = 0;
    virtual bool getPeriodicityZ() = 0;

    virtual void freeMemory() = 0;

    virtual bool nodeInCellIs(Cell& cell, char type) const = 0;

    virtual void findSparseIndices(SPtr<Grid> fineGrid) = 0;

    virtual real getFirstFluidNode(real coords[3], int direction, real startCoord) const = 0;
    virtual real getLastFluidNode(real coords[3], int direction, real startCoord) const = 0;
};

#endif

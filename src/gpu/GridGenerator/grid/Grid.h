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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz, Martin Sch�nherr
//=======================================================================================
#ifndef GRID_H
#define GRID_H

#include "gpu/GridGenerator/global.h"

#include "gpu/GridGenerator/geometries/Vertex/Vertex.h"

#include "gpu/GridGenerator/grid/Cell.h"

class TriangularMesh;
struct Vertex;
struct Triangle;
class GridInterface;
class Object;
class BoundingBox;
enum class SideType;

class Grid
{
public:
    virtual ~Grid() = default;

    virtual SPtr<const Object> getObject() const = 0;

    virtual real getDelta() const = 0;
    virtual uint getSparseSize() const = 0;
    virtual uint getSize() const = 0;

    virtual real getStartX() const = 0;
    virtual real getStartY() const = 0;
    virtual real getStartZ() const = 0;

    virtual real getEndX() const = 0;
    virtual real getEndY() const = 0;
    virtual real getEndZ() const = 0;

    virtual Vertex getMinimumOnNode(Vertex exact) const = 0;
    virtual Vertex getMaximumOnNode(Vertex exact) const = 0;

    virtual uint getNumberOfNodesX() const = 0;
    virtual uint getNumberOfNodesY() const = 0;
    virtual uint getNumberOfNodesZ() const = 0;

    virtual uint getNumberOfNodesCF() const = 0;
    virtual uint getNumberOfNodesFC() const = 0;

    virtual int getSparseIndex(uint matrixIndex) const      = 0;
    virtual char getFieldEntry(uint matrixIndex) const = 0;
    virtual void setFieldEntry(uint matrixIndex, char type) = 0;

    virtual void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const = 0;
    virtual bool isSparseIndexInFluidNodeIndicesBorder(uint &sparseIndex) const = 0;

    virtual bool isStopperForBC(uint index) const = 0;

    virtual int *getNeighborsX() const = 0;
    virtual int *getNeighborsY() const = 0;
    virtual int *getNeighborsZ() const = 0;
    virtual int *getNeighborsNegative() const = 0;

    virtual uint *getCF_coarse() const = 0;
    virtual uint *getCF_fine() const   = 0;
    virtual uint *getCF_offset() const = 0;

    virtual uint *getFC_coarse() const = 0;
    virtual uint *getFC_fine() const   = 0;
    virtual uint *getFC_offset() const = 0;

    virtual real *getDistribution() const = 0;
    virtual const std::vector<int> &getDirection() const = 0;
    virtual int getStartDirection() const = 0;
    virtual int getEndDirection() const = 0;

    virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, uint *geo) const = 0;

    virtual void transIndexToCoords(uint index, real &x, real &y, real &z) const = 0;
    virtual uint transCoordToIndex(const real &x, const real &y, const real &z) const = 0;

    virtual void inital(const SPtr<Grid> fineGrid, uint numberOfLayers) = 0;
    
    virtual void setOddStart(bool xOddStart, bool yOddStart, bool zOddStart) = 0;

    virtual void findGridInterface(SPtr<Grid> grid) = 0;

    virtual void repairGridInterfaceOnMultiGPU(SPtr<Grid> fineGrid) = 0;

    virtual void limitToSubDomain(SPtr<BoundingBox> subDomainBox) = 0;

    virtual void enableFindSolidBoundaryNodes() = 0;
    virtual void enableComputeQs()              = 0;

    virtual void mesh(TriangularMesh &geometry) = 0;
    virtual void mesh(Object *object)           = 0;

    virtual void closeNeedleCells()         = 0;
    virtual void closeNeedleCellsThinWall() = 0;

    virtual void findQs(Object *object) = 0;

    virtual void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) = 0;
    virtual void setPeriodicityX(bool periodicity) = 0;
    virtual void setPeriodicityY(bool periodicity) = 0;
    virtual void setPeriodicityZ(bool periodicity) = 0;

    virtual bool getPeriodicityX() const = 0;
    virtual bool getPeriodicityY() const = 0;
    virtual bool getPeriodicityZ() const = 0;

    virtual void setPeriodicBoundaryShiftsOnXinY(real shift) = 0;
    virtual void setPeriodicBoundaryShiftsOnXinZ(real shift) = 0;
    virtual void setPeriodicBoundaryShiftsOnYinX(real shift) = 0;
    virtual void setPeriodicBoundaryShiftsOnYinZ(real shift) = 0;
    virtual void setPeriodicBoundaryShiftsOnZinX(real shift) = 0;
    virtual void setPeriodicBoundaryShiftsOnZinY(real shift) = 0;

    virtual void setEnableFixRefinementIntoTheWall(bool enableFixRefinementIntoTheWall) = 0;

    virtual void freeMemory() = 0;

    virtual bool nodeInCellIs(Cell& cell, char type) const = 0;

    virtual void findSparseIndices(SPtr<Grid> fineGrid) = 0;

    virtual real getFirstFluidNode(real coords[3], int direction, real startCoord) const = 0;
    virtual real getLastFluidNode(real coords[3], int direction, real startCoord) const = 0;

    virtual uint getNumberOfSolidBoundaryNodes() const                          = 0;
    virtual void setNumberOfSolidBoundaryNodes(uint numberOfSolidBoundaryNodes) = 0;

    virtual real getQValue(const uint index, const uint dir) const = 0;
    virtual uint getQPatch(const uint index) const                 = 0;

    virtual void setInnerRegionFromFinerGrid(bool innerRegionFromFinerGrid) = 0;

    virtual void setNumberOfLayers(uint numberOfLayers) = 0;

    virtual void findCommunicationIndices(int direction, SPtr<BoundingBox> subDomainBox, bool doShift) = 0;

    virtual uint getNumberOfSendNodes(int direction)    = 0;
    virtual uint getNumberOfReceiveNodes(int direction) = 0;

    virtual bool isSendNode(int index) const                = 0;
    virtual bool isReceiveNode(int index) const             = 0;
    virtual uint getSendIndex(int direction, uint index)    = 0;
    virtual uint getReceiveIndex(int direction, uint index) = 0;

    virtual void repairCommunicationIndices(int direction) = 0;

    virtual bool nodeHasBC(uint index) const = 0;

    virtual std::vector<SideType> getBCAlreadySet() = 0;
    virtual void addBCalreadySet(SideType side) = 0;

    // needed for CUDA Streams 
    virtual void findFluidNodeIndices(bool onlyBulk) = 0;
    virtual uint getNumberOfFluidNodes() const = 0;
    virtual void getFluidNodeIndices(uint *fluidNodeIndices) const = 0;

    virtual void findFluidNodeIndicesBorder() = 0;
    virtual uint getNumberOfFluidNodesBorder() const = 0;
    virtual void getFluidNodeIndicesBorder(uint *fluidNodeIndicesBorder) const = 0;

    virtual void addFluidNodeIndicesMacroVars(std::vector<uint> _fluidNodeIndicesMacroVars) = 0;
    virtual void addFluidNodeIndicesApplyBodyForce(std::vector<uint> _fluidNodeIndicesApplyBodyForce) = 0;
    virtual void addFluidNodeIndicesAllFeatures(std::vector<uint> _fluidNodeIndicesAllFeatures) = 0;
    virtual void sortFluidNodeIndicesMacroVars() = 0;
    virtual void sortFluidNodeIndicesApplyBodyForce() = 0;
    virtual void sortFluidNodeIndicesAllFeatures() = 0;

    virtual uint getNumberOfFluidNodeIndicesMacroVars() const = 0;
    virtual uint getNumberOfFluidNodeIndicesApplyBodyForce() const = 0;
    virtual uint getNumberOfFluidNodeIndicesAllFeatures() const = 0; 
    virtual void getFluidNodeIndicesMacroVars(uint *fluidNodeIndicesMacroVars) const = 0;
    virtual void getFluidNodeIndicesApplyBodyForce(uint *fluidNodeIndicesApplyBodyForce) const = 0;
    virtual void getFluidNodeIndicesAllFeatures(uint *fluidNodeIndicesAllFeatures) const = 0;
};

#endif

//! \}

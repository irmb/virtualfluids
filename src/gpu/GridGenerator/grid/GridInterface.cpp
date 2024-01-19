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
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "GridInterface.h"

#include <iostream>
#include <cstring>

#include "grid/distributions/D3Q27.h"
#include "grid/GridImp.h"
#include "grid/Field.h"
#include "grid/NodeValues.h"

#include "lbm/constants/D3Q27.h"

using namespace vf::lbm;
using namespace vf::gpu;

GridInterface::GridInterface()
{

}

GridInterface::~GridInterface()
{

}


void GridInterface::findInterfaceCF(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->getField().isFluid(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid)
        return;

    const uint indexOnFineGridCF = getCoarseToFineIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridCF == INVALID_INDEX)
        return;

    const bool fineGridNodeIsFluid = fineGrid->getField().isFluid(indexOnFineGridCF);
    if (!fineGridNodeIsFluid)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    for(const auto dir : coarseGrid->distribution)
    {
        const bool isFineGridNeighborInvalid = isNeighborFineInvalid(x + dir[0] * coarseGrid->getDelta(), y + dir[1] * coarseGrid->getDelta(), z + dir[2] * coarseGrid->getDelta(), coarseGrid, fineGrid);
        if(isFineGridNeighborInvalid)
        {
            cf.coarse[cf.numberOfEntries] = this->findOffsetCF(indexOnCoarseGrid, coarseGrid, cf.numberOfEntries);
            cf.fine[cf.numberOfEntries]   = indexOnFineGridCF;

            cf.numberOfEntries++;

            coarseGrid->setNonStopperOutOfGridCellTo(indexOnCoarseGrid, FLUID_CFC);
            fineGrid->setNonStopperOutOfGridCellTo(indexOnFineGridCF, FLUID_CFF);
            break;
        }
    }
}


void GridInterface::findBoundaryGridInterfaceCF(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsBoundaryStopper = coarseGrid->getField().is(indexOnCoarseGrid, STOPPER_OUT_OF_GRID_BOUNDARY);
    if (!nodeOnCoarseGridIsBoundaryStopper)
        return;

    const uint indexOnFineGridCF = getCoarseToFineIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridCF == INVALID_INDEX)
        return;

    const bool fineGridNodeIsBoundaryStopper = fineGrid->getField().is(indexOnFineGridCF, STOPPER_OUT_OF_GRID_BOUNDARY);
    if (!fineGridNodeIsBoundaryStopper)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    for(const auto dir : coarseGrid->distribution)
    {
        const bool isFineGridNeighborInvalid = isNeighborFineInvalid(x + dir[0] * coarseGrid->getDelta(), y + dir[1] * coarseGrid->getDelta(), z + dir[2] * coarseGrid->getDelta(), coarseGrid, fineGrid);
        if(isFineGridNeighborInvalid)
        {
            cf.coarse[cf.numberOfEntries] = this->findOffsetCF(indexOnCoarseGrid, coarseGrid, cf.numberOfEntries);
            cf.fine[cf.numberOfEntries]   = indexOnFineGridCF;

            cf.numberOfEntries++;

            coarseGrid->setNonStopperOutOfGridCellTo(indexOnCoarseGrid, FLUID_CFC);
            fineGrid->setNonStopperOutOfGridCellTo(indexOnFineGridCF, FLUID_CFF);
            break;
        }
    }
}


void GridInterface::findInterfaceFC(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->getField().isFluid(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsCoarseToFine = coarseGrid->getField().isCoarseToFineNode(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid || nodeOnCoarseGridIsCoarseToFine)
        return;

    const uint indexOnFineGridFC = getFineToCoarseIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridFC == INVALID_INDEX)
        return;

    const bool fineGridNodeIsFluid = fineGrid->getField().isFluid(indexOnFineGridFC);
    if (!fineGridNodeIsFluid)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    for (const auto dir : coarseGrid->distribution)
    {
        const uint neighborIndex = coarseGrid->transCoordToIndex(x + dir[0] * coarseGrid->getDelta(), y + dir[1] * coarseGrid->getDelta(), z + dir[2] * coarseGrid->getDelta());
        if (neighborIndex != INVALID_INDEX)
        {
            const bool neighborBelongsToCoarseToFineInterpolationCell = coarseGrid->getField().isCoarseToFineNode(neighborIndex);
            if (neighborBelongsToCoarseToFineInterpolationCell)
            {
                fc.coarse[fc.numberOfEntries] = indexOnCoarseGrid;
                fc.fine[fc.numberOfEntries] = this->findOffsetFC(indexOnFineGridFC, fineGrid, fc.numberOfEntries);

                fc.numberOfEntries++;

                fineGrid->setNonStopperOutOfGridCellTo(indexOnFineGridFC, FLUID_FCF);
                coarseGrid->getField().setFieldEntry(indexOnCoarseGrid, FLUID_FCC);
                break;
            }
        }
    }
}

void GridInterface::findOverlapStopper(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->getField().isFluid(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsCoarseToFine = coarseGrid->getField().isCoarseToFineNode(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsFineToCoarse = coarseGrid->getField().isFineToCoarseNode(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid || nodeOnCoarseGridIsCoarseToFine || nodeOnCoarseGridIsFineToCoarse)
        return;

    const int indexOnFineGridFC = getFineToCoarseIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridFC == -1)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    bool neighborBelongsToFineToCoarseInterpolationCell = false;
    for (const auto dir : coarseGrid->distribution)
    {
        //if (dir[0] > 0 || dir[1] > 0 || dir[2] > 0)  //only Esoteric Twist stopper, not perfectly implemented
        //    continue;                                   //should not be here, should be made conditional

        const uint neighborIndex = coarseGrid->transCoordToIndex(x + dir[0] * coarseGrid->getDelta(), y + dir[1] * coarseGrid->getDelta(), z + dir[2] * coarseGrid->getDelta());
        neighborBelongsToFineToCoarseInterpolationCell = neighborIndex != INVALID_INDEX ? coarseGrid->getField().isFineToCoarseNode(neighborIndex) : false;
        if (neighborBelongsToFineToCoarseInterpolationCell)
        {
            coarseGrid->getField().setFieldEntryToStopperCoarseUnderFine(indexOnCoarseGrid);
            break;
        }
    }

    //should be inside of fine grid and can be deleted
    if(!neighborBelongsToFineToCoarseInterpolationCell && (fineGrid->getField().isInvalidSolid(indexOnFineGridFC) || 
                                                           fineGrid->getField().isFluid(indexOnFineGridFC) ||
                                                           fineGrid->getField().is(indexOnFineGridFC, STOPPER_SOLID) || 
                                                           fineGrid->getField().is(indexOnFineGridFC, BC_SOLID)))
        coarseGrid->getField().setFieldEntryToInvalidCoarseUnderFine(indexOnCoarseGrid);
}

 void GridInterface::findInvalidBoundaryNodes(const uint& indexOnCoarseGrid, GridImp* coarseGrid)
{
     if( !coarseGrid->getField().is(indexOnCoarseGrid, STOPPER_OUT_OF_GRID_BOUNDARY ) ) return;

     if( !coarseGrid->hasNeighborOfType(indexOnCoarseGrid, FLUID) &&
         !coarseGrid->hasNeighborOfType(indexOnCoarseGrid, FLUID_CFC) &&
         !coarseGrid->hasNeighborOfType(indexOnCoarseGrid, FLUID_CFF) )
         coarseGrid->getField().setFieldEntryToInvalidCoarseUnderFine(indexOnCoarseGrid);
}

bool GridInterface::isNeighborFineInvalid(real x, real y, real z, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    const uint neighbor = coarseGrid->transCoordToIndex(x, y, z);

    if( neighbor == INVALID_INDEX )
        return false;

    if( (neighbor != INVALID_INDEX) && (coarseGrid->getField().isStopperOutOfGrid(neighbor) || coarseGrid->getField().is(neighbor, STOPPER_OUT_OF_GRID_BOUNDARY)) )
        return false;

    const uint indexOnFineGrid = getCoarseToFineIndexOnFineGrid(neighbor, coarseGrid, fineGrid);
    if (indexOnFineGrid == INVALID_INDEX)
        return true;

    return fineGrid->getField().isInvalidOutOfGrid(indexOnFineGrid) || fineGrid->getField().isStopperOutOfGrid(indexOnFineGrid);
}

uint GridInterface::getCoarseToFineIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    if( indexOnCoarseGrid == INVALID_INDEX )
        return INVALID_INDEX;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);
    const real xFine = x + (fineGrid->getDelta() * 0.5);
    const real yFine = y + (fineGrid->getDelta() * 0.5);
    const real zFine = z + (fineGrid->getDelta() * 0.5);

    return fineGrid->transCoordToIndex(xFine, yFine, zFine);
}

uint GridInterface::getFineToCoarseIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);
    const real xFine = x - (fineGrid->getDelta() * 0.5);
    const real yFine = y - (fineGrid->getDelta() * 0.5);
    const real zFine = z - (fineGrid->getDelta() * 0.5);

    return fineGrid->transCoordToIndex(xFine, yFine, zFine);
}

void GridInterface::findForGridInterfaceSparseIndexCF(GridImp* coarseGrid, GridImp* fineGrid, uint index)
{
    findSparseIndex(cf.coarse, coarseGrid, index);
    findSparseIndex(cf.fine, fineGrid, index);

    if( cf.coarse[index] == 686916 )
        printf("%d ===> %d \n", cf.coarse[index], cf.fine[index]);
}

void GridInterface::findForGridInterfaceSparseIndexFC(GridImp* coarseGrid, GridImp* fineGrid, uint index)
{
    findSparseIndex(fc.coarse, coarseGrid, index);
    findSparseIndex(fc.fine, fineGrid, index);
}

void GridInterface::repairGridInterfaceOnMultiGPU(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid)
{
    {
        std::vector<uint> tmpCFC;
        std::vector<uint> tmpCFF;
        std::vector<uint> tmpCFOffset;

        for (uint index = 0; index < cf.numberOfEntries; index++) {

            real x, y, z;
            coarseGrid->transIndexToCoords(this->cf.coarse[index], x, y, z);
            Cell cell(x, y, z, coarseGrid->getDelta());

            if (coarseGrid->cellContainsOnly(cell, FLUID_CFC)) {
                tmpCFC.push_back     (this->cf.coarse[index]);
                tmpCFF.push_back     (this->cf.fine[index]);
                tmpCFOffset.push_back(this->cf.offset[index]);
            }
        }

        delete[] cf.coarse;
        delete[] cf.fine;
        delete[] cf.offset;

        cf.numberOfEntries = (uint)tmpCFC.size();

        cf.coarse = new uint[cf.numberOfEntries];
        cf.fine   = new uint[cf.numberOfEntries];
        cf.offset = new uint[cf.numberOfEntries];

        memcpy(cf.coarse, tmpCFC.data()     , sizeof(uint)*cf.numberOfEntries);
        memcpy(cf.fine  , tmpCFF.data()     , sizeof(uint)*cf.numberOfEntries);
        memcpy(cf.offset, tmpCFOffset.data(), sizeof(uint)*cf.numberOfEntries);
    }

    {
        std::vector<uint> tmpFCF;
        std::vector<uint> tmpFCC;
        std::vector<uint> tmpFCOffset;

        for (uint index = 0; index < fc.numberOfEntries; index++) {

            real x, y, z;
            fineGrid->transIndexToCoords(this->fc.fine[index], x, y, z);
            Cell cell(x, y, z, fineGrid->getDelta());

            if (fineGrid->cellContainsOnly(cell, FLUID_FCF)) {
                tmpFCF.push_back     (this->fc.fine[index]);
                tmpFCC.push_back     (this->fc.coarse[index]);
                tmpFCOffset.push_back(this->fc.offset[index]);
            }
        }
        
        delete[] fc.fine;
        delete[] fc.coarse;
        delete[] fc.offset;

        fc.numberOfEntries = (uint)tmpFCC.size();
        
        fc.fine   = new uint[fc.numberOfEntries];
        fc.coarse = new uint[fc.numberOfEntries];
        fc.offset = new uint[fc.numberOfEntries];
        
        memcpy(fc.fine  , tmpFCF.data()     , sizeof(uint)*fc.numberOfEntries);
        memcpy(fc.coarse, tmpFCC.data()     , sizeof(uint)*fc.numberOfEntries);
        memcpy(fc.offset, tmpFCOffset.data(), sizeof(uint)*fc.numberOfEntries);
    }
}

void GridInterface::findSparseIndex(uint* indices, GridImp* grid, uint index)
{
    const uint matrixIndex = indices[index];
    const uint sparseIndex = grid->getSparseIndex(matrixIndex);
    indices[index] = sparseIndex;
}

uint GridInterface::findOffsetCF(const uint& indexOnCoarseGrid, GridImp* coarseGrid, uint interfaceIndex)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    Cell cell(x, y, z, coarseGrid->getDelta());

    if( coarseGrid->cellContainsOnly( cell, FLUID, FLUID_CFC ) ){
        this->cf.offset[ interfaceIndex ] = dir::d000;
        return indexOnCoarseGrid;
    }

    uint dirIndex = 0;
    for(const auto dir : coarseGrid->distribution){
    
        Cell neighborCell( x + dir[0] * coarseGrid->getDelta(), 
                           y + dir[1] * coarseGrid->getDelta(), 
                           z + dir[2] * coarseGrid->getDelta(), 
                           coarseGrid->getDelta() );

        if( coarseGrid->cellContainsOnly( neighborCell, FLUID, FLUID_CFC ) ){
            this->cf.offset[ interfaceIndex ] = dirIndex;

            return coarseGrid->transCoordToIndex( x + dir[0] * coarseGrid->getDelta(),
                                                  y + dir[1] * coarseGrid->getDelta(),
                                                  z + dir[2] * coarseGrid->getDelta() );
        }
    
        dirIndex++;
    }

    // this point should never be reached
    return indexOnCoarseGrid;
}

uint GridInterface::findOffsetFC(const uint& indexOnFineGrid, GridImp* fineGrid, uint interfaceIndex)
{
    real x, y, z;
    fineGrid->transIndexToCoords(indexOnFineGrid, x, y, z);

    Cell cell(x, y, z, fineGrid->getDelta());

    if( fineGrid->cellContainsOnly( cell, FLUID, FLUID_FCF ) ){
        this->fc.offset[ interfaceIndex ] = dir::d000;
        return indexOnFineGrid;
    }

    uint dirIndex = 0;
    for(const auto dir : fineGrid->distribution){
    
        Cell neighborCell( x + dir[0] * fineGrid->getDelta(), 
                           y + dir[1] * fineGrid->getDelta(), 
                           z + dir[2] * fineGrid->getDelta(), 
                           fineGrid->getDelta() );

        if( fineGrid->cellContainsOnly( neighborCell, FLUID, FLUID_CFC ) ){
            this->fc.offset[interfaceIndex] = dirIndex;

            return fineGrid->transCoordToIndex(x + dir[0] * fineGrid->getDelta(),
                                               y + dir[1] * fineGrid->getDelta(),
                                               z + dir[2] * fineGrid->getDelta());
        }
    
        dirIndex++;
    }

    // this point should never be reached
    return indexOnFineGrid;
}

void GridInterface::print() const
{
    printf("Grid Interface - CF nodes: %d, FC nodes: %d\n", cf.numberOfEntries, fc.numberOfEntries);
}

//! \}

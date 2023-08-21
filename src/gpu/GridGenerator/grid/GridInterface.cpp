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
//! \file GridInterface.cpp
//! \ingroup grid
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

void GridInterface::findInterfaceBaseToNested(const uint& indexOnBaseGrid, GridImp* baseGrid, GridImp* nestedGrid, bool isRotatingGrid)
{
    const bool nodeOnBaseGridIsFluid = baseGrid->getField().isFluid(indexOnBaseGrid);
    if (!nodeOnBaseGridIsFluid)
        return;

    const uint indexOnNestedGridBN = getBaseToNestedIndexOnFineGrid(indexOnBaseGrid, baseGrid, nestedGrid);
    if (indexOnNestedGridBN == INVALID_INDEX)
        return;

    const bool nestedGridNodeIsFluid = nestedGrid->getField().isFluid(indexOnNestedGridBN);
    if (!nestedGridNodeIsFluid)
        return;

    real x, y, z;
    baseGrid->transIndexToCoords(indexOnBaseGrid, x, y, z);

    for(const auto dir : baseGrid->distribution)
    {
        const bool isNestedGridNeighborInvalid = isNeighborNestedInvalid(x + dir[0] * baseGrid->getDelta(), y + dir[1] * baseGrid->getDelta(), z + dir[2] * baseGrid->getDelta(), baseGrid, nestedGrid);
        if(isNestedGridNeighborInvalid)
        {
            bn.base[bn.numberOfEntries] = this->findOffsetBaseToNested(indexOnBaseGrid, baseGrid, bn.numberOfEntries);
            bn.nested[bn.numberOfEntries]   = indexOnNestedGridBN;

            bn.numberOfEntries++;

            baseGrid->setNonStopperOutOfGridCellTo(indexOnBaseGrid, FLUID_BNB);

            if (isRotatingGrid) {
                nestedGrid->getField().setFieldEntry(indexOnNestedGridBN, FLUID_BNN);
            } else {
                nestedGrid->setNonStopperOutOfGridCellTo(indexOnNestedGridBN, FLUID_CFF);
            }

            break;
        }
    }
}


void GridInterface::findBoundaryGridInterfaceCF(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsBoundaryStopper = coarseGrid->getField().is(indexOnCoarseGrid, STOPPER_OUT_OF_GRID_BOUNDARY);
    if (!nodeOnCoarseGridIsBoundaryStopper)
        return;

    const uint indexOnFineGridCF = getBaseToNestedIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridCF == INVALID_INDEX)
        return;

    const bool fineGridNodeIsBoundaryStopper = fineGrid->getField().is(indexOnFineGridCF, STOPPER_OUT_OF_GRID_BOUNDARY);
    if (!fineGridNodeIsBoundaryStopper)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    for(const auto dir : coarseGrid->distribution)
    {
        const bool isFineGridNeighborInvalid = isNeighborNestedInvalid(x + dir[0] * coarseGrid->getDelta(), y + dir[1] * coarseGrid->getDelta(), z + dir[2] * coarseGrid->getDelta(), coarseGrid, fineGrid);
        if(isFineGridNeighborInvalid)
        {
			bn.base[bn.numberOfEntries] = this->findOffsetBaseToNested(indexOnCoarseGrid, coarseGrid, bn.numberOfEntries);
			bn.nested[bn.numberOfEntries]   = indexOnFineGridCF;

            bn.numberOfEntries++;

            coarseGrid->setNonStopperOutOfGridCellTo(indexOnCoarseGrid, FLUID_CFC);
            fineGrid->setNonStopperOutOfGridCellTo(indexOnFineGridCF, FLUID_CFF);

            break;
        }
    }
}

uint GridInterface::findNeighborIndex(real coordX, real coordY, real coordZ, GridImp *grid, Direction dir)
{
    return grid->transCoordToIndex(coordX + dir[0] * grid->getDelta(), coordY + dir[1] * grid->getDelta(),
                                   coordZ + dir[2] * grid->getDelta());
}

void GridInterface::findInterfaceFC(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->getField().isFluid(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsCoarseToFine = coarseGrid->getField().isBaseToNestedNode(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid || nodeOnCoarseGridIsCoarseToFine)
        return;

    const uint indexOnFineGridFC = getNestedToBaseIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridFC == INVALID_INDEX)
        return;

    const bool fineGridNodeIsFluid = fineGrid->getField().isFluid(indexOnFineGridFC);
    if (!fineGridNodeIsFluid)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    for (const auto dir : coarseGrid->distribution)
    {
        const uint neighborIndex = findNeighborIndex(x, y, z, coarseGrid, dir);
		if (neighborIndex != INVALID_INDEX)
		{
			const bool neighborBelongsToCoarseToFineInterpolationCell = coarseGrid->getField().isBaseToNestedNode(neighborIndex);
			if (neighborBelongsToCoarseToFineInterpolationCell)
			{
				nb.base[nb.numberOfEntries] = indexOnCoarseGrid;
				nb.nested[nb.numberOfEntries] = this->findOffsetNestedToBase(indexOnFineGridFC, fineGrid, nb.numberOfEntries);

				nb.numberOfEntries++;

				fineGrid->setNonStopperOutOfGridCellTo(indexOnFineGridFC, FLUID_FCF);
				coarseGrid->getField().setFieldEntry(indexOnCoarseGrid, FLUID_FCC);
				break;
			}
		}
    }
}

void GridInterface::findInterfaceNestedToBaseWithGap(uint indexOnBaseGrid, GridImp *baseGrid, GridImp *nestedGrid)
{
    const bool nodeOnBaseGridIsFluid = baseGrid->getField().isFluid(indexOnBaseGrid);
    const bool nodeOnBaseGridIsBaseToNested = baseGrid->getField().isBaseToNestedNode(indexOnBaseGrid);
    if (!nodeOnBaseGridIsFluid || nodeOnBaseGridIsBaseToNested) return;

    const uint indexOnNestedForNestedToBase = getNestedToBaseIndexOnFineGrid(indexOnBaseGrid, baseGrid, nestedGrid);

    real x, y, z;
    baseGrid->transIndexToCoords(indexOnBaseGrid, x, y, z);

    for (const auto dir : baseGrid->distribution) {
        const uint neighborIndex = findNeighborIndex(x, y, z, baseGrid, dir);
        if (neighborIndex == INVALID_INDEX) return;
        const bool neighborIsInterpolationGapNode = baseGrid->getField().isInterpolationGapNode(neighborIndex);

        if (neighborIsInterpolationGapNode) {
            nb.base[nb.numberOfEntries] = indexOnBaseGrid;
            nb.nested[nb.numberOfEntries] =
                this->findOffsetNestedToBase(indexOnNestedForNestedToBase, nestedGrid, nb.numberOfEntries);
            nb.numberOfEntries++;

            baseGrid->getField().setFieldEntry(indexOnBaseGrid, FLUID_NBB);
            nestedGrid->setNonStopperOutOfGridCellTo(indexOnNestedForNestedToBase, FLUID_NBN);
            break;
        }
    }
}

void GridInterface::findInterpolationGapOnBaseGrid(const uint &indexOnBaseGrid, GridImp *baseGrid, GridImp *nestedGrid)
{
    const bool nodeIsFluid = baseGrid->getField().isFluid(indexOnBaseGrid);
    const bool nodeIsBaseToNested = baseGrid->getField().isBaseToNestedNode(indexOnBaseGrid);
    if (!nodeIsFluid || nodeIsBaseToNested) return;

    const uint indexOnFineGridNestedToBase = getNestedToBaseIndexOnFineGrid(indexOnBaseGrid, baseGrid, nestedGrid);
    if (indexOnFineGridNestedToBase == INVALID_INDEX) return;

    const bool nestedGridNodeIsFluid = nestedGrid->getField().isFluid(indexOnFineGridNestedToBase);
    if (!nestedGridNodeIsFluid) return;

    real x, y, z;
    baseGrid->transIndexToCoords(indexOnBaseGrid, x, y, z);

    for (const auto dir : baseGrid->distribution) {
        const uint neighborIndex = findNeighborIndex(x, y, z, baseGrid, dir);
        if (neighborIndex == INVALID_INDEX) continue;
        const bool neighborIsBaseToNested = baseGrid->getField().isBaseToNestedNode(neighborIndex);

        if (neighborIsBaseToNested) {
            baseGrid->getField().setFieldEntry(indexOnBaseGrid, INTERPOLATION_GAP);
            break;
        }
    }
}

void GridInterface::findSecondInterpolationGapOnBaseGrid(const uint &indexOnBaseGrid, GridImp *baseGrid, GridImp *nestedGrid)
{
    const bool nodeIsFluid = baseGrid->getField().isFluid(indexOnBaseGrid);
    const bool nodeIsBaseToNested = baseGrid->getField().isBaseToNestedNode(indexOnBaseGrid);
    const bool nodeIsGap = baseGrid->getField().is(indexOnBaseGrid, INTERPOLATION_GAP);
    if (!nodeIsFluid || nodeIsBaseToNested || nodeIsGap) return;

    const uint indexOnFineGridNestedToBase = getNestedToBaseIndexOnFineGrid(indexOnBaseGrid, baseGrid, nestedGrid);
    if (indexOnFineGridNestedToBase == INVALID_INDEX) return;

    const bool nestedGridNodeIsFluid = nestedGrid->getField().isFluid(indexOnFineGridNestedToBase);
    if (!nestedGridNodeIsFluid) return;

    real x, y, z;
    baseGrid->transIndexToCoords(indexOnBaseGrid, x, y, z);

    for (const auto dir : baseGrid->distribution) {
        const uint neighborIndex = findNeighborIndex(x, y, z, baseGrid, dir);
        if (neighborIndex == INVALID_INDEX) continue;
        const bool neighborIsGap = baseGrid->getField().is(neighborIndex, INTERPOLATION_GAP);

        if (neighborIsGap) {
            baseGrid->getField().setFieldEntry(indexOnBaseGrid, INTERPOLATION_GAP2);
            break;
        }
    }
}

void GridInterface::findInterpolationGapOnNestedGrid(const uint &indexOnNestedGrid, GridImp *nestedGrid)
{
    const bool nodeIsFluid = nestedGrid->getField().isFluid(indexOnNestedGrid);
    const bool nodeIsBaseToNested = nestedGrid->getField().is(indexOnNestedGrid, FLUID_BNN);
    if (!nodeIsFluid || nodeIsBaseToNested) return;

    real x, y, z;
    nestedGrid->transIndexToCoords(indexOnNestedGrid, x, y, z);

    for (const auto dir : nestedGrid->distribution) {
        const uint neighborIndex = findNeighborIndex(x, y, z, nestedGrid, dir);
        if (neighborIndex == INVALID_INDEX) continue;
        const bool neighborIsToBaseToNested = nestedGrid->getField().is(neighborIndex, FLUID_BNN);

        if (neighborIsToBaseToNested) {
            nestedGrid->getField().setFieldEntry(indexOnNestedGrid, INTERPOLATION_GAP);
            break;
        }
    }
}

void GridInterface::findOverlapStopper(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->getField().isFluid(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsCoarseToFine = coarseGrid->getField().isBaseToNestedNode(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsFineToCoarse = coarseGrid->getField().isFineToCoarseNode(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid || nodeOnCoarseGridIsCoarseToFine || nodeOnCoarseGridIsFineToCoarse)
        return;

    const int indexOnFineGridFC = getNestedToBaseIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridFC == -1)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    bool neighborBelongsToFineToCoarseInterpolationCell = false;
    for (const auto dir : coarseGrid->distribution)
    {
        //if (dir[0] > 0 || dir[1] > 0 || dir[2] > 0)  //only Esoteric Twist stopper, not perfectly implemented
        //    continue;								   //should not be here, should be made conditional

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

bool GridInterface::isNeighborNestedInvalid(real x, real y, real z, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    const uint neighbor = coarseGrid->transCoordToIndex(x, y, z);

    if( neighbor == INVALID_INDEX )
        return false;

    if( (neighbor != INVALID_INDEX) && (coarseGrid->getField().isStopperOutOfGrid(neighbor) || coarseGrid->getField().is(neighbor, STOPPER_OUT_OF_GRID_BOUNDARY)) )
        return false;

    const uint indexOnFineGrid = getBaseToNestedIndexOnFineGrid(neighbor, coarseGrid, fineGrid);
    if (indexOnFineGrid == INVALID_INDEX)
        return true;

    return fineGrid->getField().isInvalidOutOfGrid(indexOnFineGrid) || fineGrid->getField().isStopperOutOfGrid(indexOnFineGrid);
}

uint GridInterface::getBaseToNestedIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
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

uint GridInterface::getNestedToBaseIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
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
    findSparseIndex(bn.base, coarseGrid, index);
    findSparseIndex(bn.nested, fineGrid, index);

    if( bn.base[index] == 686916 )
        printf("%d ===> %d \n", bn.base[index], bn.nested[index]);
}

void GridInterface::findForGridInterfaceSparseIndexFC(GridImp* coarseGrid, GridImp* fineGrid, uint index)
{
    findSparseIndex(nb.base, coarseGrid, index);
    findSparseIndex(nb.nested, fineGrid, index);
}

void GRIDGENERATOR_EXPORT GridInterface::repairGridInterfaceOnMultiGPU(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid)
{
    {
        std::vector<uint> tmpCFC;
        std::vector<uint> tmpCFF;
        std::vector<uint> tmpCFOffset;

        for (uint index = 0; index < bn.numberOfEntries; index++) {

            real x, y, z;
            coarseGrid->transIndexToCoords(this->bn.base[index], x, y, z);
            Cell cell(x, y, z, coarseGrid->getDelta());

            if (coarseGrid->cellContainsOnly(cell, FLUID_CFC)) {
                tmpCFC.push_back     (this->bn.base[index]);
                tmpCFF.push_back     (this->bn.nested[index]);
                tmpCFOffset.push_back(this->bn.offset[index]);
            }
        }

        delete[] bn.base;
        delete[] bn.nested;
        delete[] bn.offset;

        bn.numberOfEntries = (uint)tmpCFC.size();

        bn.base = new uint[bn.numberOfEntries];
        bn.nested   = new uint[bn.numberOfEntries];
        bn.offset = new uint[bn.numberOfEntries];

        memcpy(bn.base, tmpCFC.data()     , sizeof(uint)*bn.numberOfEntries);
        memcpy(bn.nested  , tmpCFF.data()     , sizeof(uint)*bn.numberOfEntries);
        memcpy(bn.offset, tmpCFOffset.data(), sizeof(uint)*bn.numberOfEntries);
    }

    {
        std::vector<uint> tmpFCF;
        std::vector<uint> tmpFCC;
        std::vector<uint> tmpFCOffset;

        for (uint index = 0; index < nb.numberOfEntries; index++) {

            real x, y, z;
            fineGrid->transIndexToCoords(this->nb.nested[index], x, y, z);
            Cell cell(x, y, z, fineGrid->getDelta());

            if (fineGrid->cellContainsOnly(cell, FLUID_FCF)) {
                tmpFCF.push_back     (this->nb.nested[index]);
                tmpFCC.push_back     (this->nb.base[index]);
                tmpFCOffset.push_back(this->nb.offset[index]);
            }
        }
        
        delete[] nb.nested;
        delete[] nb.base;
        delete[] nb.offset;

        nb.numberOfEntries = (uint)tmpFCC.size();
        
        nb.nested   = new uint[nb.numberOfEntries];
        nb.base = new uint[nb.numberOfEntries];
        nb.offset = new uint[nb.numberOfEntries];
        
        memcpy(nb.nested  , tmpFCF.data()     , sizeof(uint)*nb.numberOfEntries);
        memcpy(nb.base, tmpFCC.data()     , sizeof(uint)*nb.numberOfEntries);
        memcpy(nb.offset, tmpFCOffset.data(), sizeof(uint)*nb.numberOfEntries);
    }
}

void GridInterface::findSparseIndex(uint* indices, GridImp* grid, uint index)
{
    const uint matrixIndex = indices[index];
    const uint sparseIndex = grid->getSparseIndex(matrixIndex);
    indices[index] = sparseIndex;
}

uint GridInterface::findOffsetBaseToNested(const uint& indexOnCoarseGrid, GridImp* coarseGrid, uint interfaceIndex)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    Cell cell(x, y, z, coarseGrid->getDelta());

    if( coarseGrid->cellContainsOnly( cell, FLUID, FLUID_CFC ) ){
        this->bn.offset[ interfaceIndex ] = dir::DIR_000;
        return indexOnCoarseGrid;
    }

    uint dirIndex = 0;
    for(const auto dir : coarseGrid->distribution){
    
        Cell neighborCell( x + dir[0] * coarseGrid->getDelta(), 
                           y + dir[1] * coarseGrid->getDelta(), 
                           z + dir[2] * coarseGrid->getDelta(), 
                           coarseGrid->getDelta() );

        if( coarseGrid->cellContainsOnly( neighborCell, FLUID, FLUID_CFC ) ){
            this->bn.offset[ interfaceIndex ] = dirIndex;

			return coarseGrid->transCoordToIndex( x + dir[0] * coarseGrid->getDelta(),
				                                  y + dir[1] * coarseGrid->getDelta(),
				                                  z + dir[2] * coarseGrid->getDelta() );
        }
    
        dirIndex++;
    }

	// this point should never be reached
	return indexOnCoarseGrid;
}

uint GridInterface::findOffsetNestedToBase(const uint& indexOnFineGrid, GridImp* fineGrid, uint interfaceIndex)
{
    real x, y, z;
    fineGrid->transIndexToCoords(indexOnFineGrid, x, y, z);

    Cell cell(x, y, z, fineGrid->getDelta());

    if( fineGrid->cellContainsOnly( cell, FLUID, FLUID_FCF ) ){
        this->nb.offset[ interfaceIndex ] = dir::DIR_000;
        return indexOnFineGrid;
    }

    uint dirIndex = 0;
    for(const auto dir : fineGrid->distribution){
    
        Cell neighborCell( x + dir[0] * fineGrid->getDelta(), 
                           y + dir[1] * fineGrid->getDelta(), 
                           z + dir[2] * fineGrid->getDelta(), 
                           fineGrid->getDelta() );

        if( fineGrid->cellContainsOnly( neighborCell, FLUID, FLUID_CFC ) ){
			this->nb.offset[interfaceIndex] = dirIndex;

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
    printf("Grid Interface - CF nodes: %d, FC nodes: %d\n", bn.numberOfEntries, nb.numberOfEntries);
}

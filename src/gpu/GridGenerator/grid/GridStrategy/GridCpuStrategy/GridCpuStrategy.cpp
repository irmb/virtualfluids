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
//! \file GridCpuStrategy.cpp
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "GridCpuStrategy.h"

#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <vector>
#include <iostream>

#include "grid/distributions/Distribution.h"
#include "grid/GridImp.h"
#include "grid/NodeValues.h"

void GridCpuStrategy::allocateGridMemory(SPtr<GridImp> grid)
{
    grid->neighborIndexX        = new int[grid->size];
    grid->neighborIndexY        = new int[grid->size];
    grid->neighborIndexZ        = new int[grid->size];
    grid->neighborIndexNegative = new int[grid->size];

    grid->sparseIndices = new int[grid->size];

	grid->qIndices = new uint[grid->size];
	for (size_t i = 0; i < grid->size; i++) 
		grid->qIndices[i] = INVALID_INDEX;
}

void GridCpuStrategy::initalNodesToOutOfGrid(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (int index = 0; index < grid->size; index++)
        grid->initalNodeToOutOfGrid(index);
}

void GridCpuStrategy::allocateFieldMemory(Field* field)
{
    field->field = new char[field->size];
}


void GridCpuStrategy::findInnerNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (int index = 0; index < grid->size; index++)
        grid->findInnerNode(index);
}

void GridCpuStrategy::findEndOfGridStopperNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
	for (int index = 0; index < grid->size; index++)
		grid->findEndOfGridStopperNode(index);
}

void GridCpuStrategy::findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Find sparse indices...";

    coarseGrid->updateSparseIndices();
    findForNeighborsNewIndices(coarseGrid);
    if (fineGrid)
    {
        fineGrid->updateSparseIndices();
    }

    const uint newGridSize = coarseGrid->getSparseSize();
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "... done. new size: " << newGridSize << ", delete nodes:" << coarseGrid->getSize() - newGridSize << "\n";
}


void GridCpuStrategy::findForNeighborsNewIndices(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (int index = 0; index < grid->getSize(); index++)
        grid->setNeighborIndices(index);
}

void GridCpuStrategy::freeFieldMemory(Field* field)
{
    delete[] field->field;
}

void GridCpuStrategy::freeMemory(SPtr<GridImp> grid)
{
    if( grid->neighborIndexX        != nullptr ) { delete[] grid->neighborIndexX;        grid->neighborIndexX        = nullptr; }
    if( grid->neighborIndexY        != nullptr ) { delete[] grid->neighborIndexY;        grid->neighborIndexY        = nullptr; }
    if( grid->neighborIndexZ        != nullptr ) { delete[] grid->neighborIndexZ;        grid->neighborIndexZ        = nullptr; }
    if( grid->neighborIndexNegative != nullptr ) { delete[] grid->neighborIndexNegative; grid->neighborIndexNegative = nullptr; }
    if( grid->sparseIndices         != nullptr ) { delete[] grid->sparseIndices;         grid->sparseIndices         = nullptr; }
	if( grid->qIndices              != nullptr ) { delete[] grid->qIndices;              grid->qIndices              = nullptr; }
	if( grid->qValues               != nullptr ) { delete[] grid->qValues;               grid->qValues               = nullptr; }
	if( grid->qPatches              != nullptr ) { delete[] grid->qPatches;              grid->qPatches              = nullptr; }
}


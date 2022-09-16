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
//! \file GridInterface.h
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef GRID_INTERFACE_H
#define GRID_INTERFACE_H

#include "gpu/GridGenerator/global.h"

class GridImp;

class GridInterface
{
public:
    GRIDGENERATOR_EXPORT GridInterface();
    GRIDGENERATOR_EXPORT ~GridInterface();

    void GRIDGENERATOR_EXPORT findInterfaceCF(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid);
    void GRIDGENERATOR_EXPORT findBoundaryGridInterfaceCF(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid);


    void GRIDGENERATOR_EXPORT findInterfaceCF_GKS(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid);

    void GRIDGENERATOR_EXPORT findInterfaceFC(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid);
    void GRIDGENERATOR_EXPORT findOverlapStopper(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid);
    
    void GRIDGENERATOR_EXPORT findInvalidBoundaryNodes(const uint& indexOnCoarseGrid, GridImp* coarseGrid);

    void GRIDGENERATOR_EXPORT findForGridInterfaceSparseIndexCF(GridImp* coarseGrid, GridImp* fineGrid, uint index);
    void GRIDGENERATOR_EXPORT findForGridInterfaceSparseIndexFC(GridImp* coarseGrid, GridImp* fineGrid, uint index);

    void GRIDGENERATOR_EXPORT repairGridInterfaceOnMultiGPU(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid);

    void GRIDGENERATOR_EXPORT print() const;

    struct Interface
    {
        uint *fine, *coarse;
        uint numberOfEntries = 0;
        uint *offset;
    } fc{}, cf{};


private:
    uint getCoarseToFineIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid);
    bool isNeighborFineInvalid(real x, real y, real z, const GridImp* coarseGrid, const GridImp* fineGrid);

    uint getFineToCoarseIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid);

    static void findSparseIndex(uint* indices, GridImp* grid, uint index);

    uint findOffsetCF( const uint& indexOnCoarseGrid, GridImp* coarseGrid, uint interfaceIndex );

    uint findOffsetFC( const uint& indexOnCoarseGrid, GridImp* coarseGrid, uint interfaceIndex );
};


#endif
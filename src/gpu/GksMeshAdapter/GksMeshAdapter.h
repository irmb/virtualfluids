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
//! \file GksMeshAdapter.h
//! \ingroup GksMeshAdapter
//! \author Stephan Lenz
//=======================================================================================
#ifndef GKS_MESH_ADAPTER_H
#define GKS_MESH_ADAPTER_H

#include <memory>
#include <array>
#include <vector>

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "MeshCell.h"
#include "MeshFace.h"

#include "VirtualFluidsDefinitions.h"

#include "GksMeshAdapter_export.h"

class MultipleGridBuilder;

class GKSMESHADAPTER_EXPORT GksMeshAdapter{

public:

    SPtr<MultipleGridBuilder> gridBuilder;

    uint numberOfLevels;

    double dxCoarse;

    //////////////////////////////////////////////////////////////////////////

    std::vector<Vec3> nodes;

    //////////////////////////////////////////////////////////////////////////
    //
    //    C e l l    s o r t i n g :
    //
    //  | Level 0                   | Level 1                   | Level 2                   |
    //  | FluidCells   | GhostCells | FluidCells   | GhostCells | FluidCells   | GhostCells | 
    //
    //  GhostCells: not included in Cell update, i.e. BoundaryCells and FCC-Cells
    //  FluidCells: all other, including CFF-Cells
    //
    std::vector<MeshCell> cells;

    std::vector<uint> numberOfCellsPerLevel;
    std::vector<uint> numberOfBulkCellsPerLevel;
    std::vector<uint> startOfCellsPerLevel;

    //////////////////////////////////////////////////////////////////////////
    //
    //    F a c e    s o r t i n g :
    //
    //  | Level 0                              | Level 1                              | Level 2                              |
    //  | x-normal   | y-normal   | z-normal   | x-normal   | y-normal   | z-normal   | x-normal   | y-normal   | z-normal   |
    //
    std::vector<MeshFace> faces;

    std::vector<uint> numberOfFacesPerLevelXYZ;
    std::vector<uint> startOfFacesPerLevelXYZ;

    std::vector<uint> numberOfInnerFacesPerLevel;

    //////////////////////////////////////////////////////////////////////////
    // 
    // Connectivity from LBM grid to GKS Mesh
    //
    //    cellIdx = gridToMesh[ level ][ gridIdx ]
    //
    std::vector< std::vector<uint> > gridToMesh;

    //////////////////////////////////////////////////////////////////////////

    std::vector< uint_2 > periodicBoundaryNeighbors;

    //////////////////////////////////////////////////////////////////////////

public:

    GksMeshAdapter( SPtr<MultipleGridBuilder> gridBuilder );

    void inputGrid();


    void findCellToCellConnectivity();

    void countCells();

    void partitionCells();

    void refreshCellConnectivity(const std::vector<uint>& idxMap);

    void generateNodes();

    void computeCellGeometry();

    void generateFaces();

    void sortFaces();

    void countFaces();
    
    void findPeriodicBoundaryNeighbors();

    //////////////////////////////////////////////////////////////////////////

    double getDx(uint level);
};


#endif
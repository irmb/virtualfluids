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
//! \file MeshCell.h
//! \ingroup GksMeshAdapter
//! \author Stephan Lenz
//=======================================================================================
#ifndef MESH_CELL_H
#define MESH_CELL_H

#include <array>

#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/ArrayTypes.h"

#include "GksMeshAdapter_export.h"

struct GKSMESHADAPTER_EXPORT MeshCell{

    uint level;
    uint gridIdx;

    //////////////////////////////////////////////////////////////////////////

    uint_8 cellToNode;           uint_12 cellToEdgeNode;         uint_6 cellToFaceNode;

    // for sorting see LBM f numbering
    //                                   8                              
    //     4 o--------o 0            o--------o                    o--------o 
    //      /|       /|            7/|      4/|                   /|  4    /|  
    //     / |      / |            /3| 11   / |                  / |    2 / |  
    //  6 o--------o 2|           o--------o  |0                o--------o  |  
    //    |  |     |  |           |  |  10 |  |                 |1 |     | 0|  
    //    |5 o-----|--o 1        1|  o-----|--o                 |  o-----|--o  
    //    | /      | /            |5/     2| /6                 | / 3    | /   
    //    |/       |/             |/   9   |/                   |/    5  |/    
    //  7 o--------o 3            o--------o                    o--------o      
    //
    //  z | / y
    //    |/
    //    +---- x
    //     

    //////////////////////////////////////////////////////////////////////////

    uint_27 cellToCell;

    // for sorting see LBM f numbering

    //////////////////////////////////////////////////////////////////////////

    uint_8 children;

    // for sorting see cellToNode

    uint  parent;

    //////////////////////////////////////////////////////////////////////////

    Vec3   cellCenter;

    //////////////////////////////////////////////////////////////////////////

    // order: +x, -x, +y, -y, +z, -z (see cellToCell)
    bool_6 faceExists;
    
    bool isGhostCell;   // this denotes cells that do not have all neighbors, excluding coarse ghost cells

    bool isWall;

    bool isFluxBC;

    bool isInsulated;

    bool isRecvCell;

    char type;

    MeshCell();

    bool isCoarseGhostCell();

    bool isFineGhostCell();
};


#endif
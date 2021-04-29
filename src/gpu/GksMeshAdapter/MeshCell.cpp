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
//! \file MeshCell.cpp
//! \ingroup GksMeshAdapter
//! \author Stephan Lenz
//=======================================================================================
#include "MeshCell.h"

#include "GridGenerator/grid/NodeValues.h"

using namespace vf::gpu;

MeshCell::MeshCell(){

    level   = INVALID_INDEX;
    gridIdx = INVALID_INDEX;


    for( uint& index : this->cellToNode     ) index = INVALID_INDEX;
    for( uint& index : this->cellToEdgeNode ) index = INVALID_INDEX;
    for( uint& index : this->cellToFaceNode ) index = INVALID_INDEX;
    for( uint& index : this->cellToCell     ) index = INVALID_INDEX;
    for( uint& index : this->children       ) index = INVALID_INDEX;
    
    parent = INVALID_INDEX;

    for( bool& flag : this->faceExists    ) flag = false;

    isGhostCell = false;

    isWall = false;

    isFluxBC = false;

    isInsulated = false;

    isRecvCell = false;
}

bool MeshCell::isCoarseGhostCell()
{
    return this->type == FLUID_FCC;
}

bool MeshCell::isFineGhostCell()
{
    return this->type == FLUID_CFF;
}

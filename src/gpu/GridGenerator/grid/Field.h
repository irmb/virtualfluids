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
//! \file Field.h
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz, Martin Schoenherr
//=======================================================================================
#ifndef FIELD_H
#define FIELD_H

#include "gpu/GridGenerator/global.h"

struct Vertex;

class Field : public enableSharedFromThis<Field>
{
public:
    Field(uint size);
    Field() = default;
    void allocateMemory();
    void freeMemory();

    uint getSize() const;
    char getFieldEntry(uint index) const;

    bool is(uint index, char type) const;
    bool isCoarseToFineNode(uint index) const;
    bool isFineToCoarseNode(uint index) const;
    bool isFluid(uint index) const;
    bool isInvalidSolid(uint index) const;
    bool isQ(uint index) const;
    bool isBoundaryConditionNode(uint index) const;
    bool isInvalidCoarseUnderFine(uint index) const;
    bool isStopperOutOfGrid(uint index) const;
    bool isStopperCoarseUnderFine(uint index) const;
    bool isStopperSolid(uint index) const;
    bool isStopper(uint index) const;
    bool isInvalidOutOfGrid(uint index) const;

    void setFieldEntry(uint index, char val);
    void setFieldEntryToFluid(uint index);
    void setFieldEntryToInvalidSolid(uint index);
    void setFieldEntryToStopperOutOfGrid(uint index);
    void setFieldEntryToStopperOutOfGridBoundary(uint index);
    void setFieldEntryToStopperCoarseUnderFine(uint index);
    void setFieldEntryToInvalidCoarseUnderFine(uint index);
    void setFieldEntryToInvalidOutOfGrid(uint index);

protected:
    char *field = nullptr;
    uint size;
};

#endif

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
//=======================================================================================
#include "Field.h"

#include "grid/NodeValues.h"

using namespace vf::gpu;

Field::Field(uint size) : size(size)
{
}

void Field::allocateMemory()
{
    this->field = new char[this->size];
}

void Field::freeMemory()
{
    delete[] this->field;
}

// --------------------------------------------------------- //
//                        Getter                             //
// --------------------------------------------------------- //
uint Field::getSize() const
{
    return this->size;
}

char Field::getFieldEntry(uint index) const
{
    return this->field[index];
}

// --------------------------------------------------------- //
//                           Is                              //
// --------------------------------------------------------- //
bool Field::is(uint index, char type) const
{
    return field[index] == type;
}

bool Field::isCoarseToFineNode(uint index) const
{
    return field[index] == FLUID_CFC;
}

bool Field::isFineToCoarseNode(uint index) const
{
    return field[index] == FLUID_FCC;
}

bool Field::isFluid(uint index) const
{
    const char type = field[index];
    return type == FLUID || type == FLUID_CFC || type == FLUID_CFF || type == FLUID_FCC || type == FLUID_FCF || isBoundaryConditionNode(index);
}

bool Field::isInvalidSolid(uint index) const
{
    return field[index] == INVALID_SOLID;
}

bool Field::isInvalidOutOfGrid(uint index) const
{
    return field[index] == INVALID_OUT_OF_GRID;
}

bool Field::isInvalidCoarseUnderFine(uint index) const
{
    return field[index] == INVALID_COARSE_UNDER_FINE;
}

bool Field::isStopperOutOfGrid(uint index) const
{
    return field[index] == STOPPER_OUT_OF_GRID;
}

bool Field::isStopperCoarseUnderFine(uint index) const
{
    return field[index] == STOPPER_COARSE_UNDER_FINE;
}

bool Field::isStopperSolid(uint index) const
{
    return field[index] == STOPPER_SOLID;
}

bool Field::isStopper(uint index) const
{
    return isStopperOutOfGrid(index) || isStopperCoarseUnderFine(index) || isStopperSolid(index) || is(index, STOPPER_OUT_OF_GRID_BOUNDARY);
}

bool Field::isQ(uint index) const
{
    return field[index] == Q_DEPRECATED;
}

bool Field::isBoundaryConditionNode(uint index) const
{
    return  field[index] == BC_SOLID || field[index] == BC_OUTFLOW || field[index] == BC_VELOCITY || field[index] == BC_PRESSURE || field[index] == BC_SLIP || field[index] == BC_STRESS;
}

// --------------------------------------------------------- //
//                        Setter                             //
// --------------------------------------------------------- //
void Field::setFieldEntry(uint index, char val)
{
    this->field[index] = val;
}

void Field::setFieldEntryToFluid(uint index)
{
    this->field[index] = FLUID;
}

void Field::setFieldEntryToInvalidSolid(uint index)
{
    this->field[index] = INVALID_SOLID;
}

void Field::setFieldEntryToStopperOutOfGrid(uint index)
{
    this->field[index] = STOPPER_OUT_OF_GRID;
}

void Field::setFieldEntryToStopperOutOfGridBoundary(uint index)
{
    this->field[index] = STOPPER_OUT_OF_GRID_BOUNDARY;
}

void Field::setFieldEntryToStopperCoarseUnderFine(uint index)
{
    this->field[index] = STOPPER_COARSE_UNDER_FINE;
}

void Field::setFieldEntryToInvalidCoarseUnderFine(uint index)
{
    this->field[index] = INVALID_COARSE_UNDER_FINE;
}

void Field::setFieldEntryToInvalidOutOfGrid(uint index)
{
    this->field[index] = INVALID_OUT_OF_GRID;
}

//! \}

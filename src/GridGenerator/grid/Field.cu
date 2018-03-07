#include "Field.h"

#include <GridGenerator/grid/NodeValues.h>
#include <GridGenerator/grid/GridStrategy/GridStrategy.h>


HOST Field::Field(SPtr<GridStrategy> gridStrategy, uint size) : gridStrategy(gridStrategy), size(size)
{
    
}

Field::Field()
{
    
}

Field::~Field()
{
    
}

HOST void Field::allocateMemory()
{
    gridStrategy->allocateFieldMemory(this);
}

HOST void Field::freeMemory()
{
    gridStrategy->freeFieldMemory(this);
}

// --------------------------------------------------------- //
//                        Getter                             //
// --------------------------------------------------------- //
HOSTDEVICE uint Field::getSize() const
{
    return this->size;
}

HOSTDEVICE char Field::getFieldEntry(uint index) const
{
    return this->field[index];
}

// --------------------------------------------------------- //
//                           Is                              //
// --------------------------------------------------------- //
HOSTDEVICE bool Field::is(uint index, char type) const
{
    return field[index] == type;
}

HOSTDEVICE bool Field::isCoarseToFineNode(uint index) const
{
    return field[index] == FLUID_CFC;
}

HOSTDEVICE bool Field::isFineToCoarseNode(uint index) const
{
    return field[index] == FLUID_FCC;
}

HOSTDEVICE bool Field::isFluid(uint index) const
{
    const char type = field[index];
    return type == FLUID || type == FLUID_CFC || type == FLUID_CFF || type == FLUID_FCC || type == FLUID_FCF;
}

HOSTDEVICE bool Field::isSolid(uint index) const
{
    return field[index] == SOLID;
}

HOSTDEVICE bool Field::isOutOfGrid(uint index) const
{
    return field[index] == OUT_OF_GRID;
}

HOSTDEVICE bool Field::isInvalid(uint index) const
{
    return field[index] == INVALID_NODE;
}

HOSTDEVICE bool Field::isStopperEndOfGrid(uint index) const
{
    return field[index] == STOPPER_END_OF_GRID;
}

HOSTDEVICE bool Field::isStopperOverlapGrid(uint index) const
{
    return field[index] == STOPPER_OVERLAP_GRID;
}

HOSTDEVICE bool Field::isQ(uint index) const
{
    return field[index] == Q;
}

HOSTDEVICE bool Field::isRb(uint index) const
{
    return field[index] == VELOCITY || field[index] == PRESSURE || field[index] == NOSLIP || field[index] == SOLID;
}

// --------------------------------------------------------- //
//                        Setter                             //
// --------------------------------------------------------- //
HOSTDEVICE void Field::setFieldEntry(uint index, char val)
{
    this->field[index] = val;
}

HOSTDEVICE void Field::setFieldEntryToFluid(uint index)
{
    this->field[index] = FLUID;
}

HOSTDEVICE void Field::setFieldEntryToSolid(uint index)
{
    this->field[index] = SOLID;
}

HOSTDEVICE void Field::setFieldEntryToStopperEndOfGrid(uint index)
{
    this->field[index] = STOPPER_END_OF_GRID;
}

HOSTDEVICE void Field::setFieldEntryToStopperOverlapGrid(uint index)
{
    this->field[index] = STOPPER_OVERLAP_GRID;
}

HOSTDEVICE void Field::setFieldEntryToInvalid(uint index)
{
    this->field[index] = INVALID_NODE;
}

HOSTDEVICE void Field::setFieldEntryToOutOfGrid(uint index)
{
    this->field[index] = OUT_OF_GRID;
}

#include "Field.h"

#include "grid/NodeValues.h"
#include "grid/GridStrategy/GridStrategy.h"

using namespace VF::GPU;

CUDA_HOST Field::Field(SPtr<GridStrategy> gridStrategy, uint size) : gridStrategy(gridStrategy), size(size)
{
    
}

Field::Field()
{
    
}

Field::~Field()
{
    
}

CUDA_HOST void Field::allocateMemory()
{
    gridStrategy->allocateFieldMemory(this);
}

CUDA_HOST void Field::freeMemory()
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
    return type == FLUID || type == FLUID_CFC || type == FLUID_CFF || type == FLUID_FCC || type == FLUID_FCF || isBoundaryConditionNode(index);
}

HOSTDEVICE bool Field::isInvalidSolid(uint index) const
{
    return field[index] == INVALID_SOLID;
}

HOSTDEVICE bool Field::isInvalidOutOfGrid(uint index) const
{
    return field[index] == INVALID_OUT_OF_GRID;
}

HOSTDEVICE bool Field::isInvalidCoarseUnderFine(uint index) const
{
    return field[index] == INVALID_COARSE_UNDER_FINE;
}

HOSTDEVICE bool Field::isStopperOutOfGrid(uint index) const
{
    return field[index] == STOPPER_OUT_OF_GRID;
}

HOSTDEVICE bool Field::isStopperCoarseUnderFine(uint index) const
{
    return field[index] == STOPPER_COARSE_UNDER_FINE;
}

HOSTDEVICE bool Field::isStopperSolid(uint index) const
{
	return field[index] == STOPPER_SOLID;
}

HOSTDEVICE bool Field::isStopper(uint index) const
{
    return isStopperOutOfGrid(index) || isStopperCoarseUnderFine(index) || isStopperSolid(index) || is(index, STOPPER_OUT_OF_GRID_BOUNDARY);
}

HOSTDEVICE bool Field::isQ(uint index) const
{
    return field[index] == Q_DEPRECATED;
}

HOSTDEVICE bool Field::isBoundaryConditionNode(uint index) const
{
    return  field[index] == BC_SOLID || field[index] == BC_OUTFLOW || field[index] == BC_VELOCITY || field[index] == BC_PRESSURE || field[index] == BC_SLIP;
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

HOSTDEVICE void Field::setFieldEntryToInvalidSolid(uint index)
{
    this->field[index] = INVALID_SOLID;
}

HOSTDEVICE void Field::setFieldEntryToStopperOutOfGrid(uint index)
{
    this->field[index] = STOPPER_OUT_OF_GRID;
}

HOSTDEVICE void Field::setFieldEntryToStopperOutOfGridBoundary(uint index)
{
    this->field[index] = STOPPER_OUT_OF_GRID_BOUNDARY;
}

HOSTDEVICE void Field::setFieldEntryToStopperCoarseUnderFine(uint index)
{
    this->field[index] = STOPPER_COARSE_UNDER_FINE;
}

HOSTDEVICE void Field::setFieldEntryToInvalidCoarseUnderFine(uint index)
{
    this->field[index] = INVALID_COARSE_UNDER_FINE;
}

HOSTDEVICE void Field::setFieldEntryToInvalidOutOfGrid(uint index)
{
    this->field[index] = INVALID_OUT_OF_GRID;
}

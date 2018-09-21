#include "BoundaryCondition.h"

#include "Side.h"

bool BoundaryCondition::isSide( SideType side ) const
{
    return this->side->whoAmI() == side;
}
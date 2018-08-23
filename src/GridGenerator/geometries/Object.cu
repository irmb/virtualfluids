#include "Object.h"
#include "grid/GridImp.h"
#include "grid/GridStrategy/GridStrategy.h"

void Object::findInnerNodes(SPtr<GridImp> grid)
{
    grid->getGridStrategy()->findInnerNodes( grid );
}

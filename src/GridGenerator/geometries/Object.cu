#include "Object.h"
#include "grid/GridImp.h"
#include "grid/GridStrategy/GridStrategy.h"

void Object::findInnerNodes(SPtr<GridImp> grid)
{
    grid->getGridStrategy()->findInnerNodes( grid );
}

HOST int Object::getIntersection(const Vertex & P, const Vertex & direction, Vertex & pointOnObject, real & qVal)
{
    return 1;
}

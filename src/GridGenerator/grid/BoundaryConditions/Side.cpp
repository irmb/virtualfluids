#include "Side.h"

#include "../Grid.h"
#include "../NodeValues.h"

#include "BoundaryCondition.h"

void Side::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::string coord, real constant,
                      real startInner, real endInner, real startOuter, real endOuter)
{
    for (real v1 = startInner; v1 < endInner; v1 += grid->getDelta())
    {
        for (real v2 = startOuter; v2 < endOuter; v2 += grid->getDelta())
        {
            uint index = getIndex(grid, coord, constant, v1, v2);
            if (grid->getFieldEntry(index) == FLUID)
                boundaryCondition->indices.push_back(index);
        }
    }
}

uint Side::getIndex(SPtr<Grid> grid, std::string coord, real constant, real v1, real v2)
{
    if (coord == "x")
        return grid->transCoordToIndex(constant, v1, v2);
    if (coord == "y")
        return grid->transCoordToIndex(v1, constant, v2);
    if (coord == "z")
        return grid->transCoordToIndex(v1, v2, constant);
    return -1;
}

void MX::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityX(false);
}

void MX::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition)
{
    Side::addIndices(grid, boundaryCondition, "x", grid->getStartX() + grid->getDelta(), grid->getStartY(),
                     grid->getEndY(), grid->getStartZ(), grid->getEndZ());
}

void PX::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityX(false);
}

void PX::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition)
{
    Side::addIndices(grid, boundaryCondition, "x", grid->getEndX() - grid->getDelta(), grid->getStartY(),
                     grid->getEndY(), grid->getStartZ(), grid->getEndZ());
}

void MY::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityY(false);
}

void MY::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition)
{
    Side::addIndices(grid, boundaryCondition, "y", grid->getStartY() + grid->getDelta(), grid->getStartX(),
                     grid->getEndX(), grid->getStartZ(), grid->getEndZ());
}

void PY::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityY(false);
}

void PY::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition)
{
    Side::addIndices(grid, boundaryCondition, "y", grid->getEndY() - grid->getDelta(), grid->getStartX(),
                     grid->getEndX(), grid->getStartZ(), grid->getEndZ());
}

void MZ::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityZ(false);
}

void MZ::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition)
{
    Side::addIndices(grid, boundaryCondition, "z", grid->getStartZ() + grid->getDelta(), grid->getStartX(),
                     grid->getEndX(), grid->getStartY(), grid->getEndY());
}

void PZ::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityZ(false);
}

void PZ::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition)
{
    Side::addIndices(grid, boundaryCondition, "z", grid->getEndZ() - grid->getDelta(), grid->getStartX(),
                     grid->getEndX(), grid->getStartY(), grid->getEndY());
}

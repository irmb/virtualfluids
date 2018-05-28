#include "Side.h"

#include "../Grid.h"
#include "../NodeValues.h"

#include "BoundaryCondition.h"

void Side::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::string coord, real constant,
                      real startInner, real endInner, real startOuter, real endOuter)
{
    for (real v1 = startInner; v1 <= endInner; v1 += grid->getDelta())
    {
        for (real v2 = startOuter; v2 <= endOuter; v2 += grid->getDelta())
        {
            uint index = getIndex(grid, coord, constant, v1, v2);
            if (grid->getFieldEntry(index) == FLUID)
                boundaryCondition->indices.push_back(grid->getSparseIndex(index) + 1);
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

void Geometry::setPeriodicy(SPtr<Grid> grid)
{

}

void Geometry::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::map<SideType, bool> sideIsSet)
{

}

void MX::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityX(false);
}

void MX::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::map<SideType, bool> sideIsSet)
{
    real startInner = grid->getStartY();
    real endInner = grid->getEndY();

    real startOuter = grid->getStartZ();
    real endOuter = grid->getEndZ();

    if (sideIsSet[SideType::MY])
        startInner += grid->getDelta();

    if (sideIsSet[SideType::PY])
        endInner -= grid->getDelta();

    if (sideIsSet[SideType::MZ])
        startOuter += grid->getDelta();

    if (sideIsSet[SideType::PZ])
        endOuter -= grid->getDelta();


    Side::addIndices(grid, boundaryCondition, "x", grid->getStartX() + grid->getDelta(), startInner,
        endInner, startOuter, endOuter);
}

void PX::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityX(false);
}

void PX::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::map<SideType, bool> sideIsSet)
{
    real startInner = grid->getStartY();
    real endInner = grid->getEndY();

    real startOuter = grid->getStartZ();
    real endOuter = grid->getEndZ();

    if (sideIsSet[SideType::MY])
        startInner += grid->getDelta();

    if (sideIsSet[SideType::PY])
        endInner -= grid->getDelta();

    if (sideIsSet[SideType::MZ])
        startOuter += grid->getDelta();

    if (sideIsSet[SideType::PZ])
        endOuter -= grid->getDelta();

    Side::addIndices(grid, boundaryCondition, "x", grid->getEndX() - grid->getDelta(), startInner,
        endInner, startOuter, endOuter);
}

void MY::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityY(false);
}

void MY::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::map<SideType, bool> sideIsSet)
{
    real startInner = grid->getStartX();
    real endInner = grid->getEndX();

    real startOuter = grid->getStartZ();
    real endOuter = grid->getEndZ();

    if (sideIsSet[SideType::MX])
        startInner += grid->getDelta();

    if (sideIsSet[SideType::PX])
        endInner -= grid->getDelta();

    if (sideIsSet[SideType::MZ])
        startOuter += grid->getDelta();

    if (sideIsSet[SideType::PZ])
        endOuter -= grid->getDelta();

    Side::addIndices(grid, boundaryCondition, "y", grid->getStartY() + grid->getDelta(), startInner,
        endInner, startOuter, endOuter);
}

void PY::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityY(false);
}

void PY::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::map<SideType, bool> sideIsSet)
{
    real startInner = grid->getStartX();
    real endInner = grid->getEndX();

    real startOuter = grid->getStartZ();
    real endOuter = grid->getEndZ();

    if (sideIsSet[SideType::MX])
        startInner += grid->getDelta();

    if (sideIsSet[SideType::PX])
        endInner -= grid->getDelta();

    if (sideIsSet[SideType::MZ])
        startOuter += grid->getDelta();

    if (sideIsSet[SideType::PZ])
        endOuter -= grid->getDelta();

    Side::addIndices(grid, boundaryCondition, "y", grid->getEndY() - grid->getDelta(), startInner,
        endInner, startOuter, endOuter);
}

void MZ::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityZ(false);
}

void MZ::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::map<SideType, bool> sideIsSet)
{
    real startInner = grid->getStartX();
    real endInner = grid->getEndX();

    real startOuter = grid->getStartY();
    real endOuter = grid->getEndY();

    if (sideIsSet[SideType::MX])
        startInner += grid->getDelta();

    if (sideIsSet[SideType::PX])
        endInner -= grid->getDelta();

    if (sideIsSet[SideType::MY])
        startOuter += grid->getDelta();

    if (sideIsSet[SideType::PY])
        endOuter -= grid->getDelta();

    Side::addIndices(grid, boundaryCondition, "z", grid->getStartZ() + grid->getDelta(), startInner,
        endInner, startOuter, endOuter);
}

void PZ::setPeriodicy(SPtr<Grid> grid)
{
    grid->setPeriodicityZ(false);
}

void PZ::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::map<SideType, bool> sideIsSet)
{
    real startInner = grid->getStartX();
    real endInner = grid->getEndX();

    real startOuter = grid->getStartY();
    real endOuter = grid->getEndY();

    if (sideIsSet[SideType::MX])
        startInner += grid->getDelta();

    if (sideIsSet[SideType::PX])
        endInner -= grid->getDelta();

    if (sideIsSet[SideType::MY])
        startOuter += grid->getDelta();

    if (sideIsSet[SideType::PY])
        endOuter -= grid->getDelta();

    Side::addIndices(grid, boundaryCondition, "z", grid->getEndZ() - grid->getDelta(), startInner,
        endInner, startOuter, endOuter);
}

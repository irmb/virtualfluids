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
            const uint index = getIndex(grid, coord, constant, v1, v2);
            if (grid->getFieldEntry(index) == FLUID)
            {
                boundaryCondition->indices.push_back(grid->getSparseIndex(index) + 1);
                setPressureNeighborIndices(boundaryCondition, grid, index);
            }
        }
    }
}

void Side::setPressureNeighborIndices(SPtr<BoundaryCondition> boundaryCondition, SPtr<Grid> grid, const uint index)
{
    auto pressureBoundaryCondition = std::dynamic_pointer_cast<PressureBoundaryCondition>(boundaryCondition);
    if (pressureBoundaryCondition)
    {
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);

        real nx = x;
        real ny = y;
        real nz = z;

        if (boundaryCondition->side->getCoordinate() == X_INDEX)
            nx = boundaryCondition->side->getDirection() * grid->getDelta() + x;
        if (boundaryCondition->side->getCoordinate() == Y_INDEX)
            ny = boundaryCondition->side->getDirection() * grid->getDelta() + y;
        if (boundaryCondition->side->getCoordinate() == Z_INDEX)
            nz = boundaryCondition->side->getDirection() * grid->getDelta() + z;

        int neighborIndex = grid->transCoordToIndex(nx, ny, nz);
        int sparseIndex = grid->getSparseIndex(neighborIndex);
        pressureBoundaryCondition->neighborIndices.push_back(sparseIndex);
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


void Geometry::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::map<SideType, bool> sideIsSet)
{
    auto geometryBoundaryCondition = std::dynamic_pointer_cast<GeometryBoundaryCondition>(boundaryCondition);

    std::vector<real> qNode(27);
    bool qFound = false;

    for (int i = 0; i < grid->getSize(); i++)
    {
        for (int dir = 0; dir < grid->getEndDirection(); dir++)
        {
            const int qIndex = dir * grid->getSize() + i;
            const real q = grid->getDistribution()[qIndex];

            qNode[dir] = q;
            if (q > 0)
                qFound = true;

        }

        if (qFound)
        {
            geometryBoundaryCondition->indices.push_back(i);
            geometryBoundaryCondition->qs.push_back(qNode);
        }

        qFound = false;
    }
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

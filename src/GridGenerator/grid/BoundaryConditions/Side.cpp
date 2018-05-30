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
                grid->setFieldEntry(index, boundaryCondition->getType());
                boundaryCondition->indices.push_back(index);
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
            nx = -boundaryCondition->side->getDirection() * grid->getDelta() + x;
        if (boundaryCondition->side->getCoordinate() == Y_INDEX)
            ny = -boundaryCondition->side->getDirection() * grid->getDelta() + y;
        if (boundaryCondition->side->getCoordinate() == Z_INDEX)
            nz = -boundaryCondition->side->getDirection() * grid->getDelta() + z;

        int neighborIndex = grid->transCoordToIndex(nx, ny, nz);
        pressureBoundaryCondition->neighborIndices.push_back(neighborIndex);
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


void Geometry::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    auto geometryBoundaryCondition = std::dynamic_pointer_cast<GeometryBoundaryCondition>(boundaryCondition);

    std::vector<real> qNode(27);
    bool qFound = false;

    for (int i = 0; i < grid[level]->getSize(); i++)
    {
        if (grid[level]->getFieldEntry(i) != BC_GEOMETRY)
            continue;

        for (int dir = 0; dir < grid[level]->getEndDirection(); dir++)
        {
            const int qIndex = dir * grid[level]->getSize() + i;
            const real q = grid[level]->getDistribution()[qIndex];

            qNode[dir] = q;
            if (q > 0)
                qFound = true;
            else
                qNode[dir] = -1.0;
        }

        if (qFound)
        {
            geometryBoundaryCondition->indices.push_back(i);
            geometryBoundaryCondition->qs.push_back(qNode);
        }

        qFound = false;
    }
}



void MX::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartY();
    real endInner = grid[level]->getEndY();

    real startOuter = grid[level]->getStartZ();
    real endOuter = grid[level]->getEndZ();

    Vertex exactCoords(grid[level]->getStartX(), grid[level]->getStartY() + (grid[level]->getEndY() - grid[level]->getStartY()) / 2.0, grid[level]->getStartZ() + (grid[level]->getEndZ() - grid[level]->getStartZ()) / 2.0);
    Vertex noodCoords = grid[level]->getMaximumOnNode(exactCoords);
    real coords[3] = { noodCoords.x, noodCoords.y, noodCoords.z};
    real startCoord = grid[level]->getFirstFluidNode(coords, level, grid[level]->getStartX());

    if(!isBoundaryOnFineGrid(level, grid, startCoord))
        return;

    Side::addIndices(grid[level], boundaryCondition, "x", startCoord, startInner,
        endInner, startOuter, endOuter);
}



bool MX::isBoundaryOnFineGrid(uint level, std::vector<SPtr<Grid>> grid, real startCoord)
{
    if (level > 0) {
        real coords[3] = { grid[level - 1]->getStartX(), grid[level - 1]->getStartY() + (grid[level - 1]->getEndY() - grid[level - 1]->getStartY()) / 2.0, grid[level - 1]->getStartZ() + (grid[level - 1]->getEndZ() - grid[level - 1]->getStartZ()) / 2.0 };
        real startCoordCoarser = grid[level - 1]->getFirstFluidNode(coords, level - 1, grid[level - 1]->getStartX());
        if (startCoord > startCoordCoarser)
            return false;
    }
}


void PX::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    //real startInner = grid->getStartY();
    //real endInner = grid->getEndY();

    //real startOuter = grid->getStartZ();
    //real endOuter = grid->getEndZ();

    //real coords[3] = { grid->getEndX(), grid->getStartY() + (grid->getEndY() - grid->getStartY()) / 2.0, grid->getStartZ() + (grid->getEndZ() - grid->getStartZ()) / 2.0 };
    //real startCoord = grid->getLastFluidNode(coords, 0, grid->getEndX());

    //Side::addIndices(grid, boundaryCondition, "x", startCoord, startInner,
    //    endInner, startOuter, endOuter);
}

void MY::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    //real startInner = grid->getStartX();
    //real endInner = grid->getEndX();

    //real startOuter = grid->getStartZ();
    //real endOuter = grid->getEndZ();


    //Side::addIndices(grid, boundaryCondition, "y", grid->getStartY() + grid->getDelta(), startInner,
    //    endInner, startOuter, endOuter);
}


void PY::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    //real startInner = grid->getStartX();
    //real endInner = grid->getEndX();

    //real startOuter = grid->getStartZ();
    //real endOuter = grid->getEndZ();


    //Side::addIndices(grid, boundaryCondition, "y", grid->getEndY() - grid->getDelta(), startInner,
    //    endInner, startOuter, endOuter);
}


void MZ::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    //real startInner = grid->getStartX();
    //real endInner = grid->getEndX();

    //real startOuter = grid->getStartY();
    //real endOuter = grid->getEndY();

    //Side::addIndices(grid, boundaryCondition, "z", grid->getStartZ() + grid->getDelta(), startInner,
    //    endInner, startOuter, endOuter);
}

void PZ::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    //real startInner = grid->getStartX();
    //real endInner = grid->getEndX();

    //real startOuter = grid->getStartY();
    //real endOuter = grid->getEndY();


    //Side::addIndices(grid, boundaryCondition, "z", grid->getEndZ() - grid->getDelta(), startInner,
    //    endInner, startOuter, endOuter);
}

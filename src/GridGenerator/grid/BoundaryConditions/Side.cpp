#include "Side.h"

#include "../Grid.h"
#include "../NodeValues.h"
#include <GridGenerator/utilities/math/Math.h>

#include "BoundaryCondition.h"

void Side::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::string coord, real constant,
                      real startInner, real endInner, real startOuter, real endOuter)
{
    for (real v1 = startInner; v1 <= endInner; v1 += grid->getDelta())
    {
        for (real v2 = startOuter; v2 <= endOuter; v2 += grid->getDelta())
        {
            const uint index = getIndex(grid, coord, constant, v1, v2);

            if ((index != INVALID_INDEX) && (  grid->getFieldEntry(index) == FLUID
                                            || grid->getFieldEntry(index) == FLUID_CFC
                                            || grid->getFieldEntry(index) == FLUID_CFF
                                            || grid->getFieldEntry(index) == FLUID_FCC
                                            || grid->getFieldEntry(index) == FLUID_FCF ) )
            {
                grid->setFieldEntry(index, boundaryCondition->getType());
                boundaryCondition->indices.push_back(index);
                setPressureNeighborIndices(boundaryCondition, grid, index);

                setQs(grid, boundaryCondition, index);

                boundaryCondition->patches.push_back(0);
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

void Side::setQs(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, uint index)
{

    std::vector<real> qNode(grid->getEndDirection() + 1);

    for (int dir = 0; dir <= grid->getEndDirection(); dir++)
    {
        real x,y,z;
        grid->transIndexToCoords( index, x, y, z );

        x += grid->getDirection()[dir * DIMENSION + 0] * grid->getDelta();
        y += grid->getDirection()[dir * DIMENSION + 1] * grid->getDelta();
        z += grid->getDirection()[dir * DIMENSION + 2] * grid->getDelta();

        uint neighborIndex = grid->transCoordToIndex( x, y, z );

        if( grid->getFieldEntry(neighborIndex) == STOPPER_OUT_OF_GRID_BOUNDARY )
            qNode[dir] = 0.5;
        else
            qNode[dir] = -1.0;
    }

    boundaryCondition->qs.push_back(qNode);
}

uint Side::getIndex(SPtr<Grid> grid, std::string coord, real constant, real v1, real v2)
{
    if (coord == "x")
        return grid->transCoordToIndex(constant, v1, v2);
    if (coord == "y")
        return grid->transCoordToIndex(v1, constant, v2);
    if (coord == "z")
        return grid->transCoordToIndex(v1, v2, constant);
    return INVALID_INDEX;
}


void Geometry::addIndices(std::vector<SPtr<Grid> > grids, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    auto geometryBoundaryCondition = std::dynamic_pointer_cast<GeometryBoundaryCondition>(boundaryCondition);

    std::vector<real> qNode(grids[level]->getEndDirection() + 1);
    bool qFound = false;

    for (uint index = 0; index < grids[level]->getSize(); index++)
    {
        if (grids[level]->getFieldEntry(index) != BC_SOLID)
            continue;

        for (int dir = 0; dir <= grids[level]->getEndDirection(); dir++)
        {
			const real q = grids[level]->getQValue(index, dir);

            qNode[dir] = q;

            // also the neighbor if any Qs are required
            real x,y,z;
            grids[level]->transIndexToCoords( index, x, y, z );

            x += grids[level]->getDirection()[dir * DIMENSION + 0] * grids[level]->getDelta();
            y += grids[level]->getDirection()[dir * DIMENSION + 1] * grids[level]->getDelta();
            z += grids[level]->getDirection()[dir * DIMENSION + 2] * grids[level]->getDelta();

            uint neighborIndex = grids[level]->transCoordToIndex( x, y, z );

            if( grids[level]->getFieldEntry(neighborIndex) == STOPPER_OUT_OF_GRID_BOUNDARY )
                qNode[dir] = 0.5;
        }

        geometryBoundaryCondition->indices.push_back(index);
        geometryBoundaryCondition->qs.push_back(qNode);
        geometryBoundaryCondition->patches.push_back( grids[level]->getQPatch(index) );

        qFound = false;
    }
}



void MX::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartY();
    real endInner = grid[level]->getEndY();

    real startOuter = grid[level]->getStartZ();
    real endOuter = grid[level]->getEndZ();

    real coordinateNormal = grid[level]->getStartX() + grid[level]->getDelta();

    if( coordinateNormal > grid[0]->getStartX() + grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "x", coordinateNormal, startInner, endInner, startOuter, endOuter);

}

void PX::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartY();
    real endInner = grid[level]->getEndY();

    real startOuter = grid[level]->getStartZ();
    real endOuter = grid[level]->getEndZ();

    real coordinateNormal = grid[level]->getEndX() - grid[level]->getDelta();

    if( coordinateNormal < grid[0]->getEndX() - grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "x", coordinateNormal, startInner, endInner, startOuter, endOuter);
}

void MY::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartX();
    real endInner = grid[level]->getEndX();

    real startOuter = grid[level]->getStartZ();
    real endOuter = grid[level]->getEndZ();

    real coordinateNormal = grid[level]->getStartY() + grid[level]->getDelta();

    if( coordinateNormal > grid[0]->getStartY() + grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "y", coordinateNormal, startInner, endInner, startOuter, endOuter);
}


void PY::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartX();
    real endInner = grid[level]->getEndX();

    real startOuter = grid[level]->getStartZ();
    real endOuter = grid[level]->getEndZ();

    real coordinateNormal = grid[level]->getEndY() - grid[level]->getDelta();

    if( coordinateNormal < grid[0]->getEndY() - grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "y", coordinateNormal, startInner, endInner, startOuter, endOuter);
}


void MZ::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartX();
    real endInner = grid[level]->getEndX();

    real startOuter = grid[level]->getStartY();
    real endOuter = grid[level]->getEndY();

    real coordinateNormal = grid[level]->getStartZ() + grid[level]->getDelta();

    if( coordinateNormal > grid[0]->getStartZ() + grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "z", coordinateNormal, startInner, endInner, startOuter, endOuter);
}

void PZ::addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartX();
    real endInner = grid[level]->getEndX();

    real startOuter = grid[level]->getStartY();
    real endOuter = grid[level]->getEndY();

    real coordinateNormal = grid[level]->getEndZ() - grid[level]->getDelta();

    if( coordinateNormal < grid[0]->getEndZ() - grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "z", coordinateNormal, startInner, endInner, startOuter, endOuter);
}

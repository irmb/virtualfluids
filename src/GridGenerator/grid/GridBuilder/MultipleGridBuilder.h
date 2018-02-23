#ifndef MULTIPLE_GRID_BUILDER_H
#define MULTIPLE_GRID_BUILDER_H

#include "GridGenerator/global.h"

#include <vector>
#include <array>

#include "LevelGridBuilder.h"


#include "../GridFactory.h"

class Object;


class MultipleGridBuilder : public LevelGridBuilder
{
private:
    VF_PUBLIC MultipleGridBuilder(SPtr<GridFactory> gridFactory, Device device = Device::CPU, const std::string &d3qxx = "D3Q27");

public:
    VF_PUBLIC static SPtr<MultipleGridBuilder> makeShared(SPtr<GridFactory> gridFactory);

    VF_PUBLIC void addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta);
    VF_PUBLIC void addGrid(Object* gridShape);
    VF_PUBLIC void addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ);
    VF_PUBLIC void addFineGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, uint level);

    void addFineGridToList( uint level, real startXfine, real startYfine, real startZfine, real endXfine, real endYfine, real endZFine);

   void addIntermediateGridsToList(uint levelDifference, uint levelFine, uint nodesBetweenGrids, real startXfine, real startYfine, real startZfine, real endXfine, real endYfine, real endZfine);
    void addGridsToListIfValid(uint oldSize);

    void addGridToListIfValid(SPtr<Grid> grid);

    VF_PUBLIC uint getNumberOfLevels() const;
    VF_PUBLIC real getDelta(uint level) const;

    VF_PUBLIC real getStartX(uint level) const;
    VF_PUBLIC real getStartY(uint level) const;
    VF_PUBLIC real getStartZ(uint level) const;

    VF_PUBLIC real getEndX(uint level) const;
    VF_PUBLIC real getEndY(uint level) const;
    VF_PUBLIC real getEndZ(uint level) const;

    VF_PUBLIC std::vector<SPtr<Grid> > getGrids() const;
    VF_PUBLIC void buildGrids();

private:
    void addGridToList(SPtr<Grid> grid);
    real calculateDelta(uint level) const;
    bool coarseGridExists() const;
    bool isGridInCoarseGrid(SPtr<Grid> grid) const;
    SPtr<Grid> makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, uint level) const;
    std::array<real, 6> getStaggeredCoordinates(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const;
    std::array<real, 3> getOffset(real delta) const;
    std::vector<uint> getSpacingFactors(uint levelDifference) const;

    SPtr<Grid> makeGrid(Object* gridShape, uint level) const;
    SPtr<Grid> makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const;

    static void emitNoCoarseGridExistsWarning();
    static void emitGridIsNotInCoarseGridWarning();

    //std::vector<SPtr<Grid> > grids;

    SPtr<GridFactory> gridFactory;
 
};

#endif


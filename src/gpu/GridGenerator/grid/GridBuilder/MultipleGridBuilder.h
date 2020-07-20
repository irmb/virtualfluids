#ifndef MULTIPLE_GRID_BUILDER_H
#define MULTIPLE_GRID_BUILDER_H

#include <vector>
#include <array>

#include "Core/LbmOrGks.h"

#include "global.h"

#include "grid/GridBuilder/LevelGridBuilder.h"
#include "grid/GridFactory.h"

class Object;
class BoundingBox;

class MultipleGridBuilder : public LevelGridBuilder
{
private:
    VIRTUALFLUIDS_GPU_EXPORT MultipleGridBuilder(SPtr<GridFactory> gridFactory, Device device = Device::CPU, const std::string &d3qxx = "D3Q27");

public:
    VIRTUALFLUIDS_GPU_EXPORT static SPtr<MultipleGridBuilder> makeShared(SPtr<GridFactory> gridFactory);

    VIRTUALFLUIDS_GPU_EXPORT void addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta);
    VIRTUALFLUIDS_GPU_EXPORT void addGrid(Object* gridShape);
    VIRTUALFLUIDS_GPU_EXPORT void addGrid(Object* gridShape, uint levelFine);

    VIRTUALFLUIDS_GPU_EXPORT void addGeometry(Object* gridShape);
    VIRTUALFLUIDS_GPU_EXPORT void addGeometry(Object* solidObject, uint level);

    VIRTUALFLUIDS_GPU_EXPORT uint getNumberOfLevels() const;
    VIRTUALFLUIDS_GPU_EXPORT real getDelta(uint level) const;

    VIRTUALFLUIDS_GPU_EXPORT real getStartX(uint level) const;
    VIRTUALFLUIDS_GPU_EXPORT real getStartY(uint level) const;
    VIRTUALFLUIDS_GPU_EXPORT real getStartZ(uint level) const;

    VIRTUALFLUIDS_GPU_EXPORT real getEndX(uint level) const;
    VIRTUALFLUIDS_GPU_EXPORT real getEndY(uint level) const;
    VIRTUALFLUIDS_GPU_EXPORT real getEndZ(uint level) const;

    VIRTUALFLUIDS_GPU_EXPORT std::vector<SPtr<Grid> > getGrids() const;
    VIRTUALFLUIDS_GPU_EXPORT void buildGrids(LbmOrGks lbmOrGks, bool enableThinWalls = false);

    VIRTUALFLUIDS_GPU_EXPORT void setNumberOfLayers( uint numberOfLayersFine, uint numberOfLayersBetweenLevels );

    VIRTUALFLUIDS_GPU_EXPORT void writeGridsToVtk(const std::string& path) const;

    VIRTUALFLUIDS_GPU_EXPORT void setSubDomainBox(SPtr<BoundingBox> subDomainBox);

private:
    void addGridToList(SPtr<Grid> grid);
    real calculateDelta(uint level) const;
    bool coarseGridExists() const;
    bool isGridInCoarseGrid(SPtr<Grid> grid) const;

    void addFineGridToList(uint level, Object* gridShape);
    void addIntermediateGridsToList(uint levelDifference, uint levelFine, uint nodesBetweenGrids, Object* gridShape);
    void eraseGridsFromListIfInvalid(uint oldSize);
    void addGridToListIfValid(SPtr<Grid> grid);

	std::array<real, 6> getStaggeredCoordinates(Object* gridShape, uint level, uint levelFine, bool& xOddStart, bool& yOddStart, bool& zOddStart) const;
    std::array<real, 6> getStaggeredCoordinates(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) const;
    std::array<real, 3> getOffset(real delta) const;
    std::vector<uint> getSpacingFactors(uint levelDifference) const;

    SPtr<Grid> makeGrid(Object* gridShape, uint level, uint levelFine);
    SPtr<Grid> makeGrid(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) const;

    static void emitNoCoarseGridExistsWarning();
    static void emitGridIsNotInCoarseGridWarning();

    //std::vector<SPtr<Grid> > grids;

    SPtr<GridFactory> gridFactory;
    Object* solidObject;

    uint numberOfLayersFine;
    uint numberOfLayersBetweenLevels;

    SPtr<BoundingBox> subDomainBox;

public:

    VIRTUALFLUIDS_GPU_EXPORT void findCommunicationIndices( int direction, LbmOrGks lbmOrGks );
};

#endif


//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz, Martin Schoenherr
//=======================================================================================
#ifndef MULTIPLE_GRID_BUILDER_H
#define MULTIPLE_GRID_BUILDER_H

#include <vector>
#include <array>

#include "global.h"

#include "grid/GridBuilder/LevelGridBuilder.h"
#include "grid/GridFactory.h"
#include "grid/distributions/Distribution.h"
#include "grid/GridDimensions.h" 

class Object;
class BoundingBox;

class MultipleGridBuilder : public LevelGridBuilder
{
public:
    MultipleGridBuilder();

    virtual void addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ,
                                                    real delta);
    virtual void addCoarseGrid(const GridDimensions& gridDimensions);
    virtual void addGrid(SPtr<Object> gridShape);
    virtual void addGrid(SPtr<Object> gridShape, uint levelFine);

    virtual void addGeometry(SPtr<Object> gridShape);
    void addGeometry(SPtr<Object> solidObject, uint level);

    uint getNumberOfLevels() const;
    real getDelta(uint level) const;

    real getStartX(uint level) const;
    real getStartY(uint level) const;
    real getStartZ(uint level) const;

    real getEndX(uint level) const;
    real getEndY(uint level) const;
    real getEndZ(uint level) const;

    std::vector<SPtr<Grid> > getGrids() const;
    virtual void buildGrids(bool enableThinWalls = false);

    virtual void setNumberOfLayers(uint numberOfLayersFine, uint numberOfLayersBetweenLevels);

    void writeGridsToVtk(const std::string& path) const;

    virtual void setSubDomainBox(SPtr<BoundingBox> subDomainBox);

private:
    void addGridToList(SPtr<Grid> grid);
    real calculateDelta(uint level) const;
    bool coarseGridExists() const;
    bool isGridInCoarseGrid(SPtr<Grid> grid) const;

    void addFineGridToList(uint level, SPtr<Object> gridShape);
    void addIntermediateGridsToList(uint levelDifference, uint levelFine, uint nodesBetweenGrids, SPtr<Object>gridShape);
    void eraseGridsFromListIfInvalid(uint oldSize);
    void addGridToListIfValid(SPtr<Grid> grid);

    std::array<real, 6> getStaggeredCoordinates(SPtr<Object> gridShape, uint level, uint levelFine, bool &xOddStart, bool &yOddStart, bool &zOddStart) const;
    std::array<real, 6> getStaggeredCoordinates(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) const;
    std::array<real, 3> getOffset(real delta) const;
    std::vector<uint> getSpacingFactors(uint levelDifference) const;

    SPtr<Grid> makeGrid(SPtr<Object> gridShape, uint level, uint levelFine);
    SPtr<Grid> makeGrid(SPtr<Object> gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) const;

    static void emitNoCoarseGridExistsWarning();
    static void emitGridIsNotInCoarseGridWarning();

    SPtr<GridFactory> gridFactory;
    SPtr<Object> solidObject = nullptr;

    uint numberOfLayersFine;
    uint numberOfLayersBetweenLevels;

    SPtr<BoundingBox> subDomainBox;

public:
    virtual void findCommunicationIndices(int direction, bool doShift=false);
};

#endif


//! \}

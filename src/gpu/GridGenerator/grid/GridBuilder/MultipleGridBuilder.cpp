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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file MultipleGridBuilder.cpp
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz, Martin Sch�nherr
//=======================================================================================
#include "MultipleGridBuilder.h"

#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>

#include "utilities/math/Math.h"

#include "geometries/BoundingBox/BoundingBox.h"

#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/Grid.h"
#include "grid/GridFactory.h"

#include "io/GridVTKWriter/GridVTKWriter.h"
#include "io/STLReaderWriter/STLWriter.h"


MultipleGridBuilder::MultipleGridBuilder() : LevelGridBuilder(), numberOfLayersFine(12), numberOfLayersBetweenLevels(8), subDomainBox(nullptr)
{
    gridFactory = GridFactory::make();
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
}

void MultipleGridBuilder::addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta)
{
    boundaryConditions.push_back(std::make_shared<BoundaryConditions>());

    startX -= 0.5 * delta;
    startY -= 0.5 * delta;
    startZ -= 0.5 * delta;
    endX   += 0.5 * delta;
    endY   += 0.5 * delta;
    endZ   += 0.5 * delta;

    const auto grid = this->makeGrid(std::make_shared<Cuboid>(startX, startY, startZ, endX, endY, endZ), startX, startY, startZ, endX, endY, endZ, delta, 0);
    addGridToList(grid);
}

void MultipleGridBuilder::addCoarseGrid(const GridDimensions &gridDimensions)
{
    this->addCoarseGrid(gridDimensions.minX, gridDimensions.minY, gridDimensions.minZ, gridDimensions.maxX,
                        gridDimensions.maxY, gridDimensions.maxZ, gridDimensions.delta);
}

void MultipleGridBuilder::addGeometry(SPtr<Object> solidObject)
{
    this->solidObject = solidObject;

    for(auto bcs : boundaryConditions)
    {
        bcs->geometryBoundaryCondition = GeometryBoundaryCondition::make();
        bcs->geometryBoundaryCondition->side = SideFactory::make(SideType::GEOMETRY);
    }
}

void MultipleGridBuilder::addGeometry(SPtr<Object> solidObject, uint level)
{
    this->solidObject = solidObject;
    auto gridShape = solidObject->clone();
    gridShape->changeSizeByDelta(4.0);

    this->addGrid(gridShape, level);
}

void MultipleGridBuilder::addGrid(SPtr<Object> gridShape)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const auto grid = makeGrid(gridShape, getNumberOfLevels(), 0);

    addGridToListIfValid(grid);
}

void MultipleGridBuilder::addGrid(SPtr<Object> gridShape, uint levelFine)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    for( uint level = this->getNumberOfLevels(); level <= levelFine; level++ ){
        const auto grid = makeGrid(gridShape, level, levelFine);

        if(level != levelFine){
            grid->setInnerRegionFromFinerGrid(true);
            grid->setNumberOfLayers( this->numberOfLayersBetweenLevels );
        }
        else{
            grid->setNumberOfLayers( this->numberOfLayersFine );
        }

        grids.push_back(grid);
    }

    //////////////////////////////////////////////////////////////////////////
    // this old code by Soeren P. scaled the object
    // this did not work for concave geometries
    //////////////////////////////////////////////////////////////////////////

    //const uint nodesBetweenGrids = 12;
    //const uint levelDifference = levelFine - getNumberOfLevels();
    //const uint oldGridSize = this->getNumberOfLevels();

    //addIntermediateGridsToList(levelDifference, levelFine, nodesBetweenGrids, gridShape);
    //addFineGridToList(levelFine, gridShape->clone());


    //eraseGridsFromListIfInvalid(oldGridSize);
}

void MultipleGridBuilder::addFineGridToList(uint level, SPtr<Object> gridShape)
{
    const auto grid = makeGrid(gridShape, level, 0);
    grids.push_back(grid);
}

void MultipleGridBuilder::addIntermediateGridsToList(uint levelDifference, uint levelFine, uint nodesBetweenGrids, SPtr<Object> gridShape)
{
    if (levelDifference > 0)
    {
        auto spacings = getSpacingFactors(levelDifference);

        uint level = getNumberOfLevels();
        for (int i = levelDifference - 1; i >= 0; i--)
        {
            const real deltaToNewSize = nodesBetweenGrids * spacings[i] * calculateDelta(levelFine);
            SPtr<Object> gridShapeClone = gridShape->clone();
            gridShapeClone->changeSizeByDelta(deltaToNewSize);

            const auto grid = makeGrid(gridShapeClone, level++, 0);
            grids.push_back(grid);
        }
    }
}

std::vector<uint> MultipleGridBuilder::getSpacingFactors(uint levelDifference) const
{
    std::vector<uint> spacings(levelDifference);

    spacings[0] = uint(std::pow(2, 1));
    for (uint i = 1; i < levelDifference; i++)
        spacings[i] = spacings[i - 1] + uint(std::pow(2, i + 1));

    return spacings;
}

void MultipleGridBuilder::eraseGridsFromListIfInvalid(uint oldSize)
{
    if (!isGridInCoarseGrid(grids[oldSize]))
    {
        this->grids.erase(grids.begin() + oldSize, grids.end());
        emitGridIsNotInCoarseGridWarning();
    }
}

void MultipleGridBuilder::addGridToListIfValid(SPtr<Grid> grid)
{
    if (!isGridInCoarseGrid(grid))
        return emitGridIsNotInCoarseGridWarning();

    addGridToList(grid);
}

SPtr<Grid> MultipleGridBuilder::makeGrid(SPtr<Object> gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) const
{
    return gridFactory->makeGrid(gridShape, startX, startY, startZ, endX, endY, endZ, delta, level);
}

bool MultipleGridBuilder::coarseGridExists() const
{
    return !grids.empty();
}

SPtr<Grid> MultipleGridBuilder::makeGrid(SPtr<Object> gridShape, uint level, uint levelFine)
{
    boundaryConditions.push_back(std::make_shared<BoundaryConditions>());

    const real delta = calculateDelta(level);

    bool xOddStart = false, yOddStart = false, zOddStart = false;

    auto staggeredCoordinates = getStaggeredCoordinates(gridShape, level, levelFine, xOddStart, yOddStart, zOddStart);

    SPtr<Grid> newGrid = this->makeGrid(gridShape, staggeredCoordinates[0],
                                                   staggeredCoordinates[1],
                                                   staggeredCoordinates[2],
                                                   staggeredCoordinates[3],
                                                   staggeredCoordinates[4],
                                                   staggeredCoordinates[5], delta, level);

    newGrid->setOddStart( xOddStart, yOddStart, zOddStart );

    return newGrid;
}

real MultipleGridBuilder::calculateDelta(uint level) const
{
    real delta = this->getDelta(0);
    for (uint i = 0; i < level; i++)
        delta *= 0.5;
    return delta;
}

std::array<real, 6> MultipleGridBuilder::getStaggeredCoordinates(SPtr<Object> gridShape, uint level, uint levelFine, bool& xOddStart, bool& yOddStart, bool& zOddStart) const
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // This method computes the start and end coordinates with respect to the coarse grid
    // The following sketch visualizes this procedure for level 2:
    //
    //                          /----------------------- domain --------------------------------/ 
    //                          |                      /----------------- refinement region ------------------------/
    //                          |                      |                                        |                   |
    //                          |                      |                                        |                   |
    // Level 2:                 |                 2   2|  2   2   2   2   2   2   2   2   2   2 | 2  (2) (2) (2)    |
    // Level 1:             1   |   1       1       1  |    1       1       1       1       1   |   1      (1)      |
    // Level 0:         0       |       0              |0               0               0       |       0           |
    //                          |                      |                                        |                   |
    //                          |  out of refinement   |             inside refinement          |   out of domain   |   out of refinement
    //                          |                      |                                        |                   |
    //  Step  1) ------------------>x                  |                                        |   x<-------------------------------------------
    //  Step  2)                |    ------>x------>x------>x                                   |                   |
    //  Step  3)                |                      |  x<                                    |    >x             |
    //  Step  4)                |                 x<------                                      |      ------>x     |
    //  Step  5)                |                      |                                        | x<--x<--x<--      |
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real deltaCoarse = this->grids[level-1]->getDelta();
    const real delta = 0.5 * deltaCoarse;


    std::array<real, 6> staggeredCoordinates;

    real X1Minimum = gridShape->getX1Minimum();
    real X2Minimum = gridShape->getX2Minimum();
    real X3Minimum = gridShape->getX3Minimum();
    real X1Maximum = gridShape->getX1Maximum();
    real X2Maximum = gridShape->getX2Maximum();
    real X3Maximum = gridShape->getX3Maximum();

    if( levelFine >= level ){
        real overlap = 0.0;
        
        if(level != levelFine)
        {
            overlap += ( this->numberOfLayersBetweenLevels + 1 ) * delta;

            // use geometric series to account for overlap of higher levels:
            // overlap + overlap/2 + overlap/4 + ...
            overlap *= 2.0 * ( 1.0 - pow(0.5, levelFine - level ) );
        }

        // add overlap for finest level
        overlap += this->numberOfLayersFine * delta * pow(0.5, levelFine - level);

        X1Minimum -= overlap;
        X2Minimum -= overlap;
        X3Minimum -= overlap;
        X1Maximum += overlap;
        X2Maximum += overlap;
        X3Maximum += overlap;
    }

    // Step 1
    // go to boundary of coarse grid
    staggeredCoordinates[0] = this->grids[level-1]->getStartX();
    staggeredCoordinates[1] = this->grids[level-1]->getStartY();
    staggeredCoordinates[2] = this->grids[level-1]->getStartZ();
    staggeredCoordinates[3] = this->grids[level-1]->getEndX();
    staggeredCoordinates[4] = this->grids[level-1]->getEndY();
    staggeredCoordinates[5] = this->grids[level-1]->getEndZ();

    // Step 2
    // move to first coarse node in refinement region
    while (staggeredCoordinates[0] < X1Minimum) staggeredCoordinates[0] += deltaCoarse;
    while (staggeredCoordinates[1] < X2Minimum) staggeredCoordinates[1] += deltaCoarse;
    while (staggeredCoordinates[2] < X3Minimum) staggeredCoordinates[2] += deltaCoarse;
    while (staggeredCoordinates[3] > X1Maximum) staggeredCoordinates[3] -= deltaCoarse;
    while (staggeredCoordinates[4] > X2Maximum) staggeredCoordinates[4] -= deltaCoarse;
    while (staggeredCoordinates[5] > X3Maximum) staggeredCoordinates[5] -= deltaCoarse;

    // Step 3
    // make the grid staggered with one layer of stopper nodes on the outside
    staggeredCoordinates[0] -= 0.25 * deltaCoarse;
    staggeredCoordinates[1] -= 0.25 * deltaCoarse;
    staggeredCoordinates[2] -= 0.25 * deltaCoarse;
    staggeredCoordinates[3] += 0.25 * deltaCoarse;
    staggeredCoordinates[4] += 0.25 * deltaCoarse;
    staggeredCoordinates[5] += 0.25 * deltaCoarse;

    // Step 4
    // add two layers until the refinement region is completely inside the fine grid
    if (staggeredCoordinates[0] > X1Minimum) staggeredCoordinates[0] -= deltaCoarse;
    if (staggeredCoordinates[1] > X2Minimum) staggeredCoordinates[1] -= deltaCoarse;
    if (staggeredCoordinates[2] > X3Minimum) staggeredCoordinates[2] -= deltaCoarse;
    if (staggeredCoordinates[3] < X1Maximum) staggeredCoordinates[3] += deltaCoarse;
    if (staggeredCoordinates[4] < X2Maximum) staggeredCoordinates[4] += deltaCoarse;
    if (staggeredCoordinates[5] < X3Maximum) staggeredCoordinates[5] += deltaCoarse;

    // Step 5

    // if the mesh is bigger than the domain then it has an odd start
    if (staggeredCoordinates[0] < this->grids[level - 1]->getStartX()) xOddStart = true;
    if (staggeredCoordinates[1] < this->grids[level - 1]->getStartY()) yOddStart = true;
    if (staggeredCoordinates[2] < this->grids[level - 1]->getStartZ()) zOddStart = true;

    // if the refinement region is larger than the domain, then the start and end points are moved inwards again
    while (staggeredCoordinates[0] < this->grids[level - 1]->getStartX()) staggeredCoordinates[0] += delta;
    while (staggeredCoordinates[1] < this->grids[level - 1]->getStartY()) staggeredCoordinates[1] += delta;
    while (staggeredCoordinates[2] < this->grids[level - 1]->getStartZ()) staggeredCoordinates[2] += delta;
    while (staggeredCoordinates[3] > this->grids[level - 1]->getEndX()  ) staggeredCoordinates[3] -= delta;
    while (staggeredCoordinates[4] > this->grids[level - 1]->getEndY()  ) staggeredCoordinates[4] -= delta;
    while (staggeredCoordinates[5] > this->grids[level - 1]->getEndZ()  ) staggeredCoordinates[5] -= delta;

    return staggeredCoordinates;
}

std::array<real, 6> MultipleGridBuilder::getStaggeredCoordinates(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) const
{
    //previous version of Soeren P.    
    auto offset = getOffset(delta);

    const real startXStaggered = std::floor(startX) - offset[0];
    const real startYStaggered = std::floor(startY) - offset[1];
    const real startZStaggered = std::floor(startZ) - offset[2];

    const real endXStaggered = std::ceil(endX) + offset[0];
    const real endYStaggered = std::ceil(endY) + offset[1];
    const real endZStaggered = std::ceil(endZ) + offset[2];

    auto temporaryGrid = this->makeGrid(nullptr, startXStaggered, startYStaggered, startZStaggered, endXStaggered, endYStaggered, endZStaggered, delta, level);
    auto startStaggered = temporaryGrid->getMinimumOnNode(Vertex(startX, startY, startZ));
    auto endStaggered = temporaryGrid->getMaximumOnNode(Vertex(endX, endY, endZ));

    return std::array<real, 6>{startStaggered.x, startStaggered.y, startStaggered.z, endStaggered.x, endStaggered.y, endStaggered.z};
}

std::array<real, 3> MultipleGridBuilder::getOffset(real delta) const
{
    real offsetToNullStart = delta * 0.5;

    for (uint level = 1; level < getNumberOfLevels(); level++)
        offsetToNullStart += getDelta(level) * 0.5;

    const real offsetX = getStartX(0) - int(getStartX(0)) + offsetToNullStart;
    const real offsetY = getStartY(0) - int(getStartY(0)) + offsetToNullStart;
    const real offsetZ = getStartZ(0) - int(getStartZ(0)) + offsetToNullStart;
    return std::array<real, 3>{offsetX, offsetY, offsetZ};
}

bool MultipleGridBuilder::isGridInCoarseGrid(SPtr<Grid> grid) const
{
    return
        vf::Math::greaterEqual(grid->getStartX(), grids[0]->getStartX()) &&
        vf::Math::greaterEqual(grid->getStartY(), grids[0]->getStartY()) &&
        vf::Math::greaterEqual(grid->getStartZ(), grids[0]->getStartZ()) &&
        vf::Math::lessEqual(grid->getEndX(), grids[0]->getEndX()) &&
        vf::Math::lessEqual(grid->getEndY(), grids[0]->getEndY()) &&
        vf::Math::lessEqual(grid->getEndZ(), grids[0]->getEndZ());
}


uint MultipleGridBuilder::getNumberOfLevels() const
{
    return uint(grids.size());
}

real MultipleGridBuilder::getDelta(uint level) const
{
    if (grids.size() <= level)
        throw std::runtime_error("delta from invalid level was required.");
    return grids[level]->getDelta();
}

real MultipleGridBuilder::getStartX(uint level) const
{
    return grids[level]->getStartX();
}

real MultipleGridBuilder::getStartY(uint level) const
{
    return grids[level]->getStartY();
}

real MultipleGridBuilder::getStartZ(uint level) const
{
    return grids[level]->getStartZ();
}

real MultipleGridBuilder::getEndX(uint level) const
{
    return grids[level]->getEndX();
}

real MultipleGridBuilder::getEndY(uint level) const
{
    return grids[level]->getEndY();
}

real MultipleGridBuilder::getEndZ(uint level) const
{
    return grids[level]->getEndZ();
}

void MultipleGridBuilder::addGridToList(SPtr<Grid> grid)
{
    grids.push_back(grid);
}

std::vector<SPtr<Grid> > MultipleGridBuilder::getGrids() const
{
    return this->grids;
}

// this method initiates the actual grid generation
//
//  => before calling this MultipleGridBuilder::buildGrids(...), several options 
//     parameters and similar have to be set for the grid generation:
//      => MultipleGridBuilder::addCoarseGrid(...)                      => background grid
//      => MultipleGridBuilder::setNumberOfLayers(...)                  => additional layers outside refinement region
//      => MultipleGridBuilder::addGrid(..)                             => refinement region
//      => MultipleGridBuilder::setSubDomainBox(...)                    => sub domain for Multi GPU
//      => MultipleGridBuilder::addGeometry(...)                        => object for solid domain
//      => MultipleGridBuilder::setPeriodicBoundaryCondition(...)       => periodicity
//
//  => LBM boundary conditions are set after MultipleGridBuilder::buildGrids(...):
//      => LevelGridBuilder::setVelocityBoundaryCondition(...)
//      => LevelGridBuilder::setPressureBoundaryCondition(...)
//      => GeometryBoundaryCondition::setTangentialVelocityForPatch(...)
//      => VelocityBoundaryCondition::setVelocityProfile(...)
//
//  => Multi GPU connectivity is set after MultipleGridBuilder::buildGrids(...):
//      => MultipleGridBuilder::findCommunicationIndices(...)
//      => LevelGridBuilder::setCommunicationProcess(...)
//
void MultipleGridBuilder::buildGrids(bool enableThinWalls )
{
    //////////////////////////////////////////////////////////////////////////

    // initialize the grids:
    //      => allocate memory
    //      => set everything to INVALID_OUT_OF_GRID
    //      => find nodes inside the refinement region, either based on 
    //         an object or on the shape of the underlying fine grid
    //      => add overlap as specified by MultipleGridBuilder::setNumberOfLayers
    //      => fix odd cells, such that we obtain steps of size on the 
    //         boundary of the grid (for complete cells for interpolation)
    //      => fix refinement into the wall, i.e. make sure the refinement 
    //         goes straight into the walls
    //      => set stopper nodes to STOPPER_OUT_OF_GRID and STOPPER_OUT_OF_GRID_BOUNDARY
    //
    // the number of layers is passed to the Grid::initial(...) method as argument
    //
    // For further documentation of the grid initialization see Section 5.2.2 and
    // Figure 5.2 in the Dissertation of Stephan Lenz:
    // https://publikationsserver.tu-braunschweig.de/receive/dbbs_mods_00068716
    //
    for( int level = (int)grids.size()-1; level >= 0; level-- ) {

        VF_LOG_INFO("Start initializing level {}", level);

        // On the coarse grid every thing is Fluid (w.r.t. the refinement)
        // On the finest grid the Fluid region is defined by the Object
        // On the intermediate levels the Fluid region is defined by the fluid region of the finer level
        if(level == 0)
            grids[level]->inital( nullptr, 0 );
        else if(level == (int)grids.size() - 1)
            grids[level]->inital( nullptr, this->numberOfLayersFine );
        else
            grids[level]->inital( grids[level+1], this->numberOfLayersBetweenLevels );

        VF_LOG_INFO("Done initializing level {}", level);
    }

    //////////////////////////////////////////////////////////////////////////

    // set the solid region and find Qs
    // this is only done if a solid object was added to the GridBuilder
    //
    // For further documentation of the treatment of solid objects see Section 5.2.3
    // and figure 5.3 in the Dissertation of Stephan Lenz:
    // https://publikationsserver.tu-braunschweig.de/receive/dbbs_mods_00068716
    //
    if (solidObject)
    {
        VF_LOG_TRACE("Start with Q Computation");

        // Currently the solid object is only used on the finest grid,
        // because refinement into solid objects is not yet implemented.
        // If the solid object does not intersect the interfaces, it 
        // might be possible to have a solid object on more than the finest
        // level (See Bachelor thesis of Lennard Lux). To enable this, 
        // change the following two lines. This is not tested though!

        //for( uint level = 0; level < grids.size(); level++ )
        uint level = (uint)grids.size() - 1;
        {
            // the Grid::mesh(...) method distinguishes inside and outside regions
            // of the solid domain.:
            //      => set inner nodes to INVALID_SOLID
            //      => close needle sells
            //      => set one layer of STOPPER_SOLID nodes in the solid domain
            //      => set one layer of BC_SOLID nodes in the fluid domain
            grids[level]->mesh(solidObject.get());

            // if thin walls are enables additional BC_SOLID nodes are found by
            // Grid::findOs(...). To prevent the actual Q computation, 
            // Grid::enableFindSolidBoundaryNodes() and Grid::enableComputeQs()
            // set a flag that changes the behavior of Grid::findOs(...);
            // additionally some needle cells are closed in this process.
            if (enableThinWalls) {
                grids[level]->enableFindSolidBoundaryNodes();
                grids[level]->findQs(solidObject.get());
                grids[level]->closeNeedleCellsThinWall();
                grids[level]->enableComputeQs();
            }

            // compute the sub grid distances 
            // this works for STL and Sphere objects, but not yet for other primitives!
            grids[level]->findQs(solidObject.get());
        }

        VF_LOG_TRACE("Done with Q Computation");
    }

    //////////////////////////////////////////////////////////////////////////

    // find the interface interpolations cells
    // For further documentation see Section 5.2.4 and Figure 5.3 in the dissertation
    // of Stephan Lenz:
    // https://publikationsserver.tu-braunschweig.de/receive/dbbs_mods_00068716
    //
    for (size_t i = 0; i < grids.size() - 1; i++)
        grids[i]->findGridInterface(grids[i + 1]);

    //////////////////////////////////////////////////////////////////////////

    // set all points outside the sub domain to STOPPER_OUT_OF_GRID_BOUNDARY
    // and INVALID_OUT_OF_GRID
    if( this->subDomainBox )
        for (size_t i = 0; i < grids.size(); i++)
            grids[i]->limitToSubDomain( this->subDomainBox);

    //////////////////////////////////////////////////////////////////////////

    // shrinks the interface cell lists to the correct size
    for (size_t i = 0; i < grids.size() - 1; i++)
        grids[i]->repairGridInterfaceOnMultiGPU(grids[i + 1]);

    //////////////////////////////////////////////////////////////////////////

    // unrolls the matrix
    //      => computes the sparse indices
    //      => generates neighbor connectivity taking into account periodic boundaries
    //      => undates the interface connectivity to sparse indices (overwrites matrix indices!)
    for (size_t i = 0; i < grids.size() - 1; i++)
        grids[i]->findSparseIndices(grids[i + 1]);

    grids[grids.size() - 1]->findSparseIndices(nullptr);

    //////////////////////////////////////////////////////////////////////////
}

GRIDGENERATOR_EXPORT void MultipleGridBuilder::setNumberOfLayers(uint numberOfLayersFine, uint numberOfLayersBetweenLevels)
{
    this->numberOfLayersFine = numberOfLayersFine;
    this->numberOfLayersBetweenLevels = numberOfLayersBetweenLevels;
}

void MultipleGridBuilder::emitNoCoarseGridExistsWarning()
{
    VF_LOG_WARNING("No Coarse grid was added before. Actual Grid is not added, please create coarse grid before.");
}


void MultipleGridBuilder::emitGridIsNotInCoarseGridWarning()
{
    VF_LOG_WARNING("Grid lies not inside of coarse grid. Actual Grid is not added.");
}

void MultipleGridBuilder::findCommunicationIndices(int direction, bool doShift)
{
    VF_LOG_TRACE("Start findCommunicationIndices()");

    if( this->subDomainBox )
        for (size_t i = 0; i < grids.size(); i++)
            grids[i]->findCommunicationIndices(direction, this->subDomainBox, doShift);

    VF_LOG_TRACE("Done findCommunicationIndices()");
}

void MultipleGridBuilder::writeGridsToVtk(const std::string& path) const
{
    for(uint level = 0; level < grids.size(); level++)
    {
        std::stringstream ss;
        ss << path << level << ".vtk";

        GridVTKWriter::writeGridToVTKXML(grids[level], ss.str());

    }
}

GRIDGENERATOR_EXPORT void MultipleGridBuilder::setSubDomainBox(SPtr<BoundingBox> subDomainBox)
{
    this->subDomainBox = subDomainBox;
}


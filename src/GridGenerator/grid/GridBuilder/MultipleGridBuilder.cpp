#include "MultipleGridBuilder.h"

#include <sstream>
#include <vector>
#include <iostream>

#include "utilities/math/Math.h"
#include "../Grid.h"
#include "../GridFactory.h"

#include <VirtualFluidsBasics/utilities/logger/Logger.h>
#include "io/STLReaderWriter/STLWriter.h"
#include "io/GridVTKWriter/GridVTKWriter.h"
#include <grid/BoundaryConditions/BoundaryCondition.h>
#include <grid/BoundaryConditions/Side.h>

MultipleGridBuilder::MultipleGridBuilder(SPtr<GridFactory> gridFactory, Device device, const std::string &d3qxx) :
    LevelGridBuilder(device, d3qxx), gridFactory(gridFactory), solidObject(nullptr)
{

}

SPtr<MultipleGridBuilder> MultipleGridBuilder::makeShared(SPtr<GridFactory> gridFactory)
{
    return SPtr<MultipleGridBuilder>(new MultipleGridBuilder(gridFactory));
}

void MultipleGridBuilder::addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta)
{
    boundaryConditions.push_back(SPtr<BoundaryConditions>(new BoundaryConditions));

    const auto grid = this->makeGrid(new Cuboid(startX, startY, startZ, endX, endY, endZ), startX, startY, startZ, endX, endY, endZ, delta, 0);
    addGridToList(grid);
}

void MultipleGridBuilder::addGeometry(Object* solidObject)
{
    this->solidObject = solidObject;

    for(auto bcs : boundaryConditions)
    {
        bcs->geometryBoundaryCondition = GeometryBoundaryCondition::make();
        bcs->geometryBoundaryCondition->side = SideFactory::make(SideType::GEOMETRY);
    }
}

void MultipleGridBuilder::addGeometry(Object* solidObject, uint level)
{
    this->solidObject = solidObject;
    auto gridShape = solidObject->clone();
    gridShape->scale(4.0);

    //TriangularMesh* triangularMesh = dynamic_cast<TriangularMesh*>(gridShape);
    //STLWriter::writeSTL(triangularMesh->triangleVec, "D:/GRIDGENERATION/STL/bridge.stl");

    this->addGrid(gridShape, level);
}

void MultipleGridBuilder::addGrid(Object* gridShape)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const auto grid = makeGrid(gridShape, getNumberOfLevels(), 0);

    addGridToListIfValid(grid);
}

void MultipleGridBuilder::addGrid(Object* gridShape, uint levelFine)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    for( uint level = this->getNumberOfLevels(); level <= levelFine; level++ ){
        const auto grid = makeGrid(gridShape, level, levelFine);

        if(level != levelFine)
            grid->setInnerRegionFromFinerGrid(true);

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

void MultipleGridBuilder::addFineGridToList(uint level, Object* gridShape)
{
    const auto grid = makeGrid(gridShape, level, 0);
    grids.push_back(grid);
}

void MultipleGridBuilder::addIntermediateGridsToList(uint levelDifference, uint levelFine, uint nodesBetweenGrids, Object* gridShape)
{
    if (levelDifference > 0)
    {
        auto spacings = getSpacingFactors(levelDifference);

        // start = startFine - SUM(nodesBetweenGrids * 2^i * dxfine) 
        uint level = getNumberOfLevels();
        for (int i = levelDifference - 1; i >= 0; i--)
        {
            const real scalingFactor = nodesBetweenGrids * spacings[i] * calculateDelta(levelFine);
            Object* gridShapeClone = gridShape->clone();
            gridShapeClone->scale(scalingFactor);

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


SPtr<Grid> MultipleGridBuilder::makeGrid(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) const
{
    return gridFactory->makeGrid(gridShape, startX, startY, startZ, endX, endY, endZ, delta, level);
}

bool MultipleGridBuilder::coarseGridExists() const
{
    return !grids.empty();
}

SPtr<Grid> MultipleGridBuilder::makeGrid(Object* gridShape, uint level, uint levelFine)
{
    boundaryConditions.push_back(SPtr<BoundaryConditions>(new BoundaryConditions));

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

std::array<real, 6> MultipleGridBuilder::getStaggeredCoordinates(Object* gridShape, uint level, uint levelFine, bool& xOddStart, bool& yOddStart, bool& zOddStart) const
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //                          /----------------------- domain --------------------------------\ 
    //                          |                      /----------------- refinement region ------------------------\
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
        real overlap = 9 * delta;

        // use geometric series to account for overlapp of higher levels:
        // overlap + overlap/2 + overlap/4 + ...
        overlap *= 2.0 * ( 1.0 - pow(0.5, levelFine - level + 1 ) );

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
        throw std::exception("delta from invalid level was required.");
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

void MultipleGridBuilder::buildGrids(LbmOrGks lbmOrGks)
{
    //////////////////////////////////////////////////////////////////////////

    // orginal version with scaling the object (also use old version of MultipleGridBuilder::addGrid()
    //for (auto grid : grids)
    //    grid->inital();

     //new version with 
    for( int level = grids.size()-1; level >= 0; level-- ){

        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start initializing level " << level << "\n";

        // On the coarse grid every thing is Fluid (w.r.t. the refinement)
        // On the finest grid the Fluid region is defined by the Object
        // On the intermediate levels the Fluid region is defined by the fluid region of the finer level
        //if( level == 0 || level == grids.size()-1 )
        //    grids[level]->inital();
        //else
        //    grids[level]->inital( grids[level+1] );
        
        if( level == 0 || level == grids.size()-1 )
            grids[level]->inital( nullptr );
        else
            grids[level]->inital( grids[level+1] );

        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Done initializing level " << level << "\n";
    }

    //////////////////////////////////////////////////////////////////////////

    if (solidObject)
    {

        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start with Q Computation\n";

        grids[grids.size() - 1]->mesh(solidObject);
        grids[grids.size() - 1]->findQs(solidObject);

        //for (size_t i = 0; i < grids.size(); i++){
        //    grids[i]->mesh(solidObject);
        //    grids[i]->findQs(solidObject);
        //}

        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Done with Q Computation\n";
    }


    for (size_t i = 0; i < grids.size() - 1; i++)
        grids[i]->findGridInterface(grids[i + 1], lbmOrGks);

	if (lbmOrGks == LBM) {
		for (size_t i = 0; i < grids.size() - 1; i++)
			grids[i]->findSparseIndices(grids[i + 1]);

		grids[grids.size() - 1]->findSparseIndices(nullptr);
	}
}


void MultipleGridBuilder::emitNoCoarseGridExistsWarning()
{
    *logging::out << logging::Logger::WARNING << "No Coarse grid was added before. Actual Grid is not added, please create coarse grid before.\n";
}


void MultipleGridBuilder::emitGridIsNotInCoarseGridWarning()
{
    *logging::out << logging::Logger::WARNING << "Grid lies not inside of coarse grid. Actual Grid is not added.\n";
}

void MultipleGridBuilder::writeGridsToVtk(const std::string& path) const
{
    for(uint level = 0; level < grids.size(); level++)
    {
        std::stringstream ss;
        ss << path << level << ".vtk";

        GridVTKWriter::writeGridToVTKXML(grids[level], ss.str());
        //if( level != 0 )
        //    GridVTKWriter::writeInterpolationCellsToVTKXML(grids[level], grids[level-1], ss.str() + ".InterpolationCells");
        //else
        //    GridVTKWriter::writeInterpolationCellsToVTKXML(grids[level], nullptr       , ss.str() + ".InterpolationCells");
        //GridVTKWriter::writeSparseGridToVTK(grids[level], ss.str());
    }
}
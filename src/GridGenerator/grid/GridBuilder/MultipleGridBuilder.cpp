#include "MultipleGridBuilder.h"

#include <sstream>
#include <vector>

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

    const auto grid = makeGrid(gridShape, getNumberOfLevels());

    addGridToListIfValid(grid);
}

void MultipleGridBuilder::addGrid(Object* gridShape, uint levelFine)
{

    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const uint nodesBetweenGrids = 8;
    const uint levelDifference = levelFine - getNumberOfLevels();
    const uint oldGridSize = this->getNumberOfLevels();

    addIntermediateGridsToList(levelDifference, levelFine, nodesBetweenGrids, gridShape);
    addFineGridToList(levelFine, gridShape->clone());

    //eraseGridsFromListIfInvalid(oldGridSize);
}

void MultipleGridBuilder::addFineGridToList(uint level, Object* gridShape)
{
    const auto grid = makeGrid(gridShape, level);
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

            const auto grid = makeGrid(gridShapeClone, level++);
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

SPtr<Grid> MultipleGridBuilder::makeGrid(Object* gridShape, uint level)
{
    boundaryConditions.push_back(SPtr<BoundaryConditions>(new BoundaryConditions));

    const real delta = calculateDelta(level);

    //auto staggeredCoordinates = getStaggeredCoordinates(gridShape->getX1Minimum(), gridShape->getX2Minimum(), gridShape->getX3Minimum(), gridShape->getX1Maximum(), gridShape->getX2Maximum(), gridShape->getX3Maximum(), delta);
	auto staggeredCoordinates = getStaggeredCoordinates(gridShape, level);

	return this->makeGrid(gridShape, staggeredCoordinates[0], staggeredCoordinates[1], staggeredCoordinates[2], staggeredCoordinates[3], staggeredCoordinates[4], staggeredCoordinates[5], delta, level);
}

real MultipleGridBuilder::calculateDelta(uint level) const
{
    real delta = this->getDelta(0);
    for (uint i = 0; i < level; i++)
        delta *= 0.5;
    return delta;
}

std::array<real, 6> MultipleGridBuilder::getStaggeredCoordinates(Object* gridShape, uint level) const
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                     |
    //                                                     |
    // Level 2:                                           2   2   2   2   2   2   2   2   2   2   2   2
    // Level 1:                     1       1       1       1       1       1       1       1       1
    // Level 0:         0               0               0               0               0               0      
    //                                                     |
    //                               (coarse) outside      |      inside (fine)
    //                                                     |
    //  Step  1) ------>x                                  |
    //  Step  2)         -------------->x-------------->x-------------->x
    //  Step  3)                                                  x<----   
    //  Step  4)                                          x<------       
    //
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const real deltaCoarse = this->grids[0]->getDelta();

	std::array<real, 6> staggeredCoordinates;

	//Step 1
	staggeredCoordinates[0] = this->grids[0]->getStartX();
	staggeredCoordinates[1] = this->grids[0]->getStartY();
	staggeredCoordinates[2] = this->grids[0]->getStartZ();
	staggeredCoordinates[3] = this->grids[0]->getEndX();
	staggeredCoordinates[4] = this->grids[0]->getEndY();
	staggeredCoordinates[5] = this->grids[0]->getEndZ();

    //real X1Minimum = ( gridShape->getX1Minimum() > this->grids[0]->getStartX() ) ? ( gridShape->getX1Minimum() ) : ( this->grids[0]->getStartX() );
    //real X2Minimum = ( gridShape->getX2Minimum() > this->grids[0]->getStartY() ) ? ( gridShape->getX2Minimum() ) : ( this->grids[0]->getStartY() );
    //real X3Minimum = ( gridShape->getX3Minimum() > this->grids[0]->getStartZ() ) ? ( gridShape->getX3Minimum() ) : ( this->grids[0]->getStartZ() );
    //real X1Maximum = ( gridShape->getX1Maximum() < this->grids[0]->getEndX()   ) ? ( gridShape->getX1Maximum() ) : ( this->grids[0]->getEndX()   );
    //real X2Maximum = ( gridShape->getX2Maximum() < this->grids[0]->getEndY()   ) ? ( gridShape->getX2Maximum() ) : ( this->grids[0]->getEndY()   );
    //real X3Maximum = ( gridShape->getX3Maximum() < this->grids[0]->getEndZ()   ) ? ( gridShape->getX3Maximum() ) : ( this->grids[0]->getEndZ()   );

	//Step 2
	while (staggeredCoordinates[0] < gridShape->getX1Minimum()) staggeredCoordinates[0] += deltaCoarse;
	while (staggeredCoordinates[1] < gridShape->getX2Minimum()) staggeredCoordinates[1] += deltaCoarse;
	while (staggeredCoordinates[2] < gridShape->getX3Minimum()) staggeredCoordinates[2] += deltaCoarse;
	while (staggeredCoordinates[3] > gridShape->getX1Maximum()) staggeredCoordinates[3] -= deltaCoarse;
	while (staggeredCoordinates[4] > gridShape->getX2Maximum()) staggeredCoordinates[4] -= deltaCoarse;
	while (staggeredCoordinates[5] > gridShape->getX3Maximum()) staggeredCoordinates[5] -= deltaCoarse;

	//Step 3
	real offset = 0.0;
	for (int i = 0; i < level; i++)
		offset += 0.25 * pow(0.5, i) * deltaCoarse;

	staggeredCoordinates[0] -= offset;
	staggeredCoordinates[1] -= offset;
	staggeredCoordinates[2] -= offset;
	staggeredCoordinates[3] += offset;
	staggeredCoordinates[4] += offset;
	staggeredCoordinates[5] += offset;

	//Step 4
	while (staggeredCoordinates[0] > gridShape->getX1Minimum()) staggeredCoordinates[0] -= deltaCoarse * pow(0.5, level - 1);
	while (staggeredCoordinates[1] > gridShape->getX2Minimum()) staggeredCoordinates[1] -= deltaCoarse * pow(0.5, level - 1);
	while (staggeredCoordinates[2] > gridShape->getX3Minimum()) staggeredCoordinates[2] -= deltaCoarse * pow(0.5, level - 1);
	while (staggeredCoordinates[3] < gridShape->getX1Maximum()) staggeredCoordinates[3] += deltaCoarse * pow(0.5, level - 1);
	while (staggeredCoordinates[4] < gridShape->getX2Maximum()) staggeredCoordinates[4] += deltaCoarse * pow(0.5, level - 1);
	while (staggeredCoordinates[5] < gridShape->getX3Maximum()) staggeredCoordinates[5] += deltaCoarse * pow(0.5, level - 1);

    //Step 5
    while (staggeredCoordinates[0] < this->grids[level - 1]->getStartX()) staggeredCoordinates[0] += deltaCoarse * pow(0.5, level);
    while (staggeredCoordinates[1] < this->grids[level - 1]->getStartY()) staggeredCoordinates[1] += deltaCoarse * pow(0.5, level);
    while (staggeredCoordinates[2] < this->grids[level - 1]->getStartZ()) staggeredCoordinates[2] += deltaCoarse * pow(0.5, level);
    while (staggeredCoordinates[3] > this->grids[level - 1]->getEndX()  ) staggeredCoordinates[3] -= deltaCoarse * pow(0.5, level);
    while (staggeredCoordinates[4] > this->grids[level - 1]->getEndY()  ) staggeredCoordinates[4] -= deltaCoarse * pow(0.5, level);
    while (staggeredCoordinates[5] > this->grids[level - 1]->getEndZ()  ) staggeredCoordinates[5] -= deltaCoarse * pow(0.5, level);

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
    for (auto grid : grids)
        grid->inital();

    if (solidObject)
    {
        grids[grids.size() - 1]->mesh(solidObject);
        grids[grids.size() - 1]->findQs(solidObject);
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
        GridVTKWriter::writeSparseGridToVTK(grids[level], ss.str());
    }
}
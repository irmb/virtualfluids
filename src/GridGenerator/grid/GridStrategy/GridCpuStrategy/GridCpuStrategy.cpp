
#include "GridCpuStrategy.h"

#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <vector>

#include <GridGenerator/grid/distributions/Distribution.h>
#include <GridGenerator/grid/GridImp.h>

#include <GridGenerator/geometries/TriangularMesh/TriangularMesh.h>

#include <utilities/logger/Logger.h>
#include "grid/NodeValues.h"

#include "grid/GridInterface.h"

void GridCpuStrategy::allocateGridMemory(SPtr<GridImp> grid)
{
    grid->neighborIndexX = new int[grid->size];
    grid->neighborIndexY = new int[grid->size];
    grid->neighborIndexZ = new int[grid->size];

    grid->sparseIndices = new int[grid->size];

	grid->qIndices = new uint[grid->size];
	for (size_t i = 0; i < grid->size; i++) 
		grid->qIndices[i] = INVALID_INDEX;
	
 //   unsigned long distributionSize = grid->size * (grid->distribution.dir_end + 1);
 //   real sizeInMB = distributionSize * sizeof(real) / (1024.f*1024.f);

 //   //*logging::out << logging::Logger::LOW << "Allocating " << sizeInMB << " [MB] host memory for distributions.\n";

 //   grid->distribution.f = new real[distributionSize](); // automatic initialized with zeros
	//for (uint index = 0; index < distributionSize; index++) grid->distribution.f[index] = -1.0;
}

void GridCpuStrategy::allocateQs(SPtr<GridImp> grid)
{
	const uint numberOfQs = grid->getNumberOfSolidBoundaryNodes() * (grid->distribution.dir_end + 1);
	grid->qValues = new real[numberOfQs];
	for (size_t i = 0; i < numberOfQs; i++)
		grid->qValues[i] = -1.0;
}

void GridCpuStrategy::initalNodesToOutOfGrid(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
        grid->initalNodeToOutOfGrid(index);
}

void GridCpuStrategy::allocateFieldMemory(Field* field)
{
    field->field = new char[field->size];
}


void GridCpuStrategy::findInnerNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
        grid->findInnerNode(index);
}


void GridCpuStrategy::fixOddCells(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
        grid->fixOddCell(index);
}

}


void GridCpuStrategy::findStopperNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
        grid->findStopperNode(index);
}

void GridCpuStrategy::findEndOfGridStopperNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
	for (uint index = 0; index < grid->size; index++)
		grid->findEndOfGridStopperNode(index);
}

void GridCpuStrategy::findSolidStopperNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
	for (uint index = 0; index < grid->size; index++)
		grid->findSolidStopperNode(index);
}

void GridCpuStrategy::findBoundarySolidNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
	for (uint index = 0; index < grid->size; index++)
	{
		grid->findBoundarySolidNode(index);
	}
}

void GridCpuStrategy::mesh(SPtr<GridImp> grid, TriangularMesh &geom)
{
#pragma omp parallel for
    for (int i = 0; i < geom.size; i++)
        grid->mesh(geom.triangles[i]);
}

void GridCpuStrategy::findQs(SPtr<GridImp> grid, TriangularMesh &geom)
{
//#pragma omp parallel for
    for (int i = 0; i < geom.size; i++)
        grid->findQs(geom.triangles[i]);
}


void GridCpuStrategy::findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> fineGrid, LbmOrGks lbmOrGks)
{
    const auto coarseLevel = grid->getLevel();
    const auto fineLevel = fineGrid->getLevel();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "find interface level " << coarseLevel << " -> " << fineLevel;
    const uint oldGridSize = grid->getSparseSize();


    grid->gridInterface = new GridInterface();
    // TODO: this is stupid! concave refinements can easily have many more interface cells
    const uint sizeCF = (fineGrid->nx * fineGrid->ny + fineGrid->ny * fineGrid->nz + fineGrid->nx * fineGrid->nz);
    grid->gridInterface->cf.coarse = new uint[sizeCF];
    grid->gridInterface->cf.fine   = new uint[sizeCF];
    grid->gridInterface->cf.offset = new uint[sizeCF];
    grid->gridInterface->fc.coarse = new uint[sizeCF];
    grid->gridInterface->fc.fine   = new uint[sizeCF];
    grid->gridInterface->fc.offset = new uint[sizeCF];

    for (uint index = 0; index < grid->getSize(); index++)
        grid->findGridInterfaceCF(index, *fineGrid, lbmOrGks);

    for (uint index = 0; index < grid->getSize(); index++)
        grid->findGridInterfaceFC(index, *fineGrid);


    for (uint index = 0; index < grid->getSize(); index++)
        grid->findOverlapStopper(index, *fineGrid);

    const uint newGridSize = grid->getSparseSize();
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "  ... done. \n";
}

void GridCpuStrategy::findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Find sparse indices...";

    coarseGrid->updateSparseIndices();
    findForNeighborsNewIndices(coarseGrid);
    if (fineGrid)
    {
        fineGrid->updateSparseIndices();
        findForGridInterfaceNewIndices(coarseGrid, fineGrid);
    }

    const uint newGridSize = coarseGrid->getSparseSize();
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "... done. new size: " << newGridSize << ", delete nodes:" << coarseGrid->getSize() - newGridSize << "\n";
}


void GridCpuStrategy::findForNeighborsNewIndices(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->getSize(); index++)
        grid->setNeighborIndices(index);
}

void GridCpuStrategy::findForGridInterfaceNewIndices(SPtr<GridImp> grid, SPtr<GridImp> fineGrid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->getNumberOfNodesCF(); index++)
        grid->gridInterface->findForGridInterfaceSparseIndexCF(grid.get(), fineGrid.get(), index);

#pragma omp parallel for
    for (uint index = 0; index < grid->getNumberOfNodesFC(); index++)
        grid->gridInterface->findForGridInterfaceSparseIndexFC(grid.get(), fineGrid.get(), index);
}



void GridCpuStrategy::freeFieldMemory(Field* field)
{
    delete[] field->field;
}



void GridCpuStrategy::freeMemory(SPtr<GridImp> grid)
{
    //if(grid->gridInterface)
    //delete grid->gridInterface;

    delete[] grid->neighborIndexX;
    delete[] grid->neighborIndexY;
    delete[] grid->neighborIndexZ;
    delete[] grid->sparseIndices;
	delete[] grid->qIndices;
	delete[] grid->qValues;

    //delete[] grid->distribution.f;
}


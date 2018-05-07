
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
#include "io/GridVTKWriter/GridVTKWriter.h"
#include "numerics/geometry3d/GbTriFaceMesh3D.h"

void GridCpuStrategy::allocateGridMemory(SPtr<GridImp> grid)
{
    grid->neighborIndexX = new int[grid->size];
    grid->neighborIndexY = new int[grid->size];
    grid->neighborIndexZ = new int[grid->size];

    grid->sparseIndices = new int[grid->size];


    unsigned long distributionSize = grid->size * (grid->distribution.dir_end + 1);
    real sizeInMB = distributionSize * sizeof(real) / (1024.f*1024.f);

    //*logging::out << logging::Logger::LOW << "Allocating " << sizeInMB << " [MB] host memory for distributions.\n";

    grid->distribution.f = new real[distributionSize](); // automatic initialized with zeros
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


void GridCpuStrategy::removeOddBoundaryCellNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
        grid->removeOddBoundaryCellNode(index);
}


void GridCpuStrategy::findStopperNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
        grid->findStopperNode(index);
}

void GridCpuStrategy::mesh(SPtr<GridImp> grid, TriangularMesh &geom)
{
#pragma omp parallel for
    for (int i = 0; i < geom.size; i++)
        grid->mesh(geom.triangles[i]);
}

void GridCpuStrategy::findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> fineGrid)
{
    const auto coarseLevel = grid->getLevel(1.0);
    const auto fineLevel = fineGrid->getLevel(1.0);

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "find interface level " << coarseLevel << " -> " << fineLevel;
    const uint oldGridSize = grid->getSparseSize();


    grid->gridInterface = new GridInterface();
    const uint sizeCF = fineGrid->nx * fineGrid->ny + fineGrid->ny * fineGrid->nz + fineGrid->nx * fineGrid->nz;
    grid->gridInterface->cf.coarse = new uint[sizeCF];
    grid->gridInterface->cf.fine = new uint[sizeCF];
    grid->gridInterface->fc.coarse = new uint[sizeCF];
    grid->gridInterface->fc.fine = new uint[sizeCF];

    for (uint index = 0; index < grid->getSize(); index++)
        grid->findGridInterfaceCF(index, *fineGrid);

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

    delete[] grid->distribution.f;
}


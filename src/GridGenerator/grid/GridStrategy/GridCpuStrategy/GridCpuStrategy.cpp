#include "GridCpuStrategy.h"

#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <vector>

#include <GridGenerator/grid/distributions/Distribution.h>
#include <GridGenerator/grid/Grid.cuh>
      
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <utilities/logger/Logger.h>

void GridCpuStrategy::allocateGridMemory(SPtr<Grid> grid)
{
    grid->field = new char[grid->size];
        
    grid->neighborIndexX = new int[grid->size];
    grid->neighborIndexY = new int[grid->size];
    grid->neighborIndexZ = new int[grid->size];
        
    grid->matrixIndex = new unsigned int[grid->size];


    unsigned long distributionSize = grid->size * (grid->d.dir_end + 1);
    real sizeInMB = distributionSize * sizeof(real) / (1024.f*1024.f);

    *logging::out << logging::Logger::LOW << "Allocating " + SSTR(sizeInMB) + " [MB] host memory for distributions.\n";

    grid->d.f = new real[distributionSize](); // automatic initialized with zeros
}

void GridCpuStrategy::initalNodes(SPtr<Grid> grid)
{
#pragma omp parallel for
    for (unsigned int index = 0; index < grid->size; index++)
    {
        grid->setNeighborIndices(index);
        grid->matrixIndex[index] = index;
        grid->setFieldEntryToFluid(index);
    }
}

void GridCpuStrategy::mesh(SPtr<Grid> grid, Geometry &geom)
{
#pragma omp parallel for
    for (int i = 0; i < geom.size; i++)
        grid->meshTriangleExact(geom.triangles[i]);
}

void GridCpuStrategy::deleteSolidNodes(SPtr<Grid> grid)
{
    clock_t begin = clock();

    findInvalidNodes(grid);
    grid->removeInvalidNodes();
    findNeighborIndices(grid);

    clock_t end = clock();
    float time = (real)(real(end - begin) / CLOCKS_PER_SEC);
    *logging::out << logging::Logger::INTERMEDIATE << "time delete solid nodes: " + SSTR(time / 1000) + "sec\n";
}

void GridCpuStrategy::findInvalidNodes(SPtr<Grid> grid)
{
    bool foundInvalidNode = true;
    while (foundInvalidNode)
    {
        foundInvalidNode = false;
        for (int index = 0; index < grid->size; index++)
            grid->setInvalidNode(index, foundInvalidNode);
    }
}

void GridCpuStrategy::findNeighborIndices(SPtr<Grid> grid)
{
    for (int index = 0; index < grid->reducedSize; index++) 
        grid->findNeighborIndex(index);
}

void GridCpuStrategy::freeMemory(SPtr<Grid> grid)
{
    delete[] grid->field;

    delete[] grid->neighborIndexX;
    delete[] grid->neighborIndexY;
    delete[] grid->neighborIndexZ;
    delete[] grid->matrixIndex;

    delete[] grid->d.f;
}


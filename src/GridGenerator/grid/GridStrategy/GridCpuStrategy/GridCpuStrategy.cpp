#include "GridCpuStrategy.h"

#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <vector>

#include <GridGenerator/grid/distributions/Distribution.h>
#include <GridGenerator/grid/GridImp.cuh>
      
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <utilities/logger/Logger.h>
#include "grid/NodeValues.h"

#include "grid/GridInterface.cuh"

void GridCpuStrategy::allocateGridMemory(SPtr<GridImp> grid)
{
    grid->field = new char[grid->size];
        
    grid->neighborIndexX = new int[grid->size];
    grid->neighborIndexY = new int[grid->size];
    grid->neighborIndexZ = new int[grid->size];
        
    grid->matrixIndex = new int[grid->size];


    unsigned long distributionSize = grid->size * (grid->distribution.dir_end + 1);
    real sizeInMB = distributionSize * sizeof(real) / (1024.f*1024.f);

    *logging::out << logging::Logger::LOW << "Allocating " + SSTR(sizeInMB) + " [MB] host memory for distributions.\n";

    grid->distribution.f = new real[distributionSize](); // automatic initialized with zeros
}

void GridCpuStrategy::initalNodes(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
        grid->findInnerNode(index);

#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
        grid->findStopperNode(index);

    grid->removeInvalidNodes();
    findForNeighborsNewIndices(grid);
}

void GridCpuStrategy::mesh(SPtr<GridImp> grid, Geometry &geom)
{
#pragma omp parallel for
    for (int i = 0; i < geom.size; i++)
        grid->meshTriangle(geom.triangles[i]);
}

void GridCpuStrategy::findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> fineGrid)
{
    grid->gridInterface = new GridInterface();
    const uint sizeCF = fineGrid->nx * fineGrid->ny + fineGrid->ny * fineGrid->nz + fineGrid->nx * fineGrid->nz;
    grid->gridInterface->cf.coarse = new uint[sizeCF];
    grid->gridInterface->cf.fine = new uint[sizeCF];
    grid->gridInterface->fc.coarse = new uint[sizeCF];
    grid->gridInterface->fc.fine = new uint[sizeCF];
    grid->gridInterface->initalGridInterface(fineGrid.get());

    for (uint index = 0; index < grid->getSize(); index++)
        grid->findGridInterfaceCF(index, *fineGrid);

    for (uint index = 0; index < grid->getSize(); index++)
        grid->findGridInterfaceFC(index, *fineGrid);

    grid->removeInvalidNodes();
    findForNeighborsNewIndices(grid);
    findForGridInterfaceNewIndices(grid);
}



void GridCpuStrategy::findForNeighborsNewIndices(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->getSize(); index++)
        grid->findNeighborIndex(index);
}

void GridCpuStrategy::findForGridInterfaceNewIndices(SPtr<GridImp> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->getNumberOfNodesCF(); index++)
        grid->findForGridInterfaceNewIndexCF(index);

#pragma omp parallel for
    for (uint index = 0; index < grid->getNumberOfNodesFC(); index++)
        grid->findForGridInterfaceNewIndexFC(index);
}


void GridCpuStrategy::deleteSolidNodes(SPtr<GridImp> grid)
{
    clock_t begin = clock();

    findInvalidNodes(grid);
    grid->removeInvalidNodes();
    findForNeighborsNewIndices(grid);

    clock_t end = clock();
    real time = (real)(real(end - begin) / CLOCKS_PER_SEC);
    *logging::out << logging::Logger::INTERMEDIATE << "time delete solid nodes: " + SSTR(time / 1000) + "sec\n";
}

void GridCpuStrategy::findInvalidNodes(SPtr<GridImp> grid)
{
    bool foundInvalidNode = true;
    while (foundInvalidNode)
    {
        foundInvalidNode = false;
        for (uint index = 0; index < grid->size; index++)
            grid->setInvalidNode(index, foundInvalidNode);
    }
}



void GridCpuStrategy::freeMemory(SPtr<GridImp> grid)
{
    //if(grid->gridInterface)
        //delete grid->gridInterface;

    delete[] grid->field;

    delete[] grid->neighborIndexX;
    delete[] grid->neighborIndexY;
    delete[] grid->neighborIndexZ;
    delete[] grid->matrixIndex;

    delete[] grid->distribution.f;
}


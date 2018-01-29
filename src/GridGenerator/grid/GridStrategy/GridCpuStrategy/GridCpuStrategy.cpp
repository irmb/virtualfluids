#include "GridCpuStrategy.h"

#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <vector>

#include <GridGenerator/grid/distributions/Distribution.h>
#include <GridGenerator/grid/Grid.cuh>
      
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <utilities/logger/Logger.h>
#include "grid/NodeValues.h"

#include "grid/GridInterface.cuh"

void GridCpuStrategy::allocateGridMemory(SPtr<Grid> grid)
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

void GridCpuStrategy::initalNodes(SPtr<Grid> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
    {
        grid->matrixIndex[index] = index;
        if(grid->isEndOfGridStopper(index))
        {
            grid->setFieldEntryToSolid(index);
            grid->neighborIndexX[index] = -1;
            grid->neighborIndexY[index] = -1;
            grid->neighborIndexZ[index] = -1;
        } 
        else
        {
            grid->setFieldEntryToFluid(index);
            grid->setNeighborIndices(index);
        }
            
    }
}

void GridCpuStrategy::mesh(SPtr<Grid> grid, Geometry &geom)
{
#pragma omp parallel for
    for (int i = 0; i < geom.size; i++)
        grid->meshTriangle(geom.triangles[i]);
}

void GridCpuStrategy::removeOverlapNodes(SPtr<Grid> grid, SPtr<Grid> finerGrid)
{
    grid->gridInterface = new GridInterface(finerGrid.get());;

    createGridInterface(grid, finerGrid);
    grid->removeInvalidNodes();
    findForNeighborsNewIndices(grid);
    findForGridInterfaceNewIndices(grid);
}

void GridCpuStrategy::createGridInterface(SPtr<Grid> grid, SPtr<Grid> finerGrid)
{
    for (uint index = 0; index < grid->size; index++)
        grid->createGridInterface(index, *finerGrid.get());
}

void GridCpuStrategy::findForNeighborsNewIndices(SPtr<Grid> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->size; index++)
        grid->findNeighborIndex(index);
}

void GridCpuStrategy::findForGridInterfaceNewIndices(SPtr<Grid> grid)
{
#pragma omp parallel for
    for (uint index = 0; index < grid->gridInterface->cf.numberOfEntries; index++)
        grid->findForGridInterfaceNewIndexCF(index);

#pragma omp parallel for
    for (uint index = 0; index < grid->gridInterface->fc.numberOfEntries; index++)
        grid->findForGridInterfaceNewIndexFC(index);
}


void GridCpuStrategy::deleteSolidNodes(SPtr<Grid> grid)
{
    clock_t begin = clock();

    findInvalidNodes(grid);
    grid->removeInvalidNodes();
    findForNeighborsNewIndices(grid);

    clock_t end = clock();
    real time = (real)(real(end - begin) / CLOCKS_PER_SEC);
    *logging::out << logging::Logger::INTERMEDIATE << "time delete solid nodes: " + SSTR(time / 1000) + "sec\n";
}

void GridCpuStrategy::findInvalidNodes(SPtr<Grid> grid)
{
    bool foundInvalidNode = true;
    while (foundInvalidNode)
    {
        foundInvalidNode = false;
        for (uint index = 0; index < grid->size; index++)
            grid->setInvalidNode(index, foundInvalidNode);
    }
}



void GridCpuStrategy::freeMemory(SPtr<Grid> grid)
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


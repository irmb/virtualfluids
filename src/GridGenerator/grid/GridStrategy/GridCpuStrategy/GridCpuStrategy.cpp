
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

    *logging::out << logging::Logger::INTERMEDIATE
        << "Grid created: " << "from (" << grid->startX << ", " << grid->startY << ", " << grid->startZ << ") to (" << grid->endX << ", " << grid->endY << ", " << grid->endZ << ")\n"
        << "nodes: " << grid->nx << " x " << grid->ny << " x " << grid->nz << " = " << grid->size << "\n";
}

void GridCpuStrategy::findInnerNodes(SPtr<GridImp> grid, TriangularMesh* triangularMesh)
{
    switch (triangularMesh->getDiscretizationMethod())
    {
    case DiscretizationMethod::POINT_UNDER_TRIANGLE:
        pointUnderTriangleMethod(grid, triangularMesh);
        break;
    case DiscretizationMethod::RAYCASTING:
        rayCastingMethod(grid, triangularMesh);
        break;
    case DiscretizationMethod::POINT_IN_OBJECT:
        pointInObjectMethod(grid, triangularMesh);
    }

    removeOddBoundaryCellNodes(grid);
}

void GridCpuStrategy::pointUnderTriangleMethod(SPtr<GridImp> grid, TriangularMesh* triangularMesh)
{
#pragma omp parallel for
    for (int i = 0; i < triangularMesh->size; i++)
        grid->meshReverse(triangularMesh->triangles[i]);

    grid->findInsideNodes();

#pragma omp parallel for
    for (int i = 0; i < grid->size; i++)
        grid->setNegativeDirBorder_toFluid(i);

    *logging::out << logging::Logger::INTERMEDIATE
        << "Grid created: " << "from (" << grid->startX << ", " << grid->startY << ", " << grid->startZ << ") to (" << grid->endX << ", " << grid->endY << ", " << grid->endZ << ")\n"
        << "nodes: " << grid->nx << " x " << grid->ny << " x " << grid->nz << " = " << grid->size << "\n";
}

void GridCpuStrategy::pointInObjectMethod(SPtr<GridImp> grid, TriangularMesh* triangularMesh)
{
    auto triangles = triangularMesh->triangleVec;
    std::vector<GbTriFaceMesh3D::Vertex> *gbVertices = new std::vector<GbTriFaceMesh3D::Vertex>(triangles.size() * 3);
    std::vector<GbTriFaceMesh3D::TriFace> *gbTriangles = new std::vector<GbTriFaceMesh3D::TriFace>(triangles.size());
    for (int i = 0; i < triangles.size(); i++)
    {
        (*gbVertices)[i * 3] = GbTriFaceMesh3D::Vertex(triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
        (*gbVertices)[i * 3 + 1] = GbTriFaceMesh3D::Vertex(triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
        (*gbVertices)[i * 3 + 2] = GbTriFaceMesh3D::Vertex(triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);

        (*gbTriangles)[i] = GbTriFaceMesh3D::TriFace(i * 3, i * 3 + 1, i * 3 + 2);
    }

    GbTriFaceMesh3D mesh("stl", gbVertices, gbTriangles);

    for (int i = 0; i < grid->getSize(); i++)
    {
        real x, y, z;
        grid->transIndexToCoords(i, x, y, z);

        if (mesh.isPointInGbObject3D(x, y, z))
            grid->field.setFieldEntryToFluid(i);
    }
}

void GridCpuStrategy::rayCastingMethod(SPtr<GridImp> grid, TriangularMesh* triangularMesh)
{
    auto triangles = triangularMesh->triangleVec;
    std::vector<GbTriFaceMesh3D::Vertex> *gbVertices = new std::vector<GbTriFaceMesh3D::Vertex>(triangles.size() * 3);
    std::vector<GbTriFaceMesh3D::TriFace> *gbTriangles = new std::vector<GbTriFaceMesh3D::TriFace>(triangles.size());
    for (int i = 0; i < triangles.size(); i++)
    {
        (*gbVertices)[i * 3] = GbTriFaceMesh3D::Vertex(triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
        (*gbVertices)[i * 3 + 1] = GbTriFaceMesh3D::Vertex(triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
        (*gbVertices)[i * 3 + 2] = GbTriFaceMesh3D::Vertex(triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);
        
        (*gbTriangles)[i] = GbTriFaceMesh3D::TriFace(i * 3, i * 3 + 1, i * 3 + 2);
    }

    GbTriFaceMesh3D mesh("stl", gbVertices, gbTriangles);

    const real minXExact = triangularMesh->minmax.minX;
    const real minYExact = triangularMesh->minmax.minY;
    const real minZExact = triangularMesh->minmax.minZ;

    const real maxXExact = triangularMesh->minmax.maxX;
    const real maxYExact = triangularMesh->minmax.maxY;
    const real maxZExact = triangularMesh->minmax.maxZ;

    const auto min = grid->getMinimumOnNode(Vertex(minXExact, minYExact, minZExact));

    const real minX = min.x;
    const real minY = min.y;
    const real minZ = min.z;

    const auto max = grid->getMaximumOnNode(Vertex(maxXExact, maxYExact, maxZExact));

    const real maxX = max.x;
    const real maxY = max.y;
    const real maxZ = max.z;

    real x, y, z;

    for (int i = 0; i < grid->getSize(); i++)
    {
        grid->field.setFieldEntryToFluid(i);
    }


    int counter = 0;

    // Test line intersection
    for (z = minZ; z <= maxZ; z += grid->getDelta())
    {
        for (y = minY; y <= maxY; y += grid->getDelta())
        {
            for (x = minX; x <= maxX; x += grid->getDelta())
            {
                counter++;
                if (mesh.intersectLine((x - grid->getDelta()), y, z, x, y, z)) break;
                else grid->field.setFieldEntryToOutOfGrid(grid->transCoordToIndex(x, y, z));

            }
        }
    }

    // Test line intersection from opposite direction
    for (z = minZ; z <= maxZ; z += grid->getDelta())
    {
        for (y = minY; y <= maxY; y += grid->getDelta())
        {
            for (x = maxX; x >= minX; x -= grid->getDelta())
            {
                if (!grid->field.isOutOfGrid(grid->transCoordToIndex(x, y, z)))
                {
                    counter++;
                    if (mesh.intersectLine((x + grid->getDelta()), y, z, x, y, z)) break;
                    else grid->field.setFieldEntryToOutOfGrid(grid->transCoordToIndex(x, y, z));
                }
            }
        }
    }

    // Test line intersection
    for (z = minZ; z <= maxZ; z += grid->getDelta())
    {
        for (x = minX; x <= maxX; x += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                if (!grid->field.isOutOfGrid(grid->transCoordToIndex(x, y, z)))
                {
                    counter++;
                    if (mesh.intersectLine(x, (y - grid->getDelta()), z, x, y, z)) break;
                    else grid->field.setFieldEntryToOutOfGrid(grid->transCoordToIndex(x, y, z));
                }
            }
        }
    }

    // Test line intersection from opposite direction
    for (z = minZ; z <= maxZ; z += grid->getDelta())
    {
        for (x = minX; x <= maxX; x += grid->getDelta())
        {
            for (y = maxY; y >= minY; y -= grid->getDelta())
            {
                if (!grid->field.isOutOfGrid(grid->transCoordToIndex(x, y, z)))
                {
                    counter++;
                    if (mesh.intersectLine(x, (y + grid->getDelta()), z, x, y, z)) break;
                    else grid->field.setFieldEntryToOutOfGrid(grid->transCoordToIndex(x, y, z));
                }
            }
        }
    }

    // Test line intersection
    for (x = minX; x <= maxX; x += grid->getDelta())
    {
        for (y = minY; y <= maxY; y += grid->getDelta())
        {
            for (z = minZ; z <= maxZ; z += grid->getDelta())
            {
                if (!grid->field.isOutOfGrid(grid->transCoordToIndex(x, y, z)))
                {
                    counter++;
                    if (mesh.intersectLine(x, y, (z - grid->getDelta()), x, y, z)) break;
                    else grid->field.setFieldEntryToOutOfGrid(grid->transCoordToIndex(x, y, z));
                }
            }
        }
    }

    // Test line intersection from opposite direction
    for (x = minX; x <= maxX; x += grid->getDelta())
    {
        for (y = minY; y <= maxY; y += grid->getDelta())
        {
            for (z = maxZ; z >= minZ; z -= grid->getDelta())
            {
                if (!grid->field.isOutOfGrid(grid->transCoordToIndex(x, y, z)))
                {
                    counter++;
                    if (mesh.intersectLine(x, y, (z + grid->getDelta()), x, y, z)) break;
                    else grid->field.setFieldEntryToOutOfGrid(grid->transCoordToIndex(x, y, z));

                }
            }
        }
    }
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

    *logging::out << logging::Logger::INTERMEDIATE << "find interface level " << coarseLevel << " -> " << fineLevel;
    const uint oldGridSize = grid->getSparseSize();


    grid->gridInterface = new GridInterface();
    const uint sizeCF = fineGrid->nx * fineGrid->ny + fineGrid->ny * fineGrid->nz + fineGrid->nx * fineGrid->nz;
    grid->gridInterface->cf.coarse = new uint[sizeCF];
    grid->gridInterface->cf.fine = new uint[sizeCF];
    grid->gridInterface->fc.coarse = new uint[sizeCF];
    grid->gridInterface->fc.fine = new uint[sizeCF];

    for (uint index = 0; index < grid->getSize(); index++)
        grid->findGridInterfaceCF(index, *fineGrid);

    GridVTKWriter::writeSparseGridToVTK(fineGrid, "D:/GRIDGENERATION/CF");

    for (uint index = 0; index < grid->getSize(); index++)
        grid->findGridInterfaceFC(index, *fineGrid);


    for (uint index = 0; index < grid->getSize(); index++)
        grid->findOverlapStopper(index, *fineGrid);

    const uint newGridSize = grid->getSparseSize();
    *logging::out << logging::Logger::INTERMEDIATE << "  ... done. \n";
}

void GridCpuStrategy::findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid)
{
    *logging::out << logging::Logger::INTERMEDIATE << "Find sparse indices...";

    coarseGrid->updateSparseIndices();
    findForNeighborsNewIndices(coarseGrid);
    if (fineGrid)
    {
        fineGrid->updateSparseIndices();
        findForGridInterfaceNewIndices(coarseGrid, fineGrid);
    }

    const uint newGridSize = coarseGrid->getSparseSize();
    *logging::out << logging::Logger::INTERMEDIATE << "... done. new size: " << newGridSize << ", delete nodes:" << coarseGrid->getSize() - newGridSize << "\n";
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



void GridCpuStrategy::deleteSolidNodes(SPtr<GridImp> grid)
{
    clock_t begin = clock();

    grid->findInsideNodes();
    grid->updateSparseIndices();
    findForNeighborsNewIndices(grid);

    clock_t end = clock();
    real time = (real)(real(end - begin) / CLOCKS_PER_SEC);
    *logging::out << logging::Logger::INTERMEDIATE << "time delete solid nodes: " << time / 1000 << "sec\n";
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


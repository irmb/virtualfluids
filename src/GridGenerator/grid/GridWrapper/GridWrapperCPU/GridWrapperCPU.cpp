#include "GridWrapperCPU.h"

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <stack>
#include <vector>

#include <GridGenerator/grid/NodeValues.h>
#include <GridGenerator/grid/distributions/Distribution.h>
      
#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>
#include <GridGenerator/geometries/Triangle/Triangle.cuh>
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <utilities/logger/Logger.h>


GridWrapperCPU::GridWrapperCPU(BoundingBox<int> &box, std::string direction)
{
    int nx = box.maxX - box.minX;
    int ny = box.maxY - box.minY;
    int nz = box.maxZ - box.minZ;
    this->grid = Grid(NULL, box.minX, box.minY, box.minZ, nx, ny, nz, DistributionHelper::getDistribution(direction));
	this->allocDistribution();
	this->allocField();

	float time = this->initalUniformGrid3d();
    
	*logging::out << logging::Logger::INTERMEDIATE << "Time CPU initial field: " + SSTR(time / 1000) + "sec\n";
	*logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";
}

GridWrapperCPU::GridWrapperCPU(uint minX, uint minY, uint minZ, uint maxX, uint maxY, uint maxZ, std::string d3Qxx)
{
    int nx = maxX - minX;
    int ny = maxY - minY;
    int nz = maxZ - minZ;
    this->grid = Grid(NULL, minX, minY, minZ, nx, ny, nz, DistributionHelper::getDistribution(d3Qxx));
    this->allocDistribution();
    this->allocField();

    float time = this->initalUniformGrid3d();

    *logging::out << logging::Logger::INTERMEDIATE << "Time CPU initial field: " + SSTR(time / 1000) + "sec\n";
    *logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";
}

GridWrapperCPU::~GridWrapperCPU() 
{
    delete[] this->grid.field;

    delete[] this->grid.neighborIndexX;
    delete[] this->grid.neighborIndexY;
    delete[] this->grid.neighborIndexZ;
    delete[] this->grid.matrixIndex;

    delete[] this->grid.d.f;
}

real GridWrapperCPU::initalUniformGrid3d()
{
	time_t begin = clock();

	this->initialGridNodes();

	time_t end = clock();
	return (real)(real(end - begin) / CLOCKS_PER_SEC);
}

void GridWrapperCPU::initialGridNodes()
{
	#pragma omp parallel for
    for (unsigned int index = 0; index < grid.size; index++)
    {
        grid.setNeighborIndices(index);
        grid.matrixIndex[index] = index;
        grid.setFieldEntryToFluid(index);
    }
}

void GridWrapperCPU::meshGrid(Geometry &geom)
{
	clock_t begin = clock();
    runMeshing(geom);
	clock_t end = clock();
    float time = (real)(real(end - begin) / CLOCKS_PER_SEC);

	*logging::out << logging::Logger::INTERMEDIATE << "time grid generation: " + SSTR(time) + "s\n";
}

void GridWrapperCPU::deleteSolidNodes()
{
    clock_t begin = clock();

    findInvalidNodes();
    grid.removeInvalidNodes();
    findNeighborIndices();

    clock_t end = clock();
    float time = (real)(real(end - begin) / CLOCKS_PER_SEC);
    *logging::out << logging::Logger::INTERMEDIATE << "time delete solid nodes: " + SSTR(time / 1000) + "sec\n";
}

void GridWrapperCPU::allocDistribution()
{
	unsigned long distributionSize = grid.size * (grid.d.dir_end + 1);
	real sizeInMB = distributionSize * sizeof(real) / (1024.f*1024.f);

	*logging::out << logging::Logger::INTERMEDIATE << "Allocating " + SSTR(sizeInMB) + " [MB] host memory for distributions.\n";

	this->grid.d.f = new real[distributionSize](); // automatic initialized with zeros
}

void GridWrapperCPU::allocField()
{
	*logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";
	*logging::out << logging::Logger::INTERMEDIATE << "Initial field with fluid: \n";
	*logging::out << logging::Logger::INTERMEDIATE << "Field offset: " + SSTR(grid.startX) + ", " + SSTR(grid.startY) + ", " + SSTR(grid.startZ) + "\n";
	*logging::out << logging::Logger::INTERMEDIATE << "Field dimension: " + SSTR(grid.nx) + ", " + SSTR(grid.ny) + ", " + SSTR(grid.nz) + "\n";
	*logging::out << logging::Logger::INTERMEDIATE << "Number of Nodes: " + SSTR(grid.size) + "\n";

	this->grid.field = new char[this->grid.size];

    this->grid.neighborIndexX = new int[this->grid.size];
    this->grid.neighborIndexY = new int[this->grid.size];
    this->grid.neighborIndexZ = new int[this->grid.size];

    this->grid.matrixIndex = new unsigned int[this->grid.size];
}

void GridWrapperCPU::runMeshing(const Geometry &geom)
{   
    #pragma omp parallel for
	for (int i = 0; i < geom.size; i++)
		grid.meshTriangle(geom.triangles[i]);
}


void GridWrapperCPU::floodFill(const Vertex &startPoint)
{
    printf("start flood...\n");
	Vertex v(startPoint);
    std::stack<Vertex> floodFillStack;
    floodFillStack.push(v);

    while (!floodFillStack.empty())
	{
		v = floodFillStack.top();
        floodFillStack.pop();

        unsigned int x = (unsigned int)floor(v.x);
        unsigned int y = (unsigned int)floor(v.y);
        unsigned int z = (unsigned int)floor(v.z);

        int index = this->grid.transCoordToIndex(x, y, z);
        if (this->grid.isFluid(index))
		{
			this->grid.setFieldEntryToSolid(index);
            floodFillStack.push(Vertex(x, y, z + 1));
            floodFillStack.push(Vertex(x, y, z - 1));
            floodFillStack.push(Vertex(x, y + 1, z));
            floodFillStack.push(Vertex(x, y - 1, z));
            floodFillStack.push(Vertex(x + 1, y, z));
            floodFillStack.push(Vertex(x - 1, y, z));
        }
    }
    printf("...end flood\n");
}

void GridWrapperCPU::findInvalidNodes()
{
    bool foundInvalidNode = true;
    while (foundInvalidNode)
    {
        foundInvalidNode = false;
        for (int index = 0; index < grid.size; index++)
            grid.setInvalidNode(index, foundInvalidNode);
    }
}

void GridWrapperCPU::findNeighborIndices()
{
    for (int index = 0; index < grid.reducedSize; index++) 
        grid.findNeighborIndex(index);
}

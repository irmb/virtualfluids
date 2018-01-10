#include "GridWrapperGPU.h"

#include "time.h"
#include "GridGenerator/global.h"

#include <GridGenerator/utilities/cuda/CudaErrorCheck.cu>
#include <GridGenerator/utilities/Launchparameter/LaunchParameter.cuh>

#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <GridGenerator/grid/kernel/runGridKernelGPU.cuh>


#include <utilities/logger/Logger.h>


GridWrapperGPU::GridWrapperGPU(BoundingBox<int> &channel, std::string d3Qxx)
{
    initalField(channel, d3Qxx);
}

GridWrapperGPU::~GridWrapperGPU()
{
    delete[] this->grid.field;
    delete[] this->grid.d.f;
}

void GridWrapperGPU::meshGrid(Geometry &geom)
{
    *logging::out << logging::Logger::INTERMEDIATE << "start meshing on GPU...\n";
    allocAndCopyTrianglesToGPU(geom);

    /*---------------------------------------------------------------------------------*/
    float time = runKernelToMesh(LaunchParameter::make_1D1D_launchParameter(geom.size, 256), grid, geom);
    /*---------------------------------------------------------------------------------*/
    *logging::out << logging::Logger::INTERMEDIATE << "Time GPU build grid: " + SSTR(time / 1000) + "sec\n";
    *logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";

    freeTrianglesFromGPU(geom);
}

void GridWrapperGPU::deleteSolidNodes()
{
    float time1 = runKernelSetToInvalid(LaunchParameter::make_2D1D_launchParameter(grid.size, 512), grid);

    copyAndFreeFieldFromGPU();
    copyAndFreeMatrixIndicesFromGPU(grid.size);

    grid.removeInvalidNodes();

    allocAndCopyFieldToGPU();
    allocAndCopyMatrixIndicesToGPU();
    
    float time2 = runKernelFindIndices(LaunchParameter::make_2D1D_launchParameter(grid.reducedSize, 256), grid);
    *logging::out << logging::Logger::INTERMEDIATE << "time delete solid nodes: " + SSTR((time1+time2) / 1000) + "sec\n";
    
}

void GridWrapperGPU::floodFill(const Vertex &vec)
{

}

void GridWrapperGPU::copyDataFromGPU()
{
    copyAndFreeFieldFromGPU();
    copyAndFreeNeighborsToCPU();
    copyAndFreeMatrixIndicesFromGPU(grid.reducedSize);
    copyAndFreeDistributiondFromGPU();
}

void GridWrapperGPU::markNodesToDeleteOutsideOfGeometry()
{
    int numberOfEdgeNodes = grid.ny * grid.nz;

    /*---------------------------------------------------------------------------------*/
    float time = runKernelToMarkNodesToDeleteOutsideOfGeometry(LaunchParameter::make_1D1D_launchParameter(numberOfEdgeNodes, 256), grid);
    /*---------------------------------------------------------------------------------*/

    *logging::out << logging::Logger::INTERMEDIATE << "mark nodes to delete: " + SSTR(time / 1000) + "sec\n";
    *logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";
}

/*#################################################################################*/
/*--------------------------------private methods----------------------------------*/
/*---------------------------------------------------------------------------------*/
void GridWrapperGPU::initalField(BoundingBox<int> &box, std::string d3Qxx)
{
    int nx = box.maxX - box.minX;
    int ny = box.maxY - box.minY;
    int nz = box.maxZ - box.minZ;
    this->grid = Grid(NULL, box.minX, box.minY, box.minZ, nx, ny, nz, DistributionHelper::getDistribution(d3Qxx));
    this->allocDistribution();
    this->allocField();
    this->allocMatrixIndicesOnGPU();
    this->allocNeighborsIndices();

    /*---------------------------------------------------------------------------------*/
    float time = runKernelInitalUniformGrid3d(LaunchParameter::make_2D1D_launchParameter(grid.size, 256), grid);
    /*---------------------------------------------------------------------------------*/

    *logging::out << logging::Logger::INTERMEDIATE << "Time GPU initial field: " + SSTR(time / 1000) + "sec\n";
    *logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";
}

void GridWrapperGPU::allocField()
{
    *logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";
    *logging::out << logging::Logger::INTERMEDIATE << "Initial field with fluid: \n";
    *logging::out << logging::Logger::INTERMEDIATE << "Field offset: " + SSTR(grid.startX) + ", " + SSTR(grid.startY) + ", " + SSTR(grid.startZ) + "\n";
    *logging::out << logging::Logger::INTERMEDIATE << "Field dimension: " + SSTR(grid.nx) + ", " + SSTR(grid.ny) + ", " + SSTR(grid.nz) + "\n";
    *logging::out << logging::Logger::INTERMEDIATE << "Number of Nodes: " + SSTR(grid.size) + "\n";
    *logging::out << logging::Logger::INTERMEDIATE << "alloc on device for grid: " + SSTR(grid.size * sizeof(char) / 1000 / 1000) + " MB\n";

    char *field_d;
    CudaSafeCall(cudaMalloc(&field_d, grid.size * sizeof(char)));
    grid.field = field_d;
}

void GridWrapperGPU::allocAndCopyTrianglesToGPU(Geometry &geom)
{
    *logging::out << logging::Logger::INTERMEDIATE << "start copying triangles ...\n";
    clock_t begin = clock();
    int size_in_bytes_triangles = sizeof(Triangle)*geom.size;
    real sizeInMB = size_in_bytes_triangles / (1024.f*1024.f);

    *logging::out << logging::Logger::INTERMEDIATE << "Allocating " + SSTR(sizeInMB) + " [MB] device memory for triangles.\n\n";

    Triangle *triangles_d;
    CudaSafeCall(cudaMalloc(&triangles_d, size_in_bytes_triangles));
    CudaSafeCall(cudaMemcpy(triangles_d, geom.triangles, size_in_bytes_triangles, cudaMemcpyHostToDevice));
    geom.triangles = triangles_d;
    CudaCheckError();
    clock_t end = clock();
    real time = real(end - begin) / CLOCKS_PER_SEC;
    *logging::out << logging::Logger::INTERMEDIATE << "time copying triangles: " + SSTR(time) + "s\n";
    *logging::out << logging::Logger::INTERMEDIATE << "...copying triangles finish!\n\n";
}


void GridWrapperGPU::allocDistribution()
{
    CudaSafeCall(cudaMemcpyToSymbol(DIRECTIONS, grid.d.dirs, grid.d.dir_end * DIMENSION * sizeof(int)));
    CudaCheckError();

    unsigned long long distributionSize = grid.size * (grid.d.dir_end + 1);
    unsigned long long size_in_bytes = distributionSize * sizeof(real);
    real sizeInMB = size_in_bytes / (1024.f*1024.f);
    *logging::out << logging::Logger::INTERMEDIATE << "Allocating " + SSTR(sizeInMB) + " [MB] device memory for distributions.\n\n";

    CudaSafeCall(cudaMalloc(&grid.d.f, size_in_bytes));
    CudaCheckError();
}

void GridWrapperGPU::allocNeighborsIndices()
{
    int size_in_bytes_neighbors = grid.size * sizeof(int);
    int *neighborIndexX, *neighborIndexY, *neighborIndexZ;
    CudaSafeCall(cudaMalloc(&neighborIndexX, size_in_bytes_neighbors));
    CudaSafeCall(cudaMalloc(&neighborIndexY, size_in_bytes_neighbors));
    CudaSafeCall(cudaMalloc(&neighborIndexZ, size_in_bytes_neighbors));
    grid.neighborIndexX = neighborIndexX;
    grid.neighborIndexY = neighborIndexY;
    grid.neighborIndexZ = neighborIndexZ;
    CudaCheckError();
}

void GridWrapperGPU::allocMatrixIndicesOnGPU()
{
    int size_in_bytes_nodes_reduced = grid.size * sizeof(unsigned int);
    unsigned int* indices_reduced_d;
    CudaSafeCall(cudaMalloc(&indices_reduced_d, size_in_bytes_nodes_reduced));
    grid.matrixIndex = indices_reduced_d;
    CudaCheckError();
}


void GridWrapperGPU::allocAndCopyMatrixIndicesToGPU()
{
    int size_in_bytes_nodes_reduced = grid.reducedSize * sizeof(unsigned int);
    unsigned int* indices_reduced_d;
    CudaSafeCall(cudaMalloc(&indices_reduced_d, size_in_bytes_nodes_reduced));
    CudaSafeCall(cudaMemcpy(indices_reduced_d, grid.matrixIndex, size_in_bytes_nodes_reduced, cudaMemcpyHostToDevice));
    delete[] grid.matrixIndex;
    grid.matrixIndex = indices_reduced_d;
    CudaCheckError();
}

void GridWrapperGPU::allocAndCopyFieldToGPU()
{
    int size_in_bytes_grid = grid.size * sizeof(char);
    char* field_d;
    CudaSafeCall(cudaMalloc(&field_d, size_in_bytes_grid));
    CudaSafeCall(cudaMemcpy(field_d, grid.field, size_in_bytes_grid, cudaMemcpyHostToDevice));
    delete[] grid.field;
    grid.field = field_d;
    CudaCheckError();
}

void GridWrapperGPU::copyAndFreeFieldFromGPU()
{
    char *field_h = new char[grid.size];
    CudaSafeCall(cudaMemcpy(field_h, grid.field, grid.size * sizeof(char), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid.field));
    CudaCheckError();
    grid.field = field_h;
}

void GridWrapperGPU::copyAndFreeDistributiondFromGPU()
{
    unsigned long long distributionSize = grid.size * (grid.d.dir_end + 1);
    real *f_host = new real[distributionSize];
    CudaSafeCall(cudaMemcpy(f_host, grid.d.f, distributionSize * sizeof(real), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid.d.f));
    CudaCheckError();
    grid.d.f = f_host;
}


void GridWrapperGPU::freeTrianglesFromGPU(const Geometry &geom)
{
    CudaSafeCall(cudaFree(geom.triangles));
    CudaCheckError();
}

void GridWrapperGPU::copyAndFreeNeighborsToCPU()
{
    int size_in_bytes_neighbors = grid.size * sizeof(int);
    int *neighborIndexX_h, *neighborIndexY_h, *neighborIndexZ_h;
    neighborIndexX_h = new int[grid.size];
    neighborIndexY_h = new int[grid.size];
    neighborIndexZ_h = new int[grid.size];

    CudaSafeCall(cudaMemcpy(neighborIndexX_h, grid.neighborIndexX, size_in_bytes_neighbors, cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(neighborIndexY_h, grid.neighborIndexY, size_in_bytes_neighbors, cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(neighborIndexZ_h, grid.neighborIndexZ, size_in_bytes_neighbors, cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid.neighborIndexX));
    CudaSafeCall(cudaFree(grid.neighborIndexY));
    CudaSafeCall(cudaFree(grid.neighborIndexZ));

    grid.neighborIndexX = neighborIndexX_h;
    grid.neighborIndexY = neighborIndexY_h;
    grid.neighborIndexZ = neighborIndexZ_h;
    CudaCheckError();
}

void GridWrapperGPU::copyAndFreeMatrixIndicesFromGPU(int size)
{
    int size_in_bytes_nodes_reduced = size * sizeof(unsigned int);
    unsigned int *indices_reduced_h = new unsigned int[size];
    CudaSafeCall(cudaMemcpy(indices_reduced_h, grid.matrixIndex, size_in_bytes_nodes_reduced, cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid.matrixIndex));
    grid.matrixIndex = indices_reduced_h;
    CudaCheckError();
}




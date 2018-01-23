#include "GridGpuStrategy.h"

#include "time.h"
#include "GridGenerator/global.h"

#include <GridGenerator/utilities/cuda/CudaErrorCheck.cu>
#include <GridGenerator/utilities/Launchparameter/LaunchParameter.cuh>

#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <GridGenerator/grid/kernel/runGridKernelGPU.cuh>
#include <GridGenerator/grid/Grid.cuh>

#include <grid/distributions/Distribution.h>

#include <utilities/logger/Logger.h>
#include <helper_cuda.h>
#include "grid/GridInterface.cuh"


void GridGpuStrategy::allocateGridMemory(SPtr<Grid> grid)
{
    printCudaInformation(1);
    cudaSetDevice(1);

    this->allocDistribution(grid);
    this->allocField(grid);
    this->allocMatrixIndicesOnGPU(grid);
    this->allocNeighborsIndices(grid);
}

void GridGpuStrategy::initalNodes(SPtr<Grid> grid)
{
    float time = runKernelInitalUniformGrid3d(LaunchParameter::make_2D1D_launchParameter(grid->size, 256), *grid.get());
}

void GridGpuStrategy::mesh(SPtr<Grid> grid, Geometry &geom)
{
    *logging::out << logging::Logger::INTERMEDIATE << "start meshing on GPU...\n";
    allocAndCopyTrianglesToGPU(geom);

    /*---------------------------------------------------------------------------------*/
    float time = runKernelToMesh(LaunchParameter::make_1D1D_launchParameter(geom.size, 256), *grid.get(), geom);
    /*---------------------------------------------------------------------------------*/
    *logging::out << logging::Logger::INTERMEDIATE << "Time GPU build grid: " + SSTR(time / 1000) + "sec\n";
    *logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";

    freeTrianglesFromGPU(geom);

    copyDataFromGPU(grid);
}

void GridGpuStrategy::removeOverlapNodes(SPtr<Grid> grid, SPtr<Grid> finerGrid)
{
    float time1 = runKernelSetOverlapNodesToInvalid(LaunchParameter::make_2D1D_launchParameter(grid->size, 512), *grid.get(), *finerGrid.get());

    copyAndFreeFieldFromGPU(grid);
    copyAndFreeMatrixIndicesFromGPU(grid, grid->size);

    grid->removeInvalidNodes();

    allocAndCopyFieldToGPU(grid, grid->size);
    allocAndCopyMatrixIndicesToGPU(grid, grid->size);

    float time = runKernelFindIndices(LaunchParameter::make_2D1D_launchParameter(grid->size, 256), *grid.get());
    //*logging::out << logging::Logger::INTERMEDIATE << "time find indices: " + SSTR(time / 1000) + "sec\n";
}

void GridGpuStrategy::deleteSolidNodes(SPtr<Grid> grid)
{
    float time1 = runKernelSetToInvalid(LaunchParameter::make_2D1D_launchParameter(grid->size, 512), *grid.get());

    copyAndFreeFieldFromGPU(grid);
    copyAndFreeMatrixIndicesFromGPU(grid, grid->size);

    grid->removeInvalidNodes();

    allocAndCopyFieldToGPU(grid, grid->size);
    allocAndCopyMatrixIndicesToGPU(grid, grid->size);

    float time2 = runKernelFindIndices(LaunchParameter::make_2D1D_launchParameter(grid->size, 256), *grid.get());
    *logging::out << logging::Logger::INTERMEDIATE << "time delete solid nodes: " + SSTR((time1 + time2) / 1000) + "sec\n";
}


void GridGpuStrategy::freeMemory(SPtr<Grid> grid)
{
    delete grid->gridInterface;
    
    delete[] grid->field;

    delete[] grid->neighborIndexX;
    delete[] grid->neighborIndexY;
    delete[] grid->neighborIndexZ;
    delete[] grid->matrixIndex;

    delete[] grid->distribution.f;
}


void GridGpuStrategy::copyDataFromGPU(SPtr<Grid> grid)
{
    copyAndFreeFieldFromGPU(grid);
    copyAndFreeNeighborsToCPU(grid);
    copyAndFreeMatrixIndicesFromGPU(grid, grid->size);
    copyAndFreeDistributiondFromGPU(grid);
}
//
//void GridWrapperGPU::markNodesToDeleteOutsideOfGeometry()
//{
//    int numberOfEdgeNodes = grid.ny * grid.nz;
//
//    /*---------------------------------------------------------------------------------*/
//    float time = runKernelToMarkNodesToDeleteOutsideOfGeometry(LaunchParameter::make_1D1D_launchParameter(numberOfEdgeNodes, 256), grid);
//    /*---------------------------------------------------------------------------------*/
//
//    *logging::out << logging::Logger::INTERMEDIATE << "mark nodes to delete: " + SSTR(time / 1000) + "sec\n";
//    *logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";
//}

/*#################################################################################*/
/*--------------------------------private methods----------------------------------*/
/*---------------------------------------------------------------------------------*/
void GridGpuStrategy::allocField(SPtr<Grid> grid)
{
    *logging::out << logging::Logger::INTERMEDIATE << "alloc on device for grid field: " + SSTR(grid->size * sizeof(char) / 1000 / 1000) + " MB\n";

    char *field_d;
    CudaSafeCall(cudaMalloc(&field_d, grid->size * sizeof(char)));
    grid->field = field_d;
}

void GridGpuStrategy::allocDistribution(SPtr<Grid> grid)
{
    CudaSafeCall(cudaMemcpyToSymbol(DIRECTIONS, grid->distribution.dirs, grid->distribution.dir_end * DIMENSION * sizeof(int)));
    CudaCheckError();

    unsigned long long distributionSize = grid->size * (grid->distribution.dir_end + 1);
    unsigned long long size_in_bytes = distributionSize * sizeof(real);
    real sizeInMB = size_in_bytes / (1024.f*1024.f);
    *logging::out << logging::Logger::INTERMEDIATE << "Allocating " + SSTR(sizeInMB) + " [MB] device memory for distributions.\n\n";

    checkCudaErrors(cudaMalloc(&grid->distribution.f, size_in_bytes));
    CudaCheckError();
}

void GridGpuStrategy::allocNeighborsIndices(SPtr<Grid> grid)
{
    int size_in_bytes_neighbors = grid->size * sizeof(int);
    int *neighborIndexX, *neighborIndexY, *neighborIndexZ;
    CudaSafeCall(cudaMalloc(&neighborIndexX, size_in_bytes_neighbors));
    CudaSafeCall(cudaMalloc(&neighborIndexY, size_in_bytes_neighbors));
    CudaSafeCall(cudaMalloc(&neighborIndexZ, size_in_bytes_neighbors));
    grid->neighborIndexX = neighborIndexX;
    grid->neighborIndexY = neighborIndexY;
    grid->neighborIndexZ = neighborIndexZ;
    CudaCheckError();
}

void GridGpuStrategy::allocMatrixIndicesOnGPU(SPtr<Grid> grid)
{
    int size_in_bytes_nodes_reduced = grid->size * sizeof(int);
    int* indices_reduced_d;
    CudaSafeCall(cudaMalloc(&indices_reduced_d, size_in_bytes_nodes_reduced));
    grid->matrixIndex = indices_reduced_d;
    CudaCheckError();
}


void GridGpuStrategy::allocAndCopyTrianglesToGPU(Geometry &geom)
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

void GridGpuStrategy::freeTrianglesFromGPU(const Geometry &geom)
{
    CudaSafeCall(cudaFree(geom.triangles));
    CudaCheckError();
}

void GridGpuStrategy::allocAndCopyMatrixIndicesToGPU(SPtr<Grid> grid, const uint& size)
{
    int size_in_bytes_nodes_reduced = size * sizeof(int);
    int* indices_reduced_d;
    CudaSafeCall(cudaMalloc(&indices_reduced_d, size_in_bytes_nodes_reduced));
    CudaSafeCall(cudaMemcpy(indices_reduced_d, grid->matrixIndex, size_in_bytes_nodes_reduced, cudaMemcpyHostToDevice));
    delete[] grid->matrixIndex;
    grid->matrixIndex = indices_reduced_d;
    CudaCheckError();
}

void GridGpuStrategy::allocAndCopyFieldToGPU(SPtr<Grid> grid, const uint& size)
{
    int size_in_bytes_grid = size * sizeof(char);
    char* field_d;
    CudaSafeCall(cudaMalloc(&field_d, size_in_bytes_grid));
    CudaSafeCall(cudaMemcpy(field_d, grid->field, size_in_bytes_grid, cudaMemcpyHostToDevice));
    delete[] grid->field;
    grid->field = field_d;
    CudaCheckError();
}

void GridGpuStrategy::copyAndFreeFieldFromGPU(SPtr<Grid> grid)
{
    char *field_h = new char[grid->size];
    CudaSafeCall(cudaMemcpy(field_h, grid->field, grid->size * sizeof(char), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid->field));
    CudaCheckError();
    grid->field = field_h;
}

void GridGpuStrategy::copyAndFreeDistributiondFromGPU(SPtr<Grid> grid)
{
    unsigned long long distributionSize = grid->size * (grid->distribution.dir_end + 1);
    real *f_host = new real[distributionSize];
    CudaSafeCall(cudaMemcpy(f_host, grid->distribution.f, distributionSize * sizeof(real), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid->distribution.f));
    CudaCheckError();
    grid->distribution.f = f_host;
}



void GridGpuStrategy::copyAndFreeNeighborsToCPU(SPtr<Grid> grid)
{
    int size_in_bytes_neighbors = grid->size * sizeof(int);
    int *neighborIndexX_h, *neighborIndexY_h, *neighborIndexZ_h;
    neighborIndexX_h = new int[grid->size];
    neighborIndexY_h = new int[grid->size];
    neighborIndexZ_h = new int[grid->size];

    CudaSafeCall(cudaMemcpy(neighborIndexX_h, grid->neighborIndexX, size_in_bytes_neighbors, cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(neighborIndexY_h, grid->neighborIndexY, size_in_bytes_neighbors, cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(neighborIndexZ_h, grid->neighborIndexZ, size_in_bytes_neighbors, cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid->neighborIndexX));
    CudaSafeCall(cudaFree(grid->neighborIndexY));
    CudaSafeCall(cudaFree(grid->neighborIndexZ));

    grid->neighborIndexX = neighborIndexX_h;
    grid->neighborIndexY = neighborIndexY_h;
    grid->neighborIndexZ = neighborIndexZ_h;
    CudaCheckError();
}

void GridGpuStrategy::copyAndFreeMatrixIndicesFromGPU(SPtr<Grid> grid, int size)
{
    int size_in_bytes_nodes_reduced = size * sizeof(int);
    int *indices_reduced_h = new int[size];
    CudaSafeCall(cudaMemcpy(indices_reduced_h, grid->matrixIndex, size_in_bytes_nodes_reduced, cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid->matrixIndex));
    grid->matrixIndex = indices_reduced_h;
    CudaCheckError();
}




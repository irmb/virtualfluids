#include "GridGpuStrategy.h"

#include "time.h"
#include "GridGenerator/global.h"

#include <GridGenerator/utilities/cuda/CudaErrorCheck.cu>
#include <GridGenerator/utilities/Launchparameter/LaunchParameter.cuh>

#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <GridGenerator/grid/kernel/runGridKernelGPU.cuh>
#include <GridGenerator/grid/GridImp.h>

#include <grid/distributions/Distribution.h>

#include <utilities/logger/Logger.h>
#include <helper_cuda.h>
#include "grid/GridInterface.h"


void GridGpuStrategy::allocateGridMemory(SPtr<GridImp> grid)
{
    //printCudaInformation(0);
    cudaSetDevice(0);

    this->allocDistribution(grid);
    this->allocField(grid);
    this->allocMatrixIndicesOnGPU(grid);
    this->allocNeighborsIndices(grid);
}

void GridGpuStrategy::initalNodes(SPtr<GridImp> grid)
{
    float time = runKernelInitalUniformGrid3d(LaunchParameter::make_2D1D_launchParameter(grid->size, 256), *grid.get());
}

void GridGpuStrategy::mesh(SPtr<GridImp> grid, Geometry &geom)
{
    *logging::out << logging::Logger::INTERMEDIATE << "start meshing on GPU...\n";
    allocAndCopyTrianglesToGPU(geom);

    /*---------------------------------------------------------------------------------*/
    float time = runKernelToMesh(LaunchParameter::make_1D1D_launchParameter(geom.size, 256), *grid.get(), geom);
    /*---------------------------------------------------------------------------------*/
    *logging::out << logging::Logger::INTERMEDIATE << "Time GPU build grid: " + SSTR(time / 1000) + "sec\n";
    *logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";

    freeTrianglesFromGPU(geom);

}

void GridGpuStrategy::findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> fineGrid)
{
    copyAndFreeFieldFromGPU(grid->getField());
    copyAndFreeFieldFromGPU(fineGrid->getField());
    copyAndFreeMatrixIndicesFromGPU(grid, grid->size);
    copyAndFreeMatrixIndicesFromGPU(fineGrid, fineGrid->getSize());

    grid->gridInterface = new GridInterface();
    const uint sizeCF = fineGrid->nx * fineGrid->ny + fineGrid->ny * fineGrid->nz + fineGrid->nx * fineGrid->nz;
    grid->gridInterface->cf.coarse = new uint[sizeCF];
    grid->gridInterface->cf.fine = new uint[sizeCF];
    grid->gridInterface->fc.coarse = new uint[sizeCF];
    grid->gridInterface->fc.fine = new uint[sizeCF];

    //for (uint index = 0; index < grid->getSize(); index++)
    //    grid->findGridInterface(index, *fineGrid.get());

    uint *cfc;
    uint *cff;
    uint *fcc;
    uint *fcf;

    CudaSafeCall(cudaMalloc(&cfc, sizeCF * sizeof(uint)));
    CudaSafeCall(cudaMalloc(&cff, sizeCF * sizeof(uint)));
    CudaSafeCall(cudaMalloc(&fcc, sizeCF * sizeof(uint)));
    CudaSafeCall(cudaMalloc(&fcf, sizeCF * sizeof(uint)));

    CudaSafeCall(cudaMemcpy(cfc, grid->gridInterface->cf.coarse, grid->gridInterface->cf.numberOfEntries * sizeof(uint), cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(cff, grid->gridInterface->cf.fine, grid->gridInterface->cf.numberOfEntries * sizeof(uint), cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(fcc, grid->gridInterface->fc.coarse, grid->gridInterface->fc.numberOfEntries * sizeof(uint), cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(fcf, grid->gridInterface->fc.fine, grid->gridInterface->fc.numberOfEntries * sizeof(uint), cudaMemcpyHostToDevice));

    grid->gridInterface->cf.coarse = cfc;
    grid->gridInterface->cf.fine = cff;
    grid->gridInterface->fc.coarse = fcc;
    grid->gridInterface->fc.fine = fcf;

    GridInterface *gridInterface_d;
    CudaSafeCall(cudaMalloc(&gridInterface_d, sizeof(GridInterface)));
    CudaSafeCall(cudaMemcpy(gridInterface_d, grid->gridInterface, sizeof(GridInterface), cudaMemcpyHostToDevice));
    grid->gridInterface = gridInterface_d;
    CudaCheckError();


    grid->updateSparseIndices();

    allocAndCopyFieldToGPU(grid->getField());
    allocAndCopyFieldToGPU(fineGrid->getField());
    allocAndCopyMatrixIndicesToGPU(grid, grid->size);
    allocAndCopyMatrixIndicesToGPU(fineGrid, fineGrid->size);

    float time = runKernelToFindNeighborsNewIndices(LaunchParameter::make_2D1D_launchParameter(grid->size, 256), *grid.get());
    float time2 = runKernelToFindGridInterfaceNewIndices(LaunchParameter::make_2D1D_launchParameter(grid->size, 256), *grid.get());

    copyAndFreeGridInterfaceFromGPU(grid);


    //*logging::out << logging::Logger::INTERMEDIATE << "time find indices: " + SSTR(time / 1000) + "sec\n";
}

void GridGpuStrategy::deleteSolidNodes(SPtr<GridImp> grid)
{
    float time1 = runKernelSetToInvalid(LaunchParameter::make_2D1D_launchParameter(grid->size, 512), *grid.get());

    copyAndFreeFieldFromGPU(grid->getField());
    copyAndFreeMatrixIndicesFromGPU(grid, grid->size);

    grid->updateSparseIndices();

    allocAndCopyFieldToGPU(grid->getField());
    allocAndCopyMatrixIndicesToGPU(grid, grid->size);

    float time2 = runKernelFindIndices(LaunchParameter::make_2D1D_launchParameter(grid->size, 256), *grid.get());
    *logging::out << logging::Logger::INTERMEDIATE << "time delete solid nodes: " + SSTR((time1 + time2) / 1000) + "sec\n";
}



void GridGpuStrategy::freeMemory(SPtr<GridImp> grid)
{
    //delete[] grid->gridInterface->cf.coarse;
    //delete[] grid->gridInterface->cf.fine;

    //delete[] grid->gridInterface->fc.coarse;
    //delete[] grid->gridInterface->fc.fine;

    //delete grid->gridInterface;

    delete[] grid->neighborIndexX;
    delete[] grid->neighborIndexY;
    delete[] grid->neighborIndexZ;
    delete[] grid->sparseIndices;

    delete[] grid->distribution.f;
}



void GridGpuStrategy::copyDataFromGPU(SPtr<GridImp> grid)
{
    copyAndFreeFieldFromGPU(grid->getField());
    copyAndFreeNeighborsToCPU(grid);
    copyAndFreeMatrixIndicesFromGPU(grid, grid->size);
    copyAndFreeDistributiondFromGPU(grid);
}

void GridGpuStrategy::allocField(SPtr<GridImp> grid)
{
    //int size_in_bytes = grid->size * sizeof(char);
    //CudaSafeCall(cudaMalloc(grid->getField().field, size_in_bytes));
}


void GridGpuStrategy::allocateFieldMemory(Field* field)
{
    *logging::out << logging::Logger::INTERMEDIATE << "alloc on device for grid field: " + SSTR(field->size * sizeof(char) / 1000 / 1000) + " MB\n";

    char *field_d;
    CudaSafeCall(cudaMalloc(&field_d, field->size * sizeof(char)));
    field->field = field_d;
}

void GridGpuStrategy::freeFieldMemory(Field* field)
{
}

void GridGpuStrategy::findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid)
{
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

void GridGpuStrategy::allocDistribution(SPtr<GridImp> grid)
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

void GridGpuStrategy::allocNeighborsIndices(SPtr<GridImp> grid)
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

void GridGpuStrategy::allocMatrixIndicesOnGPU(SPtr<GridImp> grid)
{
    int size_in_bytes_nodes_reduced = grid->size * sizeof(int);
    int* indices_reduced_d;
    CudaSafeCall(cudaMalloc(&indices_reduced_d, size_in_bytes_nodes_reduced));
    grid->sparseIndices = indices_reduced_d;
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

void GridGpuStrategy::allocAndCopyMatrixIndicesToGPU(SPtr<GridImp> grid, const uint& size)
{
    int size_in_bytes_nodes_reduced = size * sizeof(int);
    int* indices_reduced_d;
    CudaSafeCall(cudaMalloc(&indices_reduced_d, size_in_bytes_nodes_reduced));
    CudaSafeCall(cudaMemcpy(indices_reduced_d, grid->sparseIndices, size_in_bytes_nodes_reduced, cudaMemcpyHostToDevice));
    delete[] grid->sparseIndices;
    grid->sparseIndices = indices_reduced_d;
    CudaCheckError();
}

void GridGpuStrategy::allocAndCopyFieldToGPU(Field& field)
{
    int size_in_bytes_grid = field.getSize() * sizeof(char);
    char* field_d;
    CudaSafeCall(cudaMalloc(&field_d, size_in_bytes_grid));
    CudaSafeCall(cudaMemcpy(field_d, field.field, size_in_bytes_grid, cudaMemcpyHostToDevice));
    delete[] field.field;
    field.field = field_d;
    CudaCheckError();
}

void GridGpuStrategy::copyAndFreeFieldFromGPU(Field& field)
{
    char *field_h = new char[field.size];
    CudaSafeCall(cudaMemcpy(field_h, field.field, field.size * sizeof(char), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(field.field));
    CudaCheckError();
    field.field = field_h;
}

void GridGpuStrategy::copyAndFreeDistributiondFromGPU(SPtr<GridImp> grid)
{
    unsigned long long distributionSize = grid->size * (grid->distribution.dir_end + 1);
    real *f_host = new real[distributionSize];
    CudaSafeCall(cudaMemcpy(f_host, grid->distribution.f, distributionSize * sizeof(real), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid->distribution.f));
    CudaCheckError();
    grid->distribution.f = f_host;
}



void GridGpuStrategy::copyAndFreeNeighborsToCPU(SPtr<GridImp> grid)
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

void GridGpuStrategy::copyAndFreeMatrixIndicesFromGPU(SPtr<GridImp> grid, int size)
{
    int size_in_bytes_nodes_reduced = size * sizeof(int);
    int *indices_reduced_h = new int[size];
    CudaSafeCall(cudaMemcpy(indices_reduced_h, grid->sparseIndices, size_in_bytes_nodes_reduced, cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(grid->sparseIndices));
    grid->sparseIndices = indices_reduced_h;
    CudaCheckError();
}


void GridGpuStrategy::copyAndFreeGridInterfaceFromGPU(SPtr<GridImp> grid)
{
    GridInterface *gridInterface = new GridInterface();
    CudaSafeCall(cudaMemcpy(gridInterface, grid->gridInterface, sizeof(GridInterface), cudaMemcpyDeviceToHost));

    uint *cfc = new uint[gridInterface->cf.numberOfEntries];
    uint *cff = new uint[gridInterface->cf.numberOfEntries];
    uint *fcc = new uint[gridInterface->fc.numberOfEntries];
    uint *fcf = new uint[gridInterface->fc.numberOfEntries];


    CudaSafeCall(cudaMemcpy(cfc, gridInterface->cf.coarse, gridInterface->cf.numberOfEntries * sizeof(uint), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(cff, gridInterface->cf.fine, gridInterface->cf.numberOfEntries * sizeof(uint), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(fcc, gridInterface->fc.coarse, gridInterface->fc.numberOfEntries * sizeof(uint), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(fcf, gridInterface->fc.fine, gridInterface->fc.numberOfEntries * sizeof(uint), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(gridInterface->cf.coarse));
    CudaSafeCall(cudaFree(gridInterface->cf.fine));
    CudaSafeCall(cudaFree(gridInterface->fc.coarse));
    CudaSafeCall(cudaFree(gridInterface->fc.fine));
    CudaSafeCall(cudaFree(grid->gridInterface));

    grid->gridInterface = gridInterface;
    grid->gridInterface->cf.coarse = cfc;
    grid->gridInterface->cf.fine = cff;
    grid->gridInterface->fc.coarse = fcc;
    grid->gridInterface->fc.fine = fcf;
    CudaCheckError();
}

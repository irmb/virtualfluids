#include "runGridKernelGPU.cuh"

#include <GridGenerator/utilities/cuda/cudaDefines.h>
#include <GridGenerator/utilities/cuda/cudaKernelCall.h>
#include <GridGenerator/utilities/Launchparameter/LaunchParameter.cuh>

#include <GridGenerator/grid/GridImp.h>
#include <GridGenerator/geometries/TriangularMesh/TriangularMesh.h>

GLOBAL void initalField(GridImp grid);
GLOBAL void runMeshing(GridImp grid, const TriangularMesh geom);

GLOBAL void findGridInterface(GridImp grid, GridImp finerGrid);
GLOBAL void runKernelTomarkNodesToDeleteOutsideOfGeometry(GridImp grid);
GLOBAL void markNodesToDeleteOutsideOfGeometry(GridImp grid);
GLOBAL void findNeighborIndicesKernel(GridImp grid);
GLOBAL void setOverlapNodesToInvalid(GridImp grid, GridImp finderGrid);

//////////////////////////////////////////////////////////////////////////

float runKernelInitalUniformGrid3d(const LaunchParameter& para, GridImp &grid)
{
    return runKernel(initalField, para, grid);
}

GLOBAL void initalField(GridImp grid)
{
    uint index = LaunchParameter::getGlobalIdx_2D_1D();
    if (index < grid.getSize())
        grid.findInnerNode(index);
}

//////////////////////////////////////////////////////////////////////////

float runKernelToMesh(const LaunchParameter& para, GridImp &grid, const TriangularMesh &geom)
{
    return runKernel(runMeshing, para, grid, geom);
}

GLOBAL void runMeshing(GridImp grid, const TriangularMesh geom)
{
    unsigned int i = LaunchParameter::getGlobalIdx_1D_1D();
    if (i < geom.size)
        grid.mesh(geom.triangles[i]);
}

//////////////////////////////////////////////////////////////////////////

float runKernelToMarkNodesToDeleteOutsideOfGeometry(const LaunchParameter& para, GridImp &grid)
{
    return runKernel(markNodesToDeleteOutsideOfGeometry, para, grid);
}

GLOBAL void markNodesToDeleteOutsideOfGeometry(GridImp grid)
{
    //int numberOfEdgeNodes = grid.ny * grid.nz;
    //unsigned int i = LaunchParameter::getGlobalIdx_1D_1D();

    //if (i < numberOfEdgeNodes)
    //{
    //    int x = 0;
    //    int y = i % grid.ny;
    //    int z = i / grid.ny;

    //    grid.setFieldEntryToSolid(grid.transCoordToIndex(x, y, z));

    //    while (grid.isFluid(i) && x < grid.nx - 1)
    //    {
    //        x += 1;
    //        grid.setFieldEntryToSolid(grid.transCoordToIndex(x, y, z));
    //    }
    //}

}


//////////////////////////////////////////////////////////////////////////

float runKernelSetOverlapNodesToInvalid(const LaunchParameter& para, GridImp &grid, GridImp &finerGrid)
{
    return runKernel(setOverlapNodesToInvalid, para, grid, finerGrid);
}

GLOBAL void setOverlapNodesToInvalid(GridImp grid, GridImp finerGrid)
{
    const uint index = LaunchParameter::getGlobalIdx_2D_1D();
    if (index < grid.getSize())
        grid.findGridInterfaceCF(index, finerGrid);
}



/*#################################################################################*/
/*---------------------------------find invalid nodes------------------------------*/
/*---------------------------------------------------------------------------------*/
GLOBAL void setInvalidNodes(GridImp grid, bool *foundInvalidNode);
/*---------------------------------------------------------------------------------*/

float runKernelSetToInvalid(const LaunchParameter& para, GridImp &grid)
{
    bool* foundInvalidNode_d;
    bool* foundInvalidNode_h = new bool();
    *foundInvalidNode_h = true;
    CudaSafeCall( cudaMalloc(&foundInvalidNode_d, sizeof(bool)));
    CudaCheckError();
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    while (*foundInvalidNode_h)
    {
        *foundInvalidNode_h = false;
        CudaSafeCall(cudaMemcpy(foundInvalidNode_d, foundInvalidNode_h, sizeof(bool), cudaMemcpyHostToDevice));
        setInvalidNodes << <para.blocks, para.threads >> >(grid, foundInvalidNode_d);
        cudaDeviceSynchronize();
        CudaCheckError();

        CudaSafeCall(cudaMemcpy(foundInvalidNode_h, foundInvalidNode_d, sizeof(bool), cudaMemcpyDeviceToHost));
    }

    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    CudaCheckError();
    cudaFree(foundInvalidNode_d);
    delete foundInvalidNode_h;

    return elapsedTime;
}


GLOBAL void setInvalidNodes(GridImp grid, bool *foundInvalidNode)
{
    uint index = LaunchParameter::getGlobalIdx_2D_1D();
    //if (index < grid.getSize())
        //grid.setInsideNode(index, *foundInvalidNode);
}


/*#################################################################################*/

float runKernelFindIndices(const LaunchParameter& para, GridImp &grid)
{
    return runKernel(findNeighborIndicesKernel, para, grid);
}

GLOBAL void findNeighborIndicesKernel(GridImp grid)
{
    uint index = LaunchParameter::getGlobalIdx_2D_1D();
    if (index < grid.getSize())
        grid.setNeighborIndices(index);
}


/*#################################################################################*/

float runKernelToFindGridInterface(const LaunchParameter& para, GridImp &grid, GridImp &finerGrid)
{
    return runKernel(findGridInterface, para, grid, finerGrid);
}

GLOBAL void findGridInterface(GridImp grid, GridImp finerGrid)
{
    uint index = LaunchParameter::getGlobalIdx_2D_1D();
    if (index < grid.getSize())
        grid.findGridInterfaceCF(index, finerGrid);
}
/*#################################################################################*/

GLOBAL void findNeighborsNewIndices(GridImp grid);
float runKernelToFindNeighborsNewIndices(const LaunchParameter& para, GridImp &grid)
{
    return runKernel(findNeighborsNewIndices, para, grid);
}

GLOBAL void findNeighborsNewIndices(GridImp grid)
{
    unsigned int index = LaunchParameter::getGlobalIdx_2D_1D();
    if (index < grid.getSize())
        grid.setNeighborIndices(index);
}
/*#################################################################################*/

GLOBAL void findGridInterfaceNewIndicesFC(GridImp grid);
GLOBAL void findGridInterfaceNewIndicesCF(GridImp grid);
float runKernelToFindGridInterfaceNewIndices(const LaunchParameter& para, GridImp &grid)
{
    runKernel(findGridInterfaceNewIndicesCF, para, grid);
    return runKernel(findGridInterfaceNewIndicesFC, para, grid);
}


GLOBAL void findGridInterfaceNewIndicesCF(GridImp grid)
{
    unsigned int index = LaunchParameter::getGlobalIdx_2D_1D();
    //if (index < grid.getNumberOfNodesCF())
        //grid.findForGridInterfaceSparseIndexCF(index);
}

GLOBAL void findGridInterfaceNewIndicesFC(GridImp grid)
{
    unsigned int index = LaunchParameter::getGlobalIdx_2D_1D();
    //if (index < grid.getNumberOfNodesFC())
    //    grid.findForGridInterfaceSparseIndexFCcoarse(index);
}

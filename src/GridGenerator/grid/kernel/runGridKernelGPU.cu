#include "runGridKernelGPU.cuh"

#include <GridGenerator/utilities/cuda/cudaDefines.h>
#include <GridGenerator/utilities/cuda/cudaKernelCall.h>
#include <GridGenerator/utilities/Launchparameter/LaunchParameter.cuh>

#include <GridGenerator/grid/Grid.cuh>
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

GLOBAL void initalField(Grid grid);
GLOBAL void runMeshing(Grid grid, const Geometry geom);
GLOBAL void runKernelTomarkNodesToDeleteOutsideOfGeometry(Grid grid);
GLOBAL void markNodesToDeleteOutsideOfGeometry(Grid grid);
GLOBAL void findNeighborIndicesKernel(Grid grid);
GLOBAL void setOverlapNodesToInvalid(Grid grid, Grid finderGrid);

//////////////////////////////////////////////////////////////////////////

float runKernelInitalUniformGrid3d(const LaunchParameter& para, Grid &grid)
{
	return runKernel(initalField, para, grid);
}

GLOBAL void initalField(Grid grid)
{
    unsigned int index = LaunchParameter::getGlobalIdx_2D_1D();
    if (index < grid.getSize())
    {
        grid.setNeighborIndices(index);
        grid.matrixIndex[index] = index;
        grid.setFieldEntryToFluid(index);
    }
}

//////////////////////////////////////////////////////////////////////////

float runKernelToMesh(const LaunchParameter& para, Grid &grid, const Geometry &geom)
{
	return runKernel(runMeshing, para, grid, geom);
}

GLOBAL void runMeshing(Grid grid, const Geometry geom)
{
	unsigned int i = LaunchParameter::getGlobalIdx_1D_1D();
	if (i < geom.size)
		grid.meshTriangle(geom.triangles[i]);
}

//////////////////////////////////////////////////////////////////////////

float runKernelToMarkNodesToDeleteOutsideOfGeometry(const LaunchParameter& para, Grid &grid)
{
	return runKernel(markNodesToDeleteOutsideOfGeometry, para, grid);
}

GLOBAL void markNodesToDeleteOutsideOfGeometry(Grid grid)
{
    int numberOfEdgeNodes = grid.ny * grid.nz;
    unsigned int i = LaunchParameter::getGlobalIdx_1D_1D();

    if (i < numberOfEdgeNodes)
    {
        int x = 0;
        int y = i % grid.ny;
        int z = i / grid.ny;

        grid.setFieldEntryToSolid(grid.transCoordToIndex(x, y, z));

        while (grid.isFluid(i) && x < grid.nx - 1)
        {
            x += 1;
            grid.setFieldEntryToSolid(grid.transCoordToIndex(x, y, z));
        }
    }

}


//////////////////////////////////////////////////////////////////////////

float runKernelSetOverlapNodesToInvalid(const LaunchParameter& para, Grid &grid, Grid &finerGrid)
{
    return runKernel(setOverlapNodesToInvalid, para, grid, finerGrid);
}

GLOBAL void setOverlapNodesToInvalid(Grid grid, Grid finerGrid)
{
    const uint index = LaunchParameter::getGlobalIdx_2D_1D();
    if (index < grid.getSize())
        grid.createGridInterface(index, finerGrid);
}



/*#################################################################################*/
/*---------------------------------find invalid nodes------------------------------*/
/*---------------------------------------------------------------------------------*/
GLOBAL void setInvalidNodes(Grid grid, bool *foundInvalidNode);
/*---------------------------------------------------------------------------------*/

float runKernelSetToInvalid(const LaunchParameter& para, Grid &grid)
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


GLOBAL void setInvalidNodes(Grid grid, bool *foundInvalidNode)
{
    unsigned int index = LaunchParameter::getGlobalIdx_2D_1D();
    if (index < grid.getSize())
        grid.setInvalidNode(index, *foundInvalidNode);
}


/*#################################################################################*/

float runKernelFindIndices(const LaunchParameter& para, Grid &grid)
{
    return runKernel(findNeighborIndicesKernel, para, grid);
}

GLOBAL void findNeighborIndicesKernel(Grid grid)
{
    unsigned int index = LaunchParameter::getGlobalIdx_2D_1D();
    if (index < grid.getSize())
        grid.findNeighborIndex(index);
}

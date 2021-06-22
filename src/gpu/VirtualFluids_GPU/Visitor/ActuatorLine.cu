#include "ActuatorLine.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

__global__ void interpolateVelocities(real* coordsX, real* coordsY, real* coordsZ, int* indices, int numberOfIndices);

void ActuatorLine::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{
    this->allocBladeCoords(cudaManager);
    this->initBladeCoords();
    this->initBoundingSphere(para, cudaManager);    
}

void ActuatorLine::visit(Parameter* para, int level, unsigned int t)
{
    if (level != this->level) return;

    this->copyBladeCoordsHtoD();

    unsigned int numberOfThreads = 128;
    int Grid = (this->numberOfIndices/ numberOfThreads)+1;
    int Grid1, Grid2;
    if (Grid>512)
    {
        Grid1 = 512;
        Grid2 = (Grid/Grid1)+1;
    }
    else
    {
        Grid1 = 1;
        Grid2 = Grid;
    }
    dim3 grid(Grid1, Grid2);
    dim3 threads(numberOfThreads, 1, 1 );

    interpolateVelocities<<< grid, threads >>>(
        para->getParD(this->level)->coordX_SP, 
        para->getParD(this->level)->coordY_SP, 
        para->getParD(this->level)->coordZ_SP,
        this->boundingSphereIndicesD,
        this->numberOfIndices);
    
}

void ActuatorLine::allocBladeCoords(CudaMemoryManager* cudaManager)
{
    checkCudaErrors( cudaMallocHost((void**) &this->bladeCoordsXH, this->mem_size_blades) );
    checkCudaErrors( cudaMallocHost((void**) &this->bladeCoordsYH, this->mem_size_blades) );
    checkCudaErrors( cudaMallocHost((void**) &this->bladeCoordsZH, this->mem_size_blades) );

    checkCudaErrors( cudaMalloc((void**) &this->bladeCoordsXD, this->mem_size_blades) );
    checkCudaErrors( cudaMalloc((void**) &this->bladeCoordsYD, this->mem_size_blades) );
    checkCudaErrors( cudaMalloc((void**) &this->bladeCoordsZD, this->mem_size_blades) );

    cudaManager->setMemsizeGPU(3.*this->mem_size_blades, false);
}

void ActuatorLine::initBladeCoords()
{   
    real dx = 0.5f*this->diameter/this->nBladeNodes;
    for(unsigned int node=0; node<this->nBladeNodes; node++)
    {
        this->bladeCoordsXH[node] = this->turbinePosX;
        this->bladeCoordsYH[node] = this->turbinePosY;
        this->bladeCoordsYH[node] = this->turbinePosZ+dx*node;

        this->bladeCoordsXH[node+this->nBladeNodes] = this->turbinePosX;
        this->bladeCoordsYH[node+this->nBladeNodes] = this->turbinePosY*0.5f*sqrt(3.0f)*dx*node;
        this->bladeCoordsYH[node+this->nBladeNodes] = this->turbinePosZ-0.5f*dx*node;

        this->bladeCoordsXH[node+2*this->nBladeNodes] = this->turbinePosX;
        this->bladeCoordsYH[node+2*this->nBladeNodes] = this->turbinePosY-0.5f*sqrt(3.0f)*dx*node;
        this->bladeCoordsYH[node+2*this->nBladeNodes] = this->turbinePosZ-0.5f*dx*node;
    }
}

void ActuatorLine::copyBladeCoordsHtoD()
{
    checkCudaErrors( cudaMemcpy(this->bladeCoordsXD, this->bladeCoordsXH, this->mem_size_blades, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(this->bladeCoordsYD, this->bladeCoordsYH, this->mem_size_blades, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(this->bladeCoordsZD, this->bladeCoordsZH, this->mem_size_blades, cudaMemcpyHostToDevice) );
}

void ActuatorLine::copyBladeCoordsDtoH()
{
    checkCudaErrors( cudaMemcpy(this->bladeCoordsXH, this->bladeCoordsXD, this->mem_size_blades, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(this->bladeCoordsYH, this->bladeCoordsYD, this->mem_size_blades, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(this->bladeCoordsZH, this->bladeCoordsZD, this->mem_size_blades, cudaMemcpyDeviceToHost) );
}

void ActuatorLine::freeBladeCoords()
{
    checkCudaErrors( cudaFree(this->bladeCoordsXD) );
    checkCudaErrors( cudaFree(this->bladeCoordsYD) );
    checkCudaErrors( cudaFree(this->bladeCoordsZD) );

    checkCudaErrors( cudaFreeHost(this->bladeCoordsXH) );
    checkCudaErrors( cudaFreeHost(this->bladeCoordsYH) );
    checkCudaErrors( cudaFreeHost(this->bladeCoordsZH) );
}
void ActuatorLine::rotateBlades(real angle)
{
    for(unsigned int node=0; node<this->nBladeNodes*this->nBlades; node++)
    {
        this->bladeCoordsYH[node] = this->bladeCoordsYH[node]*cos(angle)-this->bladeCoordsZH[node]*sin(angle);
        this->bladeCoordsZH[node] = this->bladeCoordsYH[node]*sin(angle)+this->bladeCoordsZH[node]*cos(angle);
    }
}

void ActuatorLine::initBoundingSphere(Parameter* para, CudaMemoryManager* cudaManager)
{
    // Actuator line exists only on 1 level
    std::vector<int> nodesInSphere;

    for (int j = 1; j <= para->getParH(this->level)->size_Mat_SP; j++)
    {
        const real coordX = para->getParH(this->level)->coordX_SP[j];
        const real coordY = para->getParH(this->level)->coordY_SP[j];
        const real coordZ = para->getParH(this->level)->coordZ_SP[j];
        const real dist = sqrt(pow(coordX-this->turbinePosX,2)+pow(coordY-this->turbinePosY,2)+pow(coordZ-this->turbinePosZ,2));
        if(dist < 0.6*this->diameter)
        {
            printf("indx in: %i \n", j);
            nodesInSphere.push_back(j);
        }
    }


    this->numberOfIndices = nodesInSphere.size();
    this->mem_size_boundingSphere = sizeof(int)*this->numberOfIndices;
    this->allocSphereIndices(cudaManager);
    this->boundingSphereIndicesH = nodesInSphere.data();
    this->copySphereIndices();
    
}

void ActuatorLine::allocSphereIndices(CudaMemoryManager* cudaManager)
{    
    printf("mem size sphere %i", (int)this->mem_size_boundingSphere);

    checkCudaErrors( cudaMallocHost((void**) &(this->boundingSphereIndicesH), this->mem_size_boundingSphere));
    checkCudaErrors( cudaMalloc((void**) &(this->boundingSphereIndicesD), this->mem_size_boundingSphere));
    cudaManager->setMemsizeGPU(this->mem_size_boundingSphere, false);
}

void ActuatorLine::copySphereIndices()
{
    checkCudaErrors( cudaMemcpy(this->boundingSphereIndicesD, this->boundingSphereIndicesH, this->mem_size_boundingSphere, cudaMemcpyHostToDevice) );
}

void ActuatorLine::freeSphereIndices()
{
    checkCudaErrors( cudaFreeHost(this->boundingSphereIndicesH) );
    checkCudaErrors( cudaFree(this->boundingSphereIndicesD) );
}

__global__ void interpolateVelocities(real* coordsX, real* coordsY, real* coordsZ, int* indices, int numberOfIndices)
{
    const uint x = threadIdx.x;  // Globaler x-Index 
    const uint y = blockIdx.x;   // Globaler y-Index 
    const uint z = blockIdx.y;   // Globaler z-Index 

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint index = nx*(ny*z + y) + x;

    if(index>=numberOfIndices) return;

    if(index==0||index+1==numberOfIndices)
    {
        printf("idx, x, y, z: %i, %f, %f, %f \n", indices[index], coordsX[index],coordsY[index],coordsZ[index]);
    }

}
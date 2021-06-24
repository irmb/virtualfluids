#include "ActuatorLine.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "Kernel/Utilities/CudaGrid.h"

__host__ __device__ uint find_nearest_cellTNE(uint index, 
                                              real* coordsX, real* coordsY, real* coordsZ, 
                                              real posX, real posY, real posZ, 
                                              uint* neighborsX, uint*neighborsY, uint* neighborsZ, uint* neighborsWSB);

__global__ void interpolateVelocities(real* globalCoordsX, real* globalCoordsY, real* globalCoordsZ, 
                                      uint* neighborsX, uint* neighborsY, uint* neighborsZ, 
                                      uint* neighborsWSB, 
                                      real* vx, real* vy, real* vz, 
                                      int* globalIndices, int numberOfIndices, 
                                      real* bladeCoordsX, real* bladeCoordsY, real* bladeCoordsZ, 
                                      real* bladeVelocitiesX, real* bladeVelocitiesY, real* bladeVelocitiesZ, 
                                      uint* bladeIndices, int numberOfNodes)
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint node = nx*(ny*z + y) + x;

    if(node>=numberOfNodes) return;

    real bladePosX = bladeCoordsX[node];
    real bladePosY = bladeCoordsY[node];
    real bladePosZ = bladeCoordsZ[node];

    uint old_index = bladeIndices[node];

    uint kTNE = find_nearest_cellTNE(old_index, globalCoordsX, globalCoordsY, globalCoordsZ, bladePosX, bladePosY, bladePosZ, neighborsX, neighborsY, neighborsZ, neighborsWSB);
    
    bladeIndices[node] = kTNE;
    
    real distX = bladePosX - globalCoordsX[kTNE];
    real distY = bladePosY - globalCoordsY[kTNE];
    real distZ = bladePosZ - globalCoordsZ[kTNE];

    uint kTNW = neighborsX[kTNE];
    uint kTSE = neighborsY[kTNE];
    uint kBNE = neighborsZ[kTNE];
    uint kTSW = neighborsY[kTNW];
    uint kBNW = neighborsZ[kTNW];
    uint kBSE = neighborsZ[kTSE];
    uint kBSW = neighborsZ[kTSW];

    //snaps to next, TODO interpolate
    real u_interpX = vx[kTNE];
    real u_interpY = vy[kTNE];
    real u_interpZ = vz[kTNE];

    bladeVelocitiesX[node] = u_interpX;
    bladeVelocitiesY[node] = u_interpY;
    bladeVelocitiesZ[node] = u_interpZ;

    // if(node==0 or node==90)
    // {
    //     printf("blade (%f, %f, %f), node TNE (%f, %f, %f), nodeBSW (%f, %f, %f)\n", bladePosX, bladePosY, bladePosZ, globalCoordsX[kTNE], globalCoordsY[kTNE], globalCoordsZ[kTNE], globalCoordsX[neighborsWSB[kTNE]], globalCoordsY[neighborsWSB[kTNE]], globalCoordsZ[neighborsWSB[kTNE]]);
    // }

}

__global__ void applyBodyForces(real* globalCoordsX, real* globalCoordsY, real* globalCoordsZ,
                           real* globalForcesX, real* globalForcesY, real* globalForcesZ, 
                           int* globalIndices, int numberOfIndices, 
                           real* bladeCoordsX, real* bladeCoordsY, real* bladeCoordsZ, 
                           real* bladeForcesX, real* bladeForcesY,real* bladeForcesZ,
                           uint* bladeIndices, int numberOfNodes)
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint index = nx*(ny*z + y) + x;

    if(index>=numberOfIndices) return;

    real posX = globalCoordsX[index];
    real posY = globalCoordsY[index];
    real posZ = globalCoordsZ[index];

    real forceX = 0.0f;
    real forceY = 0.0f;
    real forceZ = 0.0f;

    for( uint node=0; node<numberOfNodes; node++)
    {
        real distX = posX-bladeCoordsX[node];
        real distY = posY-bladeCoordsY[node];
        real distZ = posZ-bladeCoordsZ[node];
        real dist = sqrt(pow(distX,2.f)+pow(distY,2.f)+pow(distZ,2.f));

        forceX += bladeForcesX[node]*1/dist;
        forceY += bladeForcesY[node]*1/dist;
        forceZ += bladeForcesZ[node]*1/dist;
    }

    globalForcesX[index] = forceX;
    globalForcesY[index] = forceY;
    globalForcesZ[index] = forceZ;

}

__host__ __device__ uint find_nearest_cellTNE(uint index, 
                                              real* coordsX, real* coordsY, real* coordsZ, 
                                              real posX, real posY, real posZ, 
                                              uint* neighborsX, uint* neighborsY, uint* neighborsZ, uint* neighborsWSB){   
    uint new_index = index;

    while(coordsX[new_index] > posX && coordsY[new_index] > posY && coordsZ[new_index] > posZ ){ new_index = neighborsWSB[new_index];}

    while(coordsX[new_index] > posX && coordsY[new_index] > posY ){ new_index = neighborsZ[neighborsWSB[new_index]];}
    while(coordsX[new_index] > posX && coordsZ[new_index] > posZ ){ new_index = neighborsY[neighborsWSB[new_index]];}
    while(coordsY[new_index] > posY && coordsZ[new_index] > posZ ){ new_index = neighborsX[neighborsWSB[new_index]];}

    while(coordsX[new_index] > posX){ new_index = neighborsY[neighborsZ[neighborsWSB[new_index]]];}
    while(coordsY[new_index] > posY){ new_index = neighborsX[neighborsZ[neighborsWSB[new_index]]];}
    while(coordsZ[new_index] > posZ){ new_index = neighborsX[neighborsY[neighborsWSB[new_index]]];}

    while(coordsX[new_index] < posX){ new_index = neighborsX[new_index];}
    while(coordsY[new_index] < posY){ new_index = neighborsY[new_index];}
    while(coordsZ[new_index] < posZ){ new_index = neighborsZ[new_index];}

    return new_index;
}

void ActuatorLine::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{
    this->allocBladeCoords(cudaManager);
    this->allocBladeVelocities(cudaManager);
    this->allocBladeForces(cudaManager);
    this->allocBladeIndices(cudaManager);

    this->initBladeCoords();
    this->initBoundingSphere(para, cudaManager);   
    this->initBladeIndices(para);

    this->copyBladeIndicesHtoD();
}

void ActuatorLine::visit(Parameter* para, int level, unsigned int t)
{
    if (level != this->level) return;
    
    this->copyBladeCoordsHtoD();

    unsigned int numberOfThreads = 128;
    CudaGrid bladeGrid = CudaGrid(numberOfThreads, this->numberOfNodes);

    interpolateVelocities<<< bladeGrid.grid, bladeGrid.threads >>>(
        para->getParD(this->level)->coordX_SP, para->getParD(this->level)->coordY_SP, para->getParD(this->level)->coordZ_SP,        
        para->getParD(this->level)->neighborX_SP, para->getParD(this->level)->neighborY_SP, para->getParD(this->level)->neighborZ_SP, para->getParD(this->level)->neighborWSB_SP,
        para->getParD(this->level)->vx_SP, para->getParD(this->level)->vy_SP, para->getParD(this->level)->vz_SP,
        this->boundingSphereIndicesD, this->numberOfIndices,
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeVelocitiesXD, this->bladeVelocitiesYD, this->bladeVelocitiesZD,  
        this->bladeIndicesD, this->numberOfNodes);

    this->copyBladeVelocitiesDtoH();

    if(true)
    {
        real forceRatio = para->getVelocityRatio()*para->getDensityRatio();
        for( uint node=0; node<this->numberOfNodes; node++)
        {
            real u_rel = this->bladeVelocitiesXH[node]*para->getVelocityRatio();
            real v_rel = this->bladeVelocitiesYH[node]*para->getVelocityRatio()+this->omega*this->bladeRadii[node];

            real u_rel_sq = u_rel*u_rel+v_rel*v_rel;
            real phi = atan2(u_rel, v_rel);
            real Cl = 1.f;
            real Cd = 0.f;
            real c0 = 0.1f;
            real c = c0 * sqrt( 1- pow(2.f*(2*this->bladeRadii[node]-0.5*this->diameter)/(this->diameter), 2.f) );
            real Cn =   Cl*cos(phi)+Cd*sin(phi);
            real Ct =  -Cl*sin(phi)+Cd*cos(phi);

            real fRTZ_X = 0.5f*u_rel_sq*c*this->density*Cn;
            real fRTZ_Y = 0.5f*u_rel_sq*c*this->density*Ct;
            real fRTZ_Z = 0.0;
            real forceX = fRTZ_X;
            real forceY = fRTZ_Y;
            real forceZ = fRTZ_Z;
        
            this->bladeForcesXH[node] = forceX/forceRatio;
            this->bladeForcesYH[node] = forceY/forceRatio;
            this->bladeForcesZH[node] = forceZ/forceRatio;

            
        }
    }

    this->copyBladeForcesHtoD();

    CudaGrid sphereGrid = CudaGrid(numberOfThreads, this->numberOfIndices);

    applyBodyForces<<<sphereGrid.grid, sphereGrid.threads>>>(
        para->getParD(this->level)->coordX_SP, para->getParD(this->level)->coordY_SP, para->getParD(this->level)->coordZ_SP,        
        para->getParD(this->level)->forceX_SP, para->getParD(this->level)->forceY_SP, para->getParD(this->level)->forceZ_SP,        
        this->boundingSphereIndicesD, this->numberOfIndices,
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeForcesXD, this->bladeForcesYD, this->bladeForcesZD,  
        this->bladeIndicesD, this->numberOfNodes);

    real dazimuth = this->omega*para->getTimestepOfCoarseLevel()/pow(this->level,2.f);

    this->azimuth += dazimuth;
    this->rotateBlades(dazimuth);
    
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Blade coords

void ActuatorLine::allocBladeCoords(CudaMemoryManager* cudaManager)
{
    checkCudaErrors( cudaMallocHost((void**) &this->bladeCoordsXH, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMallocHost((void**) &this->bladeCoordsYH, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMallocHost((void**) &this->bladeCoordsZH, sizeof(real)*this->numberOfNodes) );

    checkCudaErrors( cudaMalloc((void**) &this->bladeCoordsXD, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMalloc((void**) &this->bladeCoordsYD, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMalloc((void**) &this->bladeCoordsZD, sizeof(real)*this->numberOfNodes) );

    cudaManager->setMemsizeGPU(3.*sizeof(real)*this->numberOfNodes, false);
}

void ActuatorLine::initBladeCoords()
{   
    //TODO bladeCoords for other nBlades than 3!
    real dx = 0.5f*this->diameter/this->nBladeNodes;        
    for( unsigned int blade=0; blade<this->nBlades; blade++)
    {
        real azimuth = (2*M_PI/this->nBlades*blade);
        for(unsigned int node=0; node<this->nBladeNodes; node++)
        {
            real r = dx*(node+1);
            this->bladeRadii[node] = r;
            this->bladeCoordsXH[node] = this->turbinePosX;
            this->bladeCoordsYH[node] = this->turbinePosY+r*sin(azimuth);
            this->bladeCoordsZH[node] = this->turbinePosZ+r*cos(azimuth);

        }
    }
}

void ActuatorLine::copyBladeCoordsHtoD()
{
    checkCudaErrors( cudaMemcpy(this->bladeCoordsXD, this->bladeCoordsXH, sizeof(real)*this->numberOfNodes, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(this->bladeCoordsYD, this->bladeCoordsYH, sizeof(real)*this->numberOfNodes, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(this->bladeCoordsZD, this->bladeCoordsZH, sizeof(real)*this->numberOfNodes, cudaMemcpyHostToDevice) );
}

void ActuatorLine::copyBladeCoordsDtoH()
{
    checkCudaErrors( cudaMemcpy(this->bladeCoordsXH, this->bladeCoordsXD, sizeof(real)*this->numberOfNodes, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(this->bladeCoordsYH, this->bladeCoordsYD, sizeof(real)*this->numberOfNodes, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(this->bladeCoordsZH, this->bladeCoordsZD, sizeof(real)*this->numberOfNodes, cudaMemcpyDeviceToHost) );
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Blade indices

void ActuatorLine::allocBladeIndices(CudaMemoryManager* cudaManager)
{
    checkCudaErrors( cudaMallocHost((void**) &this->bladeIndicesH, sizeof(uint)*this->numberOfNodes) );

    checkCudaErrors( cudaMalloc((void**) &this->bladeIndicesD, sizeof(uint)*this->numberOfNodes) );

    cudaManager->setMemsizeGPU(sizeof(uint)*this->numberOfNodes, false);
}

void ActuatorLine::initBladeIndices(Parameter* para)
{   

    for(unsigned int node=0; node<this->numberOfNodes; node++)
    {
        this->bladeIndicesH[node] = 1;
    }
}

void ActuatorLine::copyBladeIndicesHtoD()
{
    checkCudaErrors( cudaMemcpy(this->bladeIndicesD, this->bladeIndicesH, sizeof(uint)*this->numberOfNodes, cudaMemcpyHostToDevice) );
}

void ActuatorLine::freeBladeIndices()
{
    checkCudaErrors( cudaFree(this->bladeIndicesD) );

    checkCudaErrors( cudaFreeHost(this->bladeIndicesH) );
}


void ActuatorLine::rotateBlades(real angle)
{
    for(unsigned int node=0; node<this->nBladeNodes*this->nBlades; node++)
    {
        this->bladeCoordsYH[node] = this->bladeCoordsYH[node]*cos(angle)-this->bladeCoordsZH[node]*sin(angle);
        this->bladeCoordsZH[node] = this->bladeCoordsYH[node]*sin(angle)+this->bladeCoordsZH[node]*cos(angle);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Blade velocities

void ActuatorLine::allocBladeVelocities(CudaMemoryManager* cudaManager)
{
    checkCudaErrors( cudaMallocHost((void**) &this->bladeVelocitiesXH, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMallocHost((void**) &this->bladeVelocitiesYH, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMallocHost((void**) &this->bladeVelocitiesZH, sizeof(real)*this->numberOfNodes) );

    checkCudaErrors( cudaMalloc((void**) &this->bladeVelocitiesXD, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMalloc((void**) &this->bladeVelocitiesYD, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMalloc((void**) &this->bladeVelocitiesZD, sizeof(real)*this->numberOfNodes) );

    cudaManager->setMemsizeGPU(3.*sizeof(real)*this->numberOfNodes, false);
}

void ActuatorLine::copyBladeVelocitiesHtoD()
{
    checkCudaErrors( cudaMemcpy(this->bladeVelocitiesXD, this->bladeVelocitiesXH, sizeof(real)*this->numberOfNodes, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(this->bladeVelocitiesYD, this->bladeVelocitiesYH, sizeof(real)*this->numberOfNodes, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(this->bladeVelocitiesZD, this->bladeVelocitiesZH, sizeof(real)*this->numberOfNodes, cudaMemcpyHostToDevice) );
}

void ActuatorLine::copyBladeVelocitiesDtoH()
{
    checkCudaErrors( cudaMemcpy(this->bladeVelocitiesXH, this->bladeVelocitiesXD, sizeof(real)*this->numberOfNodes, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(this->bladeVelocitiesYH, this->bladeVelocitiesYD, sizeof(real)*this->numberOfNodes, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(this->bladeVelocitiesZH, this->bladeVelocitiesZD, sizeof(real)*this->numberOfNodes, cudaMemcpyDeviceToHost) );
}

void ActuatorLine::freeBladeVelocities()
{
    checkCudaErrors( cudaFree(this->bladeVelocitiesXD) );
    checkCudaErrors( cudaFree(this->bladeVelocitiesYD) );
    checkCudaErrors( cudaFree(this->bladeVelocitiesZD) );

    checkCudaErrors( cudaFreeHost(this->bladeVelocitiesXH) );
    checkCudaErrors( cudaFreeHost(this->bladeVelocitiesYH) );
    checkCudaErrors( cudaFreeHost(this->bladeVelocitiesZH) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Blade forces

void ActuatorLine::allocBladeForces(CudaMemoryManager* cudaManager)
{
    checkCudaErrors( cudaMallocHost((void**) &this->bladeForcesXH, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMallocHost((void**) &this->bladeForcesYH, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMallocHost((void**) &this->bladeForcesZH, sizeof(real)*this->numberOfNodes) );

    checkCudaErrors( cudaMalloc((void**) &this->bladeForcesXD, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMalloc((void**) &this->bladeForcesYD, sizeof(real)*this->numberOfNodes) );
    checkCudaErrors( cudaMalloc((void**) &this->bladeForcesZD, sizeof(real)*this->numberOfNodes) );

    cudaManager->setMemsizeGPU(3.*sizeof(real)*this->numberOfNodes, false);
}

void ActuatorLine::copyBladeForcesHtoD()
{
    checkCudaErrors( cudaMemcpy(this->bladeForcesXD, this->bladeForcesXH, sizeof(real)*this->numberOfNodes, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(this->bladeForcesYD, this->bladeForcesYH, sizeof(real)*this->numberOfNodes, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(this->bladeForcesZD, this->bladeForcesZH, sizeof(real)*this->numberOfNodes, cudaMemcpyHostToDevice) );
}

void ActuatorLine::copyBladeForcesDtoH()
{
    checkCudaErrors( cudaMemcpy(this->bladeForcesXH, this->bladeForcesXD, sizeof(real)*this->numberOfNodes, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(this->bladeForcesYH, this->bladeForcesYD, sizeof(real)*this->numberOfNodes, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(this->bladeForcesZH, this->bladeForcesZD, sizeof(real)*this->numberOfNodes, cudaMemcpyDeviceToHost) );
}

void ActuatorLine::freeBladeForces()
{
    checkCudaErrors( cudaFree(this->bladeForcesXD) );
    checkCudaErrors( cudaFree(this->bladeForcesYD) );
    checkCudaErrors( cudaFree(this->bladeForcesZD) );

    checkCudaErrors( cudaFreeHost(this->bladeForcesXH) );
    checkCudaErrors( cudaFreeHost(this->bladeForcesYH) );
    checkCudaErrors( cudaFreeHost(this->bladeForcesZH) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Sphere indices

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
        
        if(dist < 0.6*this->diameter) nodesInSphere.push_back(j);
    }

    this->numberOfIndices = nodesInSphere.size();
    this->allocSphereIndices(cudaManager);
    std::copy(nodesInSphere.begin(), nodesInSphere.end(), this->boundingSphereIndicesH);
    this->copySphereIndices();
    
}

void ActuatorLine::allocSphereIndices(CudaMemoryManager* cudaManager)
{    
    checkCudaErrors( cudaMallocHost((void**) &(this->boundingSphereIndicesH), sizeof(int)*this->numberOfIndices));
    checkCudaErrors( cudaMalloc((void**) &(this->boundingSphereIndicesD), sizeof(int)*this->numberOfIndices));
    cudaManager->setMemsizeGPU(sizeof(int)*this->numberOfIndices, false);
}

void ActuatorLine::copySphereIndices()
{
    checkCudaErrors( cudaMemcpy(this->boundingSphereIndicesD, this->boundingSphereIndicesH, sizeof(int)*this->numberOfIndices, cudaMemcpyHostToDevice) );
}

void ActuatorLine::freeSphereIndices()
{
    checkCudaErrors( cudaFreeHost(this->boundingSphereIndicesH) );
    checkCudaErrors( cudaFree(this->boundingSphereIndicesD) );
}

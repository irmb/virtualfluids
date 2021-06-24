#include "ActuatorLine.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "Kernel/Utilities/CudaGrid.h"
#include "lbm/constants/NumericConstants.h"
#include "Utilities/geometry_helper.h"

__host__ __device__ __inline__ real calc_gaussian3D(real posX, real posY, real posZ, real destX, real destY, real destZ, real epsilon)
{
    real distX = destX-posX;
    real distY = destY-posY;
    real distZ = destZ-posZ;
    real dist = sqrt(distX*distX+distY*distY+distZ*distZ);
    real oneOeps_sq = 1.f/(epsilon*epsilon);
    return oneOeps_sq*pow(vf::lbm::constant::cPi,-1.5f)*exp(-dist*dist*oneOeps_sq);
}


__host__ __device__ uint find_nearest_cellBSW(uint index, 
                                              real* coordsX, real* coordsY, real* coordsZ, 
                                              real posX, real posY, real posZ, 
                                              uint* neighborsX, uint*neighborsY, uint* neighborsZ, uint* neighborsWSB);

__global__ void interpolateVelocities(real* gridCoordsX, real* gridCoordsY, real* gridCoordsZ, 
                                      uint* neighborsX, uint* neighborsY, uint* neighborsZ, 
                                      uint* neighborsWSB, 
                                      real* vx, real* vy, real* vz, 
                                      int numberOfIndices, 
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
    // if(node==0 or node==90)
    // {
    //     printf("before: blade (%f, %f, %f), node BSW (%f, %f, %f), nodeTNE (%f, %f, %f)\n", bladePosX, bladePosY, bladePosZ, gridCoordsX[old_index], gridCoordsY[old_index], gridCoordsZ[old_index], gridCoordsX[neighborsX[old_index]], gridCoordsY[neighborsY[old_index]], gridCoordsZ[neighborsZ[old_index]]);
    // }
    uint kBSW = find_nearest_cellBSW(old_index, gridCoordsX, gridCoordsY, gridCoordsZ, bladePosX, bladePosY, bladePosZ, neighborsX, neighborsY, neighborsZ, neighborsWSB);
    
    bladeIndices[node] = kBSW;
    
    real u_interpX = 0.0;
    real u_interpY = 0.0;
    real u_interpZ = 0.0;

    bladeVelocitiesX[node] = u_interpX;
    bladeVelocitiesY[node] = u_interpY;
    bladeVelocitiesZ[node] = u_interpZ;

    // if(node==numberOfNodes-1)
    // {
    //     printf("after: blade (%f, %f, %f), node BSW (%f, %f, %f), nodeTNE (%f, %f, %f)\n", bladePosX, bladePosY, bladePosZ, gridCoordsX[kBSW], gridCoordsY[kBSW], gridCoordsZ[kBSW], gridCoordsX[neighborsX[kBSW]], gridCoordsY[neighborsY[kBSW]], gridCoordsZ[neighborsZ[kBSW]]);
    // }

}

__global__ void applyBodyForces(real* gridCoordsX, real* gridCoordsY, real* gridCoordsZ,
                           real* gridForcesX, real* gridForcesY, real* gridForcesZ, 
                           int* gridIndices, int numberOfIndices, 
                           real* bladeCoordsX, real* bladeCoordsY, real* bladeCoordsZ, 
                           real* bladeForcesX, real* bladeForcesY,real* bladeForcesZ,
                           real* bladeRadii,
                           real radius,
                           int numberOfNodes,
                           real epsilon)
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint index = nx*(ny*z + y) + x;

    if(index>=numberOfIndices) return;

    int gridIndex = gridIndices[index];

    real posX = gridCoordsX[gridIndex];
    real posY = gridCoordsY[gridIndex];
    real posZ = gridCoordsZ[gridIndex];

    real fXYZ_X = 0.0f;
    real fXYZ_Y = 0.0f;
    real fXYZ_Z = 0.0f;

    real last_r = 0.0f;
    real eta = 0.0f;
    real r = 0.0f;

    for( uint node=0; node<numberOfNodes; node++)
    {
        eta = calc_gaussian3D(posX, posY, posZ, bladeCoordsX[node], bladeCoordsY[node], bladeCoordsZ[node], epsilon);
        r = bladeRadii[node];
        fXYZ_X += bladeForcesX[node]*(r-last_r)*eta;
        fXYZ_Y += bladeForcesY[node]*(r-last_r)*eta;
        fXYZ_Z += bladeForcesZ[node]*(r-last_r)*eta;

        last_r = r;
    }

    fXYZ_X += bladeForcesX[numberOfNodes-1]*(radius-last_r)*eta;
    fXYZ_Y += bladeForcesY[numberOfNodes-1]*(radius-last_r)*eta;
    fXYZ_Z += bladeForcesZ[numberOfNodes-1]*(radius-last_r)*eta;

    gridForcesX[gridIndex] = fXYZ_X;
    gridForcesY[gridIndex] = fXYZ_Y;
    gridForcesZ[gridIndex] = fXYZ_Z;

}

void ActuatorLine::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{
    this->allocBladeRadii(cudaManager);
    this->allocBladeCoords(cudaManager);
    this->allocBladeVelocities(cudaManager);
    this->allocBladeForces(cudaManager);
    this->allocBladeIndices(cudaManager);

    this->initBladeRadii();
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
    vf::gpu::CudaGrid bladeGrid = vf::gpu::CudaGrid(numberOfThreads, this->numberOfNodes);

    interpolateVelocities<<< bladeGrid.grid, bladeGrid.threads >>>(
        para->getParD(this->level)->coordX_SP, para->getParD(this->level)->coordY_SP, para->getParD(this->level)->coordZ_SP,        
        para->getParD(this->level)->neighborX_SP, para->getParD(this->level)->neighborY_SP, para->getParD(this->level)->neighborZ_SP, para->getParD(this->level)->neighborWSB_SP,
        para->getParD(this->level)->vx_SP, para->getParD(this->level)->vy_SP, para->getParD(this->level)->vz_SP,
        this->numberOfIndices,
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeVelocitiesXD, this->bladeVelocitiesYD, this->bladeVelocitiesZD,  
        this->bladeIndicesD, this->numberOfNodes);

    this->copyBladeVelocitiesDtoH();

    if(true)
    {
        real forceRatio = para->getVelocityRatio()*para->getDensityRatio();
        for( int blade=0; blade<this->nBlades; blade++)
        {
            real localAzimuth = this->azimuth+2*blade*vf::lbm::constant::cPi/this->nBlades;
            for( uint bladeNode=0; bladeNode<this->nBladeNodes; bladeNode++)
            {
                uint node = bladeNode*blade+this->nBladeNodes;
                real uXYZ_X = this->bladeVelocitiesXH[node]/para->getVelocityRatio();
                real uXYZ_Y = this->bladeVelocitiesYH[node]/para->getVelocityRatio();
                real uXYZ_Z = this->bladeVelocitiesZH[node]/para->getVelocityRatio();

                real uRTZ_X, uRTZ_Y, uRTZ_Z;
                invRotateAboutX3D(uXYZ_X, uXYZ_Y, uXYZ_Z, localAzimuth, uRTZ_X, uRTZ_Y, uRTZ_Z);
                
                real u_rel = uRTZ_X;
                real v_rel = uRTZ_Y+this->omega*this->bladeRadiiH[node];

                real u_rel_sq = u_rel*u_rel+v_rel*v_rel;
                real phi = atan2(u_rel, v_rel);
                real Cl = 1.f;
                real Cd = 0.f;
                real c0 = 0.1f;
                real c = c0 * sqrt( 1- pow(2.f*(2*this->bladeRadiiH[node]-0.5*this->diameter)/(this->diameter), 2.f) );
                real Cn =   Cl*cos(phi)+Cd*sin(phi);
                real Ct =  -Cl*sin(phi)+Cd*cos(phi);

                real fRTZ_X = 0.5f*u_rel_sq*c*this->density*Cn;
                real fRTZ_Y = 0.5f*u_rel_sq*c*this->density*Ct;
                real fRTZ_Z = 0.0;

                real fXYZ_X, fXYZ_Y, fXYZ_Z;

                rotateAboutX3D(fRTZ_X, fRTZ_Y, fRTZ_Z, localAzimuth, fXYZ_X, fXYZ_Y, fXYZ_Z);
            
                this->bladeForcesXH[node] = fXYZ_X*forceRatio;
                this->bladeForcesYH[node] = fXYZ_Y*forceRatio;
                this->bladeForcesZH[node] = fXYZ_Z*forceRatio;
            }
        }
    }

    this->copyBladeForcesHtoD();

    vf::gpu::CudaGrid sphereGrid = vf::gpu::CudaGrid(numberOfThreads, this->numberOfIndices);

    applyBodyForces<<<sphereGrid.grid, sphereGrid.threads>>>(
        para->getParD(this->level)->coordX_SP, para->getParD(this->level)->coordY_SP, para->getParD(this->level)->coordZ_SP,        
        para->getParD(this->level)->forceX_SP, para->getParD(this->level)->forceY_SP, para->getParD(this->level)->forceZ_SP,        
        this->boundingSphereIndicesD, this->numberOfIndices,
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeForcesXD, this->bladeForcesYD, this->bladeForcesZD,
        this->bladeRadiiD,
        this->diameter*0.5f,  
        this->numberOfNodes,
        this->epsilon);

    real dazimuth = this->omega*this->delta_t;

    this->azimuth += dazimuth;
    this->rotateBlades(dazimuth);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Blade 

void ActuatorLine::allocBladeRadii(CudaMemoryManager* cudaManager)
{
    checkCudaErrors( cudaMallocHost((void**) &this->bladeRadiiH, sizeof(real)*this->nBladeNodes) );

    checkCudaErrors( cudaMalloc((void**) &this->bladeRadiiD, sizeof(real)*this->nBladeNodes) );

    cudaManager->setMemsizeGPU(sizeof(real)*this->nBladeNodes, false);
}

void ActuatorLine::initBladeRadii()
{   
    real dx = 0.5f*this->diameter/this->nBladeNodes;        

    for(unsigned int node=0; node<this->nBladeNodes; node++)
    {
        this->bladeRadiiH[node] = dx*(node+1);
    }

}

void ActuatorLine::copyBladeRadiiHtoD()
{
    checkCudaErrors( cudaMemcpy(this->bladeRadiiD, this->bladeRadiiH, sizeof(real)*this->nBladeNodes, cudaMemcpyHostToDevice) );
    }

void ActuatorLine::copyBladeRadiiDtoH()
{
    checkCudaErrors( cudaMemcpy(this->bladeRadiiH, this->bladeRadiiD, sizeof(real)*this->nBladeNodes, cudaMemcpyDeviceToHost) );
    }

void ActuatorLine::freeBladeRadii()
{
    checkCudaErrors( cudaFree(this->bladeRadiiD) );
    
    checkCudaErrors( cudaFreeHost(this->bladeRadiiH) );
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
    for( unsigned int blade=0; blade<this->nBlades; blade++)
    {
        real localAzimuth = this->azimuth*(2*vf::lbm::constant::cPi/this->nBlades*blade);
        for(unsigned int node=0; node<this->nBladeNodes; node++)
        {
            real coordY, coordZ, r, y;
            r = this->bladeRadiiH[node];
            y = 0.f;

            rotate2D(localAzimuth, y, r, coordY, coordZ);
            this->bladeCoordsXH[node+this->nBladeNodes*blade] = this->turbinePosX;
            this->bladeCoordsYH[node+this->nBladeNodes*blade] = this->turbinePosY+coordY;
            this->bladeCoordsZH[node+this->nBladeNodes*blade] = this->turbinePosZ+coordZ;
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
        real oldCoordY = this->bladeCoordsYH[node];
        real oldCoordZ = this->bladeCoordsZH[node];
        real newCoordY, newCoordZ;
        rotate2D(angle, oldCoordY, oldCoordZ, newCoordY, newCoordZ, this->turbinePosY, this->turbinePosZ);
        this->bladeCoordsYH[node] = newCoordY;
        this->bladeCoordsZH[node] = newCoordZ;
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

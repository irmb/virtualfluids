#include "ActuatorLine.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cuda/CudaGrid.h>
#include "lbm/constants/NumericConstants.h"
#include "VirtualFluids_GPU/GPU/GeometryUtils.h"

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"

__host__ __device__ __inline__ uint calcNode(uint bladeNode, uint nBladeNodes, uint blade, uint nBlades)
{
    return bladeNode+blade*nBladeNodes;
}

__host__ __device__ __inline__ void calcBladeAndBladeNode(uint node, uint& bladeNode, uint nBladeNodes, uint& blade, uint nBlades)
{
    blade = node/nBladeNodes;
    bladeNode = node % nBladeNodes;
}

__host__ __device__ __inline__ real calcGaussian3D(real posX, real posY, real posZ, real destX, real destY, real destZ, real epsilon)
{
    real distX = destX-posX;
    real distY = destY-posY;
    real distZ = destZ-posZ;
    real dist = sqrt(distX*distX+distY*distY+distZ*distZ);
    return pow(epsilon,-3)*pow(vf::lbm::constant::cPi,-1.5f)*exp(-pow(dist/epsilon,2));
}

__host__ __device__ __inline__ void rotateFromBladeToGlobal(
                            real& bladeCoordX_BF, real& bladeCoordY_BF, real& bladeCoordZ_BF, 
                            real& bladeCoordX_GF, real& bladeCoordY_GF, real& bladeCoordZ_GF,
                            real& azimuth, real& yaw)
{
    real tmpX, tmpY, tmpZ;

    rotateAboutX3D(azimuth, bladeCoordX_BF, bladeCoordY_BF, bladeCoordZ_BF, tmpX, tmpY, tmpZ);
    rotateAboutZ3D(yaw, tmpX, tmpY, tmpZ, bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF);

}

__host__ __device__ __inline__ void rotateFromGlobalToBlade(
                            real& bladeCoordX_BF, real& bladeCoordY_BF, real& bladeCoordZ_BF, 
                            real& bladeCoordX_GF, real& bladeCoordY_GF, real& bladeCoordZ_GF,
                            real& azimuth, real& yaw)
{
    real tmpX, tmpY, tmpZ;

    invRotateAboutZ3D(yaw, bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF, tmpX, tmpY, tmpZ);
    invRotateAboutX3D(azimuth, tmpX, tmpY, tmpZ, bladeCoordX_BF, bladeCoordY_BF, bladeCoordZ_BF);
}

__global__ void interpolateVelocities(real* gridCoordsX, real* gridCoordsY, real* gridCoordsZ, 
                                      uint* neighborsX, uint* neighborsY, uint* neighborsZ, uint* neighborsWSB, 
                                      real* vx, real* vy, real* vz, 
                                      real* bladeCoordsX, real* bladeCoordsY, real* bladeCoordsZ,
                                      real* bladeVelocitiesX, real* bladeVelocitiesY, real* bladeVelocitiesZ, 
                                      uint nBlades, uint nBladeNodes, 
                                      real azimuth, real yaw, real omega, 
                                      real turbPosX, real turbPosZ, real turbPosY,
                                      uint* bladeIndices, real velocityRatio)
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint node = nx*(ny*z + y) + x;

    uint bladeNode, blade;

    calcBladeAndBladeNode(node, bladeNode, nBladeNodes, blade, nBlades);

    if(node>=nBladeNodes*nBlades) return;

    real bladePosX_BF = bladeCoordsX[node];
    real bladePosY_BF = bladeCoordsY[node];
    real bladePosZ_BF = bladeCoordsZ[node];

    real bladePosX_GF, bladePosY_GF, bladePosZ_GF;

    real localAzimuth = azimuth+blade*2*vf::lbm::constant::cPi/nBlades;

    rotateFromBladeToGlobal(bladePosX_BF, bladePosY_BF, bladePosZ_BF, 
                            bladePosX_GF, bladePosY_GF, bladePosZ_GF,
                            localAzimuth, yaw);

    bladePosX_GF += turbPosX;
    bladePosY_GF += turbPosY;
    bladePosZ_GF += turbPosZ;

    uint old_index = bladeIndices[node];

    uint k, ke, kn, kt;
    uint kne, kte, ktn, ktne;

    k = findNearestCellBSW(old_index, 
                           gridCoordsX, gridCoordsY, gridCoordsZ, 
                           bladePosX_GF, bladePosY_GF, bladePosZ_GF, 
                           neighborsX, neighborsY, neighborsZ, neighborsWSB);
        
    bladeIndices[node] = k;

    getNeighborIndicesOfBSW(k, ke, kn, kt, kne, kte, ktn, ktne, neighborsX, neighborsY, neighborsZ);

    real dW, dE, dN, dS, dT, dB;

    real invDeltaX = 1.f/(gridCoordsX[ktne]-gridCoordsX[k]);
    real distX = invDeltaX*(bladePosX_GF-gridCoordsX[k]);
    real distY = invDeltaX*(bladePosY_GF-gridCoordsY[k]);
    real distZ = invDeltaX*(bladePosZ_GF-gridCoordsZ[k]);

    getInterpolationWeights(dW, dE, dN, dS, dT, dB, distX, distY, distZ);

    real bladeVelX_GF = trilinearInterpolation(dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vx)*velocityRatio;
    real bladeVelY_GF = trilinearInterpolation(dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vy)*velocityRatio;
    real bladeVelZ_GF = trilinearInterpolation(dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vz)*velocityRatio;

    real bladeVelX_BF, bladeVelY_BF, bladeVelZ_BF;

    rotateFromGlobalToBlade(bladeVelX_BF, bladeVelY_BF, bladeVelZ_BF, bladeVelX_GF, bladeVelY_GF, bladeVelZ_GF, localAzimuth, yaw);

    bladeVelocitiesX[node] = bladeVelX_BF;
    bladeVelocitiesY[node] = bladeVelY_BF+omega*bladePosZ_BF;
    bladeVelocitiesZ[node] = bladeVelZ_BF;
}


__global__ void applyBodyForces(real* gridCoordsX, real* gridCoordsY, real* gridCoordsZ,
                                real* gridForcesX, real* gridForcesY, real* gridForcesZ, 
                                real* bladeCoordsX, real* bladeCoordsY, real* bladeCoordsZ, 
                                real* bladeForcesX, real* bladeForcesY,real* bladeForcesZ,
                                uint nBlades, uint nBladeNodes,
                                real azimuth, real yaw, real omega, 
                                real turbPosX, real turbPosZ, real turbPosY,
                                real* bladeRadii, real forceRatio,
                                uint* gridIndices, uint numberOfIndices, 
                                real radius, real epsilon, real delta_x)
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

    real gridForceX = 0.0f;
    real gridForceY = 0.0f;
    real gridForceZ = 0.0f;

    real deltaXCubed = pow(delta_x,3);

    real dAzimuth = 2*vf::lbm::constant::cPi/nBlades;

    for( uint blade=0; blade<nBlades; blade++)
    { 
        real last_r = 0.0f;
        real r = 0.0f;
        real localAzimuth = azimuth+blade*dAzimuth;
        
        for( uint bladeNode=0; bladeNode<nBladeNodes; bladeNode++)
        {
            r = bladeRadii[bladeNode];

            uint node = calcNode(bladeNode, nBladeNodes, blade, nBlades);

            real bladePosX_GF, bladePosY_GF, bladePosZ_GF;

            rotateFromBladeToGlobal(bladeCoordsX[node], bladeCoordsY[node], bladeCoordsZ[node], 
                                    bladePosX_GF, bladePosY_GF, bladePosZ_GF,
                                    localAzimuth, yaw);

            bladePosX_GF += turbPosX;
            bladePosY_GF += turbPosY;
            bladePosZ_GF += turbPosZ;

            real eta = (r-last_r)*calcGaussian3D(posX, posY, posZ, bladePosX_GF, bladePosY_GF, bladePosZ_GF, epsilon)*deltaXCubed;

            real forceX_GF, forceY_GF, forceZ_GF;

            rotateFromBladeToGlobal(bladeForcesX[node], bladeForcesY[node], bladeForcesZ[node], forceX_GF, forceY_GF, forceZ_GF, localAzimuth, yaw);
            
            gridForceX += forceX_GF*eta;
            gridForceY += forceY_GF*eta;
            gridForceZ += forceZ_GF*eta;

            last_r = r;
        }         
    }

    real invForceRatio = 1.f/forceRatio;

    gridForcesX[gridIndex] = gridForceX*invForceRatio;
    gridForcesY[gridIndex] = gridForceY*invForceRatio;
    gridForcesZ[gridIndex] = gridForceZ*invForceRatio;
}


void ActuatorLine::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{
    if(!para->getIsBodyForce()) throw std::runtime_error("try to allocate ActuatorLine but BodyForce is not set in Parameter.");
    this->initBladeRadii(cudaManager);
    this->initBladeCoords(cudaManager);    
    this->initBladeIndices(para, cudaManager);
    this->initBladeVelocities(cudaManager);
    this->initBladeForces(cudaManager);    
    this->initBoundingSphere(para, cudaManager);
}


void ActuatorLine::interact(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t)
{
    if (level != this->level) return;

    cudaManager->cudaCopyBladeCoordsHtoD(this);

    uint numberOfThreads = para->getParH(level)->numberofthreads;
    vf::cuda::CudaGrid bladeGrid = vf::cuda::CudaGrid(numberOfThreads, this->numberOfNodes);

    interpolateVelocities<<< bladeGrid.grid, bladeGrid.threads >>>(
        para->getParD(this->level)->coordX_SP, para->getParD(this->level)->coordY_SP, para->getParD(this->level)->coordZ_SP,        
        para->getParD(this->level)->neighborX_SP, para->getParD(this->level)->neighborY_SP, para->getParD(this->level)->neighborZ_SP, para->getParD(this->level)->neighborWSB_SP,
        para->getParD(this->level)->vx_SP, para->getParD(this->level)->vy_SP, para->getParD(this->level)->vz_SP,
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeVelocitiesXD, this->bladeVelocitiesYD, this->bladeVelocitiesZD,  
        this->nBlades, this->nBladeNodes,
        this->azimuth, this->yaw, this->omega, 
        this->turbinePosX, this->turbinePosZ, this->turbinePosY,
        this->bladeIndicesD, para->getVelocityRatio());

    cudaManager->cudaCopyBladeVelocitiesDtoH(this);

    this->calcBladeForces();

    cudaManager->cudaCopyBladeForcesHtoD(this);

    vf::cuda::CudaGrid sphereGrid = vf::cuda::CudaGrid(numberOfThreads, this->numberOfIndices);

    applyBodyForces<<<sphereGrid.grid, sphereGrid.threads>>>(
        para->getParD(this->level)->coordX_SP, para->getParD(this->level)->coordY_SP, para->getParD(this->level)->coordZ_SP,        
        para->getParD(this->level)->forceX_SP, para->getParD(this->level)->forceY_SP, para->getParD(this->level)->forceZ_SP,        
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeForcesXD, this->bladeForcesYD, this->bladeForcesZD,
        this->nBlades, this->nBladeNodes,
        this->azimuth, this->yaw, this->omega, 
        this->turbinePosX, this->turbinePosZ, this->turbinePosY,
        this->bladeRadiiD, this->density*pow(para->getViscosityRatio(),2),
        this->boundingSphereIndicesD, this->numberOfIndices,
        this->diameter*0.5f, this->epsilon, this->delta_x);

    this->azimuth = fmod(this->azimuth+this->omega*this->delta_t,2*vf::lbm::constant::cPi);
}


void ActuatorLine::free(Parameter* para, CudaMemoryManager* cudaManager)
{
    cudaManager->cudaFreeBladeRadii(this);
    cudaManager->cudaFreeBladeCoords(this);
    cudaManager->cudaFreeBladeVelocities(this);
    cudaManager->cudaFreeBladeForces(this);
    cudaManager->cudaFreeBladeIndices(this);
    cudaManager->cudaFreeSphereIndices(this);
}


void ActuatorLine::calcForcesEllipticWing()
{
    uint node;
    real u_rel, v_rel, u_rel_sq;
    real phi;
    real Cl = 1.f;
    real Cd = 0.f;
    real c0 = 1.f;

    real c, Cn, Ct;

    for( uint blade=0; blade<this->nBlades; blade++)
    { 
        for( uint bladeNode=0; bladeNode<this->nBladeNodes; bladeNode++)
        {        
            node = calcNode(bladeNode, this->nBladeNodes, blade, this->nBlades);

            u_rel = this->bladeVelocitiesXH[node];
            v_rel = this->bladeVelocitiesYH[node];
            u_rel_sq = u_rel*u_rel+v_rel*v_rel;
            phi = atan2(u_rel, v_rel);

            c = c0 * sqrt( 1.f- pow(4.f*this->bladeRadiiH[bladeNode]/this->diameter-1.f, 2.f) );
            Cn =   Cl*cos(phi)+Cd*sin(phi);
            Ct =  -Cl*sin(phi)+Cd*cos(phi);
        
            this->bladeForcesXH[node] = -0.5f*u_rel_sq*c*this->density*Cn;
            this->bladeForcesYH[node] = -0.5f*u_rel_sq*c*this->density*Ct;
            this->bladeForcesZH[node] = 0.0f;
        }
    }
}

void ActuatorLine::calcBladeForces()
{
    this->calcForcesEllipticWing();
}

void ActuatorLine::initBladeRadii(CudaMemoryManager* cudaManager)
{   
    cudaManager->cudaAllocBladeRadii(this);

    real dx = 0.5f*this->diameter/this->nBladeNodes;  

    for(uint node=0; node<this->nBladeNodes; node++)
    {
        this->bladeRadiiH[node] = dx*(node+1);
    }
    cudaManager->cudaCopyBladeRadiiHtoD(this);
}

void ActuatorLine::initBladeCoords(CudaMemoryManager* cudaManager)
{   
    cudaManager->cudaAllocBladeCoords(this);

    for(uint blade=0; blade<this->nBlades; blade++)
    {
        for(uint bladeNode=0; bladeNode<this->nBladeNodes; bladeNode++)
        {
            uint node = calcNode(bladeNode, this->nBladeNodes, blade, this->nBlades);

            this->bladeCoordsXH[node] = 0.f;
            this->bladeCoordsYH[node] = 0.f;
            this->bladeCoordsZH[node] = this->bladeRadiiH[bladeNode];
        }
    }
    cudaManager->cudaCopyBladeCoordsHtoD(this);
}

void ActuatorLine::initBladeVelocities(CudaMemoryManager* cudaManager)
{   
    cudaManager->cudaAllocBladeVelocities(this);

    for(uint node=0; node<this->numberOfNodes; node++)
    {
        this->bladeVelocitiesXH[node] = 0.f;
        this->bladeVelocitiesYH[node] = 0.f;
        this->bladeVelocitiesZH[node] = 0.f;
    }
    cudaManager->cudaCopyBladeVelocitiesHtoD(this);
}

void ActuatorLine::initBladeForces(CudaMemoryManager* cudaManager)
{   
    cudaManager->cudaAllocBladeForces(this);

    for(uint node=0; node<this->numberOfNodes; node++)
    {
        this->bladeForcesXH[node] = 0.f;
        this->bladeForcesYH[node] = 0.f;
        this->bladeForcesZH[node] = 0.f;
    }
    cudaManager->cudaCopyBladeForcesHtoD(this);
}

void ActuatorLine::initBladeIndices(Parameter* para, CudaMemoryManager* cudaManager)
{   
    cudaManager->cudaAllocBladeIndices(this);

    for(uint node=0; node<this->numberOfNodes; node++)
    {
        this->bladeIndicesH[node] = 1;
    }
    cudaManager->cudaCopyBladeIndicesHtoD(this);
}

void ActuatorLine::initBoundingSphere(Parameter* para, CudaMemoryManager* cudaManager)
{
    // Actuator line exists only on 1 level
    std::vector<int> nodesInSphere;
    real sphereRadius = 0.5*this->diameter+4.f*this->epsilon;

    for (uint j = 1; j <= para->getParH(this->level)->size_Mat_SP; j++)
    {
        const real distX = para->getParH(this->level)->coordX_SP[j]-this->turbinePosX;
        const real distY = para->getParH(this->level)->coordY_SP[j]-this->turbinePosY;
        const real distZ = para->getParH(this->level)->coordZ_SP[j]-this->turbinePosZ;
        const real dist = sqrt(pow(distX,2)+pow(distY,2)+pow(distZ,2));
        if(dist < sphereRadius) nodesInSphere.push_back(j);
    }

    this->numberOfIndices = uint(nodesInSphere.size());
    cudaManager->cudaAllocSphereIndices(this);
    std::copy(nodesInSphere.begin(), nodesInSphere.end(), this->boundingSphereIndicesH);
    cudaManager->cudaCopySphereIndicesHtoD(this);
}

void ActuatorLine::setBladeCoords(real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ)
{ 

    for(uint node=0; node<this->numberOfNodes; node++)
    {
        this->bladeCoordsXH[node] = _bladeCoordsX[node];
        this->bladeCoordsYH[node] = _bladeCoordsY[node];
        this->bladeCoordsZH[node] = _bladeCoordsZ[node];
    }
}

void ActuatorLine::setBladeVelocities(real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ)
{ 
    for(uint node=0; node<this->numberOfNodes; node++)
    {
        this->bladeVelocitiesXH[node] = _bladeVelocitiesX[node];
        this->bladeVelocitiesYH[node] = _bladeVelocitiesY[node];
        this->bladeVelocitiesZH[node] = _bladeVelocitiesZ[node];
    }
}

void ActuatorLine::setBladeForces(real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ)
{ 
    for(uint node=0; node<this->numberOfNodes; node++)
    {
        this->bladeForcesXH[node] = _bladeForcesX[node];
        this->bladeForcesYH[node] = _bladeForcesY[node];
        this->bladeForcesZH[node] = _bladeForcesZ[node];
    }
}
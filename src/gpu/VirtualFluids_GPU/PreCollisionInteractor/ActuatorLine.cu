#include "ActuatorLine.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cuda/CudaGrid.h>
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
    bladeNode = node - blade*nBladeNodes;
}

__host__ __device__ __forceinline__ real distSqrd(real distX, real distY, real distZ)
{
    return distX*distX+distY*distY+distZ*distZ;
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
                                      real turbPosX, real turbPosY, real turbPosZ,
                                      uint* bladeIndices, real velocityRatio, real invDeltaX)
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

    real bladeCoordX_BF = bladeCoordsX[node];
    real bladeCoordY_BF = bladeCoordsY[node];
    real bladeCoordZ_BF = bladeCoordsZ[node];

    real bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF;

    real localAzimuth = azimuth+blade*c2Pi/nBlades;

    rotateFromBladeToGlobal(bladeCoordX_BF, bladeCoordY_BF, bladeCoordZ_BF, 
                            bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF,
                            localAzimuth, yaw);

    bladeCoordX_GF += turbPosX;
    bladeCoordY_GF += turbPosY;
    bladeCoordZ_GF += turbPosZ;

    uint old_index = bladeIndices[node];

    uint k, ke, kn, kt;
    uint kne, kte, ktn, ktne;

    k = findNearestCellBSW(old_index, 
                           gridCoordsX, gridCoordsY, gridCoordsZ, 
                           bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF, 
                           neighborsX, neighborsY, neighborsZ, neighborsWSB);
        
    bladeIndices[node] = k;

    getNeighborIndicesOfBSW(k, ke, kn, kt, kne, kte, ktn, ktne, neighborsX, neighborsY, neighborsZ);

    real dW, dE, dN, dS, dT, dB;

    real distX = invDeltaX*(bladeCoordX_GF-gridCoordsX[k]);
    real distY = invDeltaX*(bladeCoordY_GF-gridCoordsY[k]);
    real distZ = invDeltaX*(bladeCoordZ_GF-gridCoordsZ[k]);

    getInterpolationWeights(dW, dE, dN, dS, dT, dB, distX, distY, distZ);

    real bladeVelX_GF = trilinearInterpolation(dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vx)*velocityRatio;
    real bladeVelY_GF = trilinearInterpolation(dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vy)*velocityRatio;
    real bladeVelZ_GF = trilinearInterpolation(dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vz)*velocityRatio;

    real bladeVelX_BF, bladeVelY_BF, bladeVelZ_BF;

    rotateFromGlobalToBlade(bladeVelX_BF, bladeVelY_BF, bladeVelZ_BF, bladeVelX_GF, bladeVelY_GF, bladeVelZ_GF, localAzimuth, yaw);

    bladeVelocitiesX[node] = bladeVelX_BF;
    bladeVelocitiesY[node] = bladeVelY_BF+omega*bladeCoordZ_BF;
    bladeVelocitiesZ[node] = bladeVelZ_BF;
}


__global__ void applyBodyForces(real* gridCoordsX, real* gridCoordsY, real* gridCoordsZ,
                                real* gridForcesX, real* gridForcesY, real* gridForcesZ, 
                                real* bladeCoordsX, real* bladeCoordsY, real* bladeCoordsZ, 
                                real* bladeForcesX, real* bladeForcesY,real* bladeForcesZ,
                                uint nBlades, uint nBladeNodes,
                                real azimuth, real yaw, real omega, 
                                real turbPosX, real turbPosY, real turbPosZ,
                                uint* gridIndices, uint nIndices, 
                                real invEpsilonSqrd, real factorGaussian)
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint index = nx*(ny*z + y) + x;

    if(index>=nIndices) return;

    int gridIndex = gridIndices[index];

    real gridCoordX_RF = gridCoordsX[gridIndex] - turbPosX;
    real gridCoordY_RF = gridCoordsY[gridIndex] - turbPosY;
    real gridCoordZ_RF = gridCoordsZ[gridIndex] - turbPosZ;

    real gridForceX_RF = c0o1;
    real gridForceY_RF = c0o1;
    real gridForceZ_RF = c0o1;

    real dAzimuth = c2Pi/nBlades;

    for( uint blade=0; blade<nBlades; blade++)
    { 
        real localAzimuth = azimuth+blade*dAzimuth;
        
        for( uint bladeNode=0; bladeNode<nBladeNodes; bladeNode++)
        {
            uint node = calcNode(bladeNode, nBladeNodes, blade, nBlades);

            real bladeCoordX_RF, bladeCoordY_RF, bladeCoordZ_RF;

            rotateFromBladeToGlobal(bladeCoordsX[node], bladeCoordsY[node], bladeCoordsZ[node], 
                                    bladeCoordX_RF, bladeCoordY_RF, bladeCoordZ_RF,
                                    localAzimuth, yaw);

            real eta = factorGaussian*exp(-distSqrd(bladeCoordX_RF-gridCoordX_RF, bladeCoordY_RF-gridCoordY_RF, bladeCoordZ_RF-gridCoordZ_RF)*invEpsilonSqrd);

            real forceX_RF, forceY_RF, forceZ_RF;

            rotateFromBladeToGlobal(bladeForcesX[node], bladeForcesY[node], bladeForcesZ[node], forceX_RF, forceY_RF, forceZ_RF, localAzimuth, yaw);
            
            gridForceX_RF += forceX_RF*eta;
            gridForceY_RF += forceY_RF*eta;
            gridForceZ_RF += forceZ_RF*eta;
        }
    }

    atomicAdd(&gridForcesX[gridIndex],gridForceX_RF);
    atomicAdd(&gridForcesY[gridIndex],gridForceY_RF);
    atomicAdd(&gridForcesZ[gridIndex],gridForceZ_RF);
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

    vf::cuda::CudaGrid bladeGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->nNodes);

    interpolateVelocities<<< bladeGrid.grid, bladeGrid.threads >>>(
        para->getParD(this->level)->coordX_SP, para->getParD(this->level)->coordY_SP, para->getParD(this->level)->coordZ_SP,        
        para->getParD(this->level)->neighborX_SP, para->getParD(this->level)->neighborY_SP, para->getParD(this->level)->neighborZ_SP, para->getParD(this->level)->neighborWSB_SP,
        para->getParD(this->level)->vx_SP, para->getParD(this->level)->vy_SP, para->getParD(this->level)->vz_SP,
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeVelocitiesXD, this->bladeVelocitiesYD, this->bladeVelocitiesZD,  
        this->nBlades, this->nBladeNodes,
        this->azimuth, this->yaw, this->omega, 
        this->turbinePosX, this->turbinePosY, this->turbinePosZ,
        this->bladeIndicesD, para->getVelocityRatio(), this->invDeltaX);

    cudaManager->cudaCopyBladeVelocitiesDtoH(this);

    this->calcBladeForces();

    cudaManager->cudaCopyBladeForcesHtoD(this);

    vf::cuda::CudaGrid sphereGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->nIndices);

    applyBodyForces<<<sphereGrid.grid, sphereGrid.threads>>>(
        para->getParD(this->level)->coordX_SP, para->getParD(this->level)->coordY_SP, para->getParD(this->level)->coordZ_SP,        
        para->getParD(this->level)->forceX_SP, para->getParD(this->level)->forceY_SP, para->getParD(this->level)->forceZ_SP,        
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeForcesXD, this->bladeForcesYD, this->bladeForcesZD,
        this->nBlades, this->nBladeNodes,
        this->azimuth, this->yaw, this->omega, 
        this->turbinePosX, this->turbinePosY, this->turbinePosZ,
        this->boundingSphereIndicesD, this->nIndices,
        this->invEpsilonSqrd, this->factorGaussian);

    this->azimuth = fmod(this->azimuth+this->omega*this->deltaT,c2Pi);
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
    real Cl = c1o1;
    real Cd = c0o1;
    real c0 = c1o1;

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
            
            real tmp = c4o1*this->bladeRadiiH[bladeNode]/this->diameter-c1o1;
            c = c0 * sqrt( c1o1- tmp*tmp );
            Cn = Cl*cos(phi)+Cd*sin(phi);
            Ct = Cl*sin(phi)-Cd*cos(phi);
        
            this->bladeForcesXH[node] = -c1o2*u_rel_sq*c*this->density*Cn;
            this->bladeForcesYH[node] = -c1o2*u_rel_sq*c*this->density*Ct;
            this->bladeForcesZH[node] = c0o1;
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

    real dr = c1o2*this->diameter/this->nBladeNodes;  

    for(uint node=0; node<this->nBladeNodes; node++)
    {
        this->bladeRadiiH[node] = dr*(node+1);
    }
    cudaManager->cudaCopyBladeRadiiHtoD(this);

    real dxOPiSqrtEps = this->deltaX/(this->epsilon*sqrt(cPi));
    this->factorGaussian = dr*dxOPiSqrtEps*dxOPiSqrtEps*dxOPiSqrtEps/this->forceRatio;
}

void ActuatorLine::initBladeCoords(CudaMemoryManager* cudaManager)
{   
    cudaManager->cudaAllocBladeCoords(this);

    for(uint blade=0; blade<this->nBlades; blade++)
    {
        for(uint bladeNode=0; bladeNode<this->nBladeNodes; bladeNode++)
        {
            uint node = calcNode(bladeNode, this->nBladeNodes, blade, this->nBlades);

            this->bladeCoordsXH[node] = c0o1;
            this->bladeCoordsYH[node] = c0o1;
            this->bladeCoordsZH[node] = this->bladeRadiiH[bladeNode];
        }
    }
    cudaManager->cudaCopyBladeCoordsHtoD(this);
}

void ActuatorLine::initBladeVelocities(CudaMemoryManager* cudaManager)
{   
    cudaManager->cudaAllocBladeVelocities(this);

    for(uint node=0; node<this->nNodes; node++)
    {
        this->bladeVelocitiesXH[node] = c0o1;
        this->bladeVelocitiesYH[node] = c0o1;
        this->bladeVelocitiesZH[node] = c0o1;
    }
    cudaManager->cudaCopyBladeVelocitiesHtoD(this);
}

void ActuatorLine::initBladeForces(CudaMemoryManager* cudaManager)
{   
    cudaManager->cudaAllocBladeForces(this);

    for(uint node=0; node<this->nNodes; node++)
    {
        this->bladeForcesXH[node] = c0o1;
        this->bladeForcesYH[node] = c0o1;
        this->bladeForcesZH[node] = c0o1;
    }
    cudaManager->cudaCopyBladeForcesHtoD(this);
}

void ActuatorLine::initBladeIndices(Parameter* para, CudaMemoryManager* cudaManager)
{   
    cudaManager->cudaAllocBladeIndices(this);

    for(uint node=0; node<this->nNodes; node++)
    {
        this->bladeIndicesH[node] = 1;
    }
    cudaManager->cudaCopyBladeIndicesHtoD(this);
}

void ActuatorLine::initBoundingSphere(Parameter* para, CudaMemoryManager* cudaManager)
{
    // Actuator line exists only on 1 level
    std::vector<int> nodesInSphere;
    real sphereRadius = c1o2*this->diameter+c4o1*this->epsilon;
    real sphereRadiusSqrd = sphereRadius*sphereRadius;

    for (uint j = 1; j <= para->getParH(this->level)->size_Mat_SP; j++)
    {
        const real distX = para->getParH(this->level)->coordX_SP[j]-this->turbinePosX;
        const real distY = para->getParH(this->level)->coordY_SP[j]-this->turbinePosY;
        const real distZ = para->getParH(this->level)->coordZ_SP[j]-this->turbinePosZ;
        if(distSqrd(distX,distY,distZ) < sphereRadiusSqrd) nodesInSphere.push_back(j);
    }

    this->nIndices = uint(nodesInSphere.size());
    cudaManager->cudaAllocSphereIndices(this);
    std::copy(nodesInSphere.begin(), nodesInSphere.end(), this->boundingSphereIndicesH);
    cudaManager->cudaCopySphereIndicesHtoD(this);
}

void ActuatorLine::setBladeCoords(real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ)
{ 

    for(uint node=0; node<this->nNodes; node++)
    {
        this->bladeCoordsXH[node] = _bladeCoordsX[node];
        this->bladeCoordsYH[node] = _bladeCoordsY[node];
        this->bladeCoordsZH[node] = _bladeCoordsZ[node];
    }
}

void ActuatorLine::setBladeVelocities(real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ)
{ 
    for(uint node=0; node<this->nNodes; node++)
    {
        this->bladeVelocitiesXH[node] = _bladeVelocitiesX[node];
        this->bladeVelocitiesYH[node] = _bladeVelocitiesY[node];
        this->bladeVelocitiesZH[node] = _bladeVelocitiesZ[node];
    }
}

void ActuatorLine::setBladeForces(real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ)
{ 
    for(uint node=0; node<this->nNodes; node++)
    {
        this->bladeForcesXH[node] = _bladeForcesX[node];
        this->bladeForcesYH[node] = _bladeForcesY[node];
        this->bladeForcesZH[node] = _bladeForcesZ[node];
    }
}
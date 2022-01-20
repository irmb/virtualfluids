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
#include "Kernel/Utilities/DistributionHelper.cuh"
#include "lbm/MacroscopicQuantities.h"

__host__ __device__ __inline__ real calcGaussian3D(real posX, real posY, real posZ, real destX, real destY, real destZ, real epsilon)
{
    real distX = destX-posX;
    real distY = destY-posY;
    real distZ = destZ-posZ;
    real dist = sqrt(distX*distX+distY*distY+distZ*distZ);
    return pow(epsilon,-3)*pow(vf::lbm::constant::cPi,-1.5f)*exp(-pow(dist/epsilon,2));
}


__global__ void interpolateVelocities(real* gridCoordsX, real* gridCoordsY, real* gridCoordsZ, 
                                      uint* neighborX, uint* neighborY, uint* neighborZ, 
                                      uint* neighborsWSB, 
                                      real* distributions, 
                                      uint size_Mat, bool isEvenTimestep,
                                      uint numberOfIndices, 
                                      real* bladeCoordsX, real* bladeCoordsY, real* bladeCoordsZ, 
                                      real* bladeVelocitiesX, real* bladeVelocitiesY, real* bladeVelocitiesZ, 
                                      uint* bladeIndices, uint numberOfNodes)
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

    uint k, ke, kn, kt;
    uint kne, kte, ktn, ktne;

    k = findNearestCellBSW(old_index, 
                           gridCoordsX, gridCoordsY, gridCoordsZ, 
                           bladePosX, bladePosY, bladePosZ, 
                           neighborX, neighborY, neighborZ, neighborsWSB);
        
    bladeIndices[node] = k;

    getNeighborIndicesOfBSW(k, ke, kn, kt, kne, kte, ktn, ktne, neighborX, neighborY, neighborZ);

    real dW, dE, dN, dS, dT, dB;

    real invDeltaX = 1.f/(gridCoordsX[ktne]-gridCoordsX[k]);
    real distX = invDeltaX*(gridCoordsX[ktne]-bladePosX);
    real distY = invDeltaX*(gridCoordsY[ktne]-bladePosY);
    real distZ = invDeltaX*(gridCoordsZ[ktne]-bladePosZ);

    getInterpolationWeights(dW, dE, dN, dS, dT, dB, 
                            distX, distY, distZ);

    real vx[8], vy[8], vz[8];
    uint ks[8] = {k, ke, kn, kt, kne, kte, ktn, ktne};

    for(uint i=0; i<8; i++)
    {
        vf::gpu::DistributionWrapper distr_wrapper(distributions, size_Mat, isEvenTimestep, ks[i], neighborX, neighborY, neighborZ);
        const auto& distribution = distr_wrapper.distribution;

        real rho = vf::lbm::getDensity(distribution.f);
        vx[i] = vf::lbm::getCompressibleVelocityX1(distribution.f, rho);
        vy[i] = vf::lbm::getCompressibleVelocityX2(distribution.f, rho);
        vz[i] = vf::lbm::getCompressibleVelocityX3(distribution.f, rho);
    }
    
    bladeVelocitiesX[node] = trilinearInterpolation(dW, dE, dN, dS, dT, dB, 0, 1, 2, 3, 4, 5, 6, 7, vx);
    bladeVelocitiesY[node] = trilinearInterpolation(dW, dE, dN, dS, dT, dB, 0, 1, 2, 3, 4, 5, 6, 7, vy);
    bladeVelocitiesZ[node] = trilinearInterpolation(dW, dE, dN, dS, dT, dB, 0, 1, 2, 3, 4, 5, 6, 7, vz);

}


__global__ void applyBodyForces(real* gridCoordsX, real* gridCoordsY, real* gridCoordsZ,
                           real* gridForcesX, real* gridForcesY, real* gridForcesZ, 
                           uint* gridIndices, uint numberOfIndices, 
                           real* bladeCoordsX, real* bladeCoordsY, real* bladeCoordsZ, 
                           real* bladeForcesX, real* bladeForcesY,real* bladeForcesZ,
                           real* bladeRadii,
                           real radius,
                           uint nBlades, uint nBladeNodes,
                           real epsilon, real delta_x)
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

    real eta = 0.0f;

    real delta_x_cubed = pow(delta_x,3);

    for( uint blade=0; blade<nBlades; blade++)
    {    
        real last_r = 0.0f;
        real r = 0.0f;

        for( uint bladeNode=0; bladeNode<nBladeNodes; bladeNode++)
        {
            int node = bladeNode+blade*nBladeNodes;
            eta = calcGaussian3D(posX, posY, posZ, bladeCoordsX[node], bladeCoordsY[node], bladeCoordsZ[node], epsilon)*delta_x_cubed;
            r = bladeRadii[bladeNode];

            fXYZ_X += bladeForcesX[node]*(r-last_r)*eta;
            fXYZ_Y += bladeForcesY[node]*(r-last_r)*eta;
            fXYZ_Z += bladeForcesZ[node]*(r-last_r)*eta;

            last_r = r;
        }    

        fXYZ_X += bladeForcesX[nBladeNodes-1]*(radius-last_r)*eta;
        fXYZ_Y += bladeForcesY[nBladeNodes-1]*(radius-last_r)*eta;
        fXYZ_Z += bladeForcesZ[nBladeNodes-1]*(radius-last_r)*eta;
    }
    atomicAdd(&gridForcesX[gridIndex], fXYZ_X*invForceRatio);
    atomicAdd(&gridForcesY[gridIndex], fXYZ_Y*invForceRatio);
    atomicAdd(&gridForcesZ[gridIndex], fXYZ_Z*invForceRatio);
}


void ActuatorLine::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{
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
        para->getParD(this->level)->d0SP.f[0], para->getParD(this->level)->size_Mat_SP, para->getParD(this->level)->evenOrOdd,
        this->numberOfIndices,
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeVelocitiesXD, this->bladeVelocitiesYD, this->bladeVelocitiesZD,  
        this->bladeIndicesD, this->numberOfNodes);

    cudaManager->cudaCopyBladeVelocitiesDtoH(this);

    if(true)
    {
        this->calcForcesEllipticWing(para);
    }

    cudaManager->cudaCopyBladeForcesHtoD(this);

    vf::cuda::CudaGrid sphereGrid = vf::cuda::CudaGrid(numberOfThreads, this->numberOfIndices);

    applyBodyForces<<<sphereGrid.grid, sphereGrid.threads>>>(
        para->getParD(this->level)->coordX_SP, para->getParD(this->level)->coordY_SP, para->getParD(this->level)->coordZ_SP,        
        para->getParD(this->level)->forceX_SP, para->getParD(this->level)->forceY_SP, para->getParD(this->level)->forceZ_SP,        
        this->boundingSphereIndicesD, this->numberOfIndices,
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeForcesXD, this->bladeForcesYD, this->bladeForcesZD,
        this->bladeRadiiD,
        this->diameter*0.5f,  
        this->nBlades, this->nBladeNodes,
        this->epsilon, this->delta_x);

    real dazimuth = this->omega*this->delta_t;

    this->azimuth += dazimuth;
    this->rotateBlades(dazimuth);
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


void ActuatorLine::calcForcesEllipticWing(Parameter* para)
{
    real localAzimuth;
    uint node;
    real uXYZ_X, uXYZ_Y, uXYZ_Z;
    real uRTZ_X, uRTZ_Y, uRTZ_Z;
    real fXYZ_X, fXYZ_Y, fXYZ_Z;
    real fRTZ_X, fRTZ_Y, fRTZ_Z;
    real r;
    real u_rel, v_rel, u_rel_sq;
    real phi;
    real Cl = 1.f;
    real Cd = 0.f;
    real c0 = 1.f;

    real c, Cn, Ct;

    real forceRatio = this->density*pow(this->delta_x,4)*pow(this->delta_t,-2);

    for( uint blade=0; blade<this->nBlades; blade++)
    {
        localAzimuth = this->azimuth+2*blade*vf::lbm::constant::cPi/this->nBlades;
        for( uint bladeNode=0; bladeNode<this->nBladeNodes; bladeNode++)
        {
            node = bladeNode+blade*this->nBladeNodes;
            uXYZ_X = this->bladeVelocitiesXH[node]*para->getVelocityRatio();
            uXYZ_Y = this->bladeVelocitiesYH[node]*para->getVelocityRatio();
            uXYZ_Z = this->bladeVelocitiesZH[node]*para->getVelocityRatio();

            invRotateAboutX3D(localAzimuth, uXYZ_X, uXYZ_Y, uXYZ_Z, uRTZ_X, uRTZ_Y, uRTZ_Z);
            r = this->bladeRadiiH[bladeNode];

            u_rel = uRTZ_X;
            v_rel = uRTZ_Y+this->omega*r;
            u_rel_sq = u_rel*u_rel+v_rel*v_rel;
            phi = atan2(u_rel, v_rel);

            c = c0 * sqrt( 1.f- pow(4.f*r/this->diameter-1.f, 2.f) );
            Cn =   Cl*cos(phi)+Cd*sin(phi);
            Ct =  -Cl*sin(phi)+Cd*cos(phi);

            fRTZ_X = 0.5f*u_rel_sq*c*this->density*Cn;
            fRTZ_Y = 0.5f*u_rel_sq*c*this->density*Ct;
            fRTZ_Z = 0.0;

            rotateAboutX3D(localAzimuth, fRTZ_X, fRTZ_Y, fRTZ_Z, fXYZ_X, fXYZ_Y, fXYZ_Z);
        
            this->bladeForcesXH[node] = fXYZ_X/forceRatio;
            this->bladeForcesYH[node] = fXYZ_Y/forceRatio;
            this->bladeForcesZH[node] = fXYZ_Z/forceRatio;
        }
    }
}

void ActuatorLine::rotateBlades(real angle)
{
    for(uint node=0; node<this->nBladeNodes*this->nBlades; node++)
    {
        real oldCoordX = this->bladeCoordsXH[node];
        real oldCoordY = this->bladeCoordsYH[node];
        real oldCoordZ = this->bladeCoordsZH[node];

        real newCoordX, newCoordY, newCoordZ;
        rotateAboutX3D(angle, oldCoordX, oldCoordY, oldCoordZ, newCoordX, newCoordY, newCoordZ, this->turbinePosX, this->turbinePosY, this->turbinePosZ);
        
        this->bladeCoordsXH[node] = newCoordX;
        this->bladeCoordsYH[node] = newCoordY;
        this->bladeCoordsZH[node] = newCoordZ;
    }
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

    for( uint blade=0; blade<this->nBlades; blade++)
    {
        real localAzimuth = this->azimuth+(2*vf::lbm::constant::cPi/this->nBlades)*blade;
        for(uint node=0; node<this->nBladeNodes; node++)
        {
            real coordX, coordY, coordZ;
            real x,y,z;
            x = 0.f;
            y = 0.f;
            z = this->bladeRadiiH[node];
            rotateAboutX3D(localAzimuth, x, y, z, coordX, coordY, coordZ);
            this->bladeCoordsXH[node+this->nBladeNodes*blade] = coordX+this->turbinePosX;
            this->bladeCoordsYH[node+this->nBladeNodes*blade] = coordY+this->turbinePosY;
            this->bladeCoordsZH[node+this->nBladeNodes*blade] = coordZ+this->turbinePosZ;
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

    real* coordsX = para->getParH(this->level)->coordX_SP;
    real* coordsY = para->getParH(this->level)->coordY_SP;
    real* coordsZ = para->getParH(this->level)->coordZ_SP;

    for(uint node=0; node<this->numberOfNodes; node++)
    {
        this->bladeIndicesH[node] = findNearestCellBSW(1, coordsX, coordsY, coordsZ, 
                                                       this->bladeCoordsXH[node], this->bladeCoordsYH[node], this->bladeCoordsZH[node],
                                                       para->getParH(this->level)->neighborX_SP, para->getParH(this->level)->neighborY_SP, para->getParH(this->level)->neighborZ_SP,
                                                       para->getParH(this->level)->neighborWSB_SP);
        
    }
    cudaManager->cudaCopyBladeIndicesHtoD(this);
}

void ActuatorLine::initBoundingSphere(Parameter* para, CudaMemoryManager* cudaManager)
{
    // Actuator line exists only on 1 level
    std::vector<int> nodesInSphere;

    for (uint j = 1; j <= para->getParH(this->level)->size_Mat_SP; j++)
    {
        const real coordX = para->getParH(this->level)->coordX_SP[j];
        const real coordY = para->getParH(this->level)->coordY_SP[j];
        const real coordZ = para->getParH(this->level)->coordZ_SP[j];
        const real dist = sqrt(pow(coordX-this->turbinePosX,2)+pow(coordY-this->turbinePosY,2)+pow(coordZ-this->turbinePosZ,2));
        
        if(dist < 0.6*this->diameter) nodesInSphere.push_back(j);
    }

    this->numberOfIndices = uint(nodesInSphere.size());
    cudaManager->cudaAllocSphereIndices(this);
    std::copy(nodesInSphere.begin(), nodesInSphere.end(), this->boundingSphereIndicesH);
    cudaManager->cudaCopySphereIndicesHtoD(this);
}
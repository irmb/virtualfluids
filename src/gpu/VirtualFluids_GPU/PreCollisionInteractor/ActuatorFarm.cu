//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file ActuatorFarm.cu
//! \ingroup PreCollisionInteractor
//! \author Henrik Asmuth, Henry Korb
//======================================================================================
#include "ActuatorFarm.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "cuda/CudaGrid.h"
#include "VirtualFluids_GPU/GPU/GeometryUtils.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"

#include "Parameter/Parameter.h"
#include "Parameter/CudaStreamManager.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"
#include "lbm/constants/NumericConstants.h"
#include "logger/Logger.h"

using namespace vf::lbm::constant;


__host__ __device__ __inline__ uint calcNode(uint bladeNode, uint numberOfBladeNodes, uint blade, uint numberOfBlades, uint turbine, uint numberOfTurbines)
{

    return bladeNode+numberOfBladeNodes*(blade+numberOfBlades*turbine);
}

__host__ __device__ __inline__ void calcTurbineBladeAndBladeNode(uint node, uint& bladeNode, uint numberOfBladeNodes, uint& blade, uint numberOfBlades, uint& turbine, uint numberOfTurbines)
{
    turbine = node/(numberOfBladeNodes*numberOfBlades);
    uint x_off = turbine*numberOfBladeNodes*numberOfBlades;
    blade = (node - x_off)/numberOfBlades;
    uint y_off = numberOfBladeNodes*blade+x_off;
    bladeNode = (node - y_off);
}

__host__ __device__ __forceinline__ real distSqrd(real distX, real distY, real distZ)
{
    return distX*distX+distY*distY+distZ*distZ;
}

void swapArrays(real* &arr1, real* &arr2)
{
    real* tmp = arr1;
    arr1 = arr2;
    arr2 = tmp;
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
                                      uint numberOfTurbines, uint numberOfBlades, uint numberOfBladeNodes, 
                                      real* azimuths, real* yaws, real* omegas, 
                                      real* turbPosX, real* turbPosY, real* turbPosZ,
                                      uint* bladeIndices, real velocityRatio, real invDeltaX)
{

    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = vf::gpu::getNodeIndex();

    if(nodeIndex>=numberOfBladeNodes*numberOfBlades*numberOfTurbines) return;

    uint turbine, bladeNode, blade;

    calcTurbineBladeAndBladeNode(nodeIndex, bladeNode, numberOfBladeNodes, blade, numberOfBlades, turbine, numberOfTurbines);

    real bladeCoordX_BF = bladeCoordsX[nodeIndex];
    real bladeCoordY_BF = bladeCoordsY[nodeIndex];
    real bladeCoordZ_BF = bladeCoordsZ[nodeIndex];

    real bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF;

    real localAzimuth = azimuths[turbine]+blade*c2Pi/numberOfBlades;
    real yaw = yaws[turbine];


    rotateFromBladeToGlobal(bladeCoordX_BF, bladeCoordY_BF, bladeCoordZ_BF, 
                            bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF,
                            localAzimuth, yaw);

    bladeCoordX_GF += turbPosX[turbine];
    bladeCoordY_GF += turbPosY[turbine];
    bladeCoordZ_GF += turbPosZ[turbine];

    uint k, ke, kn, kt;
    uint kne, kte, ktn, ktne;

    k = findNearestCellBSW(bladeIndices[nodeIndex], 
                           gridCoordsX, gridCoordsY, gridCoordsZ, 
                           bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF, 
                           neighborsX, neighborsY, neighborsZ, neighborsWSB);
        
    bladeIndices[nodeIndex] = k;

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

    rotateFromGlobalToBlade(bladeVelX_BF, bladeVelY_BF, bladeVelZ_BF, 
                            bladeVelX_GF, bladeVelY_GF, bladeVelZ_GF, 
                            localAzimuth, yaw);

    bladeVelocitiesX[nodeIndex] = bladeVelX_BF;
    bladeVelocitiesY[nodeIndex] = bladeVelY_BF+omegas[turbine]*bladeCoordZ_BF;
    bladeVelocitiesZ[nodeIndex] = bladeVelZ_BF;
}


__global__ void applyBodyForces(real* gridCoordsX, real* gridCoordsY, real* gridCoordsZ,
                                real* gridForcesX, real* gridForcesY, real* gridForcesZ, 
                                real* bladeCoordsX, real* bladeCoordsY, real* bladeCoordsZ, 
                                real* bladeForcesX, real* bladeForcesY, real* bladeForcesZ,
                                const uint numberOfTurbines, const uint numberOfBlades, const uint numberOfBladeNodes,
                                real* azimuths, real* yaws, real* diameters,
                                real* turbPosX, real* turbPosY, real* turbPosZ,
                                uint* gridIndices, uint nIndices, 
                                const real invEpsilonSqrd, const real factorGaussian)
{

    const uint index = vf::gpu::getNodeIndex();

    if(index>=nIndices) return;


    uint gridIndex = gridIndices[index];

    real gridCoordX_GF = gridCoordsX[gridIndex];
    real gridCoordY_GF = gridCoordsY[gridIndex];
    real gridCoordZ_GF = gridCoordsZ[gridIndex];

    real gridForceX_RF = c0o1;
    real gridForceY_RF = c0o1;
    real gridForceZ_RF = c0o1;

    real dAzimuth = c2Pi/numberOfBlades;

    for(uint turbine = 0; turbine<numberOfTurbines; turbine++)
    {
        real radius = c1o2*diameters[turbine];
        real gridCoordX_RF = gridCoordX_GF - turbPosX[turbine];
        real gridCoordY_RF = gridCoordY_GF - turbPosY[turbine];
        real gridCoordZ_RF = gridCoordZ_GF - turbPosZ[turbine];

        if(distSqrd(gridCoordX_RF, gridCoordY_RF, gridCoordZ_RF)*invEpsilonSqrd > radius*radius*invEpsilonSqrd+c7o1)
            continue;

        real azimuth = azimuths[turbine];
        real yaw = yaws[turbine];

        for( uint blade=0; blade<numberOfBlades; blade++)
        { 
            real localAzimuth = azimuth+blade*dAzimuth;


            real gridCoordX_BF, gridCoordY_BF, gridCoordZ_BF;

            rotateFromGlobalToBlade(gridCoordX_BF, gridCoordY_BF, gridCoordZ_BF,
                                    gridCoordX_RF, gridCoordY_RF, gridCoordZ_RF,
                                    localAzimuth, yaw);
            
            uint node;
            uint nextNode = calcNode(0, numberOfBladeNodes, blade, numberOfBlades, turbine, numberOfTurbines);

            real last_z = c0o1;
            real current_z = c0o1;
            real next_z = bladeCoordsZ[nextNode];

            real x, y, dz, eta, forceX_RF, forceY_RF, forceZ_RF;

            for( uint bladeNode=0; bladeNode<numberOfBladeNodes-1; bladeNode++)
            {
                node = nextNode;
                nextNode = calcNode(bladeNode+1, numberOfBladeNodes, blade, numberOfBlades, turbine, numberOfTurbines);

                x = bladeCoordsX[node];
                y = bladeCoordsY[node];
                last_z = current_z;
                current_z = next_z;
                next_z = bladeCoordsZ[nextNode];

                dz = c1o2*(next_z-last_z);

                eta = dz*factorGaussian*exp(-distSqrd(x-gridCoordX_BF, y-gridCoordY_BF, current_z-gridCoordZ_BF)*invEpsilonSqrd);
                rotateFromBladeToGlobal(bladeForcesX[node], bladeForcesY[node], bladeForcesZ[node], 
                                        forceX_RF, forceY_RF, forceZ_RF, 
                                        localAzimuth, yaw);
                                        
                gridForceX_RF += forceX_RF*eta;
                gridForceY_RF += forceY_RF*eta;
                gridForceZ_RF += forceZ_RF*eta;
            }

            //Handle last node separately

            node = nextNode;

            x = bladeCoordsX[node];
            y = bladeCoordsY[node];
            last_z = current_z;
            current_z = next_z;

            dz = c1o2*(radius-last_z);

            eta = dz*factorGaussian*exp(-distSqrd(x-gridCoordX_BF, y-gridCoordY_BF, current_z-gridCoordZ_BF)*invEpsilonSqrd);

            rotateFromBladeToGlobal(bladeForcesX[node], bladeForcesY[node], bladeForcesZ[node], 
                                    forceX_RF, forceY_RF, forceZ_RF, 
                                    localAzimuth, yaw);
                
            gridForceX_RF += forceX_RF*eta;
            gridForceY_RF += forceY_RF*eta;
            gridForceZ_RF += forceZ_RF*eta;
        }
    }

    gridForcesX[gridIndex] += gridForceX_RF;
    gridForcesY[gridIndex] += gridForceY_RF;
    gridForcesZ[gridIndex] += gridForceZ_RF;
}

void ActuatorFarm::addTurbine(real posX, real posY, real posZ, real diameter, real omega, real azimuth, real yaw, std::vector<real> bladeRadii)
{
    preInitPosX.push_back(posX);
    preInitPosY.push_back(posY);
    preInitPosZ.push_back(posZ);
    preInitOmegas.push_back(omega);
    preInitAzimuths.push_back(azimuth);
    preInitYaws.push_back(yaw);
    preInitDiameters.push_back(diameter);
    preInitBladeRadii.push_back(bladeRadii);
}

void ActuatorFarm::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaMemoryManager)
{
    if(!para->getIsBodyForce()) throw std::runtime_error("try to allocate ActuatorFarm but BodyForce is not set in Parameter.");
    this->forceRatio = para->getForceRatio();
    this->initTurbineGeometries(cudaMemoryManager);
    this->initBladeCoords(cudaMemoryManager);    
    this->initBladeIndices(para, cudaMemoryManager);
    this->initBladeVelocities(cudaMemoryManager);
    this->initBladeForces(cudaMemoryManager);    
    this->initBoundingSpheres(para, cudaMemoryManager);
    this->streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::ActuatorFarm);
}

void ActuatorFarm::interact(Parameter* para, CudaMemoryManager* cudaMemoryManager, int level, unsigned int t)
{
    if (level != this->level) return;

    cudaStream_t stream = para->getStreamManager()->getStream(CudaStreamIndex::ActuatorFarm, this->streamIndex);

    if(useHostArrays) cudaMemoryManager->cudaCopyBladeCoordsHtoD(this);

    vf::cuda::CudaGrid bladeGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->numberOfNodes);

    interpolateVelocities<<< bladeGrid.grid, bladeGrid.threads, 0, stream >>>(
        para->getParD(this->level)->coordinateX, para->getParD(this->level)->coordinateY, para->getParD(this->level)->coordinateZ,        
        para->getParD(this->level)->neighborX, para->getParD(this->level)->neighborY, para->getParD(this->level)->neighborZ, para->getParD(this->level)->neighborInverse,
        para->getParD(this->level)->velocityX, para->getParD(this->level)->velocityY, para->getParD(this->level)->velocityZ,
        this->bladeCoordsXDCurrentTimestep, this->bladeCoordsYDCurrentTimestep, this->bladeCoordsZDCurrentTimestep,  
        this->bladeVelocitiesXDCurrentTimestep, this->bladeVelocitiesYDCurrentTimestep, this->bladeVelocitiesZDCurrentTimestep,  
        this->numberOfTurbines, this->numberOfBlades, this->numberOfBladeNodes,
        this->azimuthsD, this->yawsD, this->omegasD, 
        this->turbinePosXD, this->turbinePosYD, this->turbinePosZD,
        this->bladeIndicesD, para->getVelocityRatio(), this->invDeltaX);

    cudaStreamSynchronize(stream);
    if(useHostArrays) cudaMemoryManager->cudaCopyBladeVelocitiesDtoH(this);
    this->calcBladeForces();
    this->swapDeviceArrays();

    if(useHostArrays) cudaMemoryManager->cudaCopyBladeForcesHtoD(this);

    vf::cuda::CudaGrid sphereGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->numberOfIndices);

    applyBodyForces<<<sphereGrid.grid, sphereGrid.threads, 0, stream>>>(
        para->getParD(this->level)->coordinateX, para->getParD(this->level)->coordinateY, para->getParD(this->level)->coordinateZ,        
        para->getParD(this->level)->forceX_SP, para->getParD(this->level)->forceY_SP, para->getParD(this->level)->forceZ_SP,        
        this->bladeCoordsXDCurrentTimestep, this->bladeCoordsYDCurrentTimestep, this->bladeCoordsZDCurrentTimestep,  
        this->bladeForcesXDCurrentTimestep, this->bladeForcesYDCurrentTimestep, this->bladeForcesZDCurrentTimestep,
        this->numberOfTurbines, this->numberOfBlades, this->numberOfBladeNodes,
        this->azimuthsD, this->yawsD, this->diametersD,
        this->turbinePosXD, this->turbinePosYD, this->turbinePosZD,
        this->boundingSphereIndicesD, this->numberOfIndices,
        this->invEpsilonSqrd, this->factorGaussian);
    cudaMemoryManager->cudaCopyBladeOrientationsHtoD(this);
    cudaStreamSynchronize(stream);
}


void ActuatorFarm::free(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    cudaMemoryManager->cudaFreeBladeGeometries(this);
    cudaMemoryManager->cudaFreeBladeOrientations(this);
    cudaMemoryManager->cudaFreeBladeCoords(this);
    cudaMemoryManager->cudaFreeBladeVelocities(this);
    cudaMemoryManager->cudaFreeBladeForces(this);
    cudaMemoryManager->cudaFreeBladeIndices(this);
    cudaMemoryManager->cudaFreeSphereIndices(this);
}


void ActuatorFarm::calcForcesEllipticWing()
{
    real u_rel, v_rel, u_rel_sq;
    real phi;
    real Cl = c1o1;
    real Cd = c0o1;
    real c0 = 20*c1o10;
    real c, Cn, Ct;
    for(uint turbine=0; turbine<this->numberOfTurbines; turbine++)
    {
        real diameter = this->diametersH[turbine];
        for( uint blade=0; blade<this->numberOfBlades; blade++)
        { 
            for( uint bladeNode=0; bladeNode<this->numberOfBladeNodes; bladeNode++)
            {        
                uint node = calcNode(bladeNode, this->numberOfBladeNodes, blade, this->numberOfBlades, turbine, this->numberOfTurbines);

                u_rel = this->bladeVelocitiesXH[node];
                v_rel = this->bladeVelocitiesYH[node];
                u_rel_sq = u_rel*u_rel+v_rel*v_rel;
                phi = atan2(u_rel, v_rel);
                
                real tmp = c4o1*this->bladeRadiiH[bladeNode]/diameter-c1o1;
                c = c0 * sqrt( c1o1- tmp*tmp );
                Cn = Cl*cos(phi)+Cd*sin(phi);
                Ct = Cl*sin(phi)-Cd*cos(phi);
                real fx = c1o2*u_rel_sq*c*this->density*Cn;
                real fy = c1o2*u_rel_sq*c*this->density*Ct;
                this->bladeForcesXH[node] = -fx;
                this->bladeForcesYH[node] = -fy;
                this->bladeForcesZH[node] = c0o1;
                // printf("u %f v %f fx %f fy %f \n", u_rel, v_rel, fx, fy);
            }
        }
        azimuthsH[turbine] = azimuthsH[turbine]+deltaT*omegasH[turbine];
    }
}

void ActuatorFarm::calcBladeForces()
{
    this->calcForcesEllipticWing();
}

void ActuatorFarm::getTaggedFluidNodes(Parameter *para, GridProvider* gridProvider)
{
    std::vector<uint> indicesInSphere(this->boundingSphereIndicesH, this->boundingSphereIndicesH+this->numberOfIndices);
    gridProvider->tagFluidNodeIndices(indicesInSphere, CollisionTemplate::AllFeatures, this->level);
}   


void ActuatorFarm::initTurbineGeometries(CudaMemoryManager* cudaMemoryManager)
{
    this->numberOfTurbines = uint(this->preInitDiameters.size());
    this->numberOfNodes = numberOfTurbines*numberOfBladeNodes*numberOfBlades;

    cudaMemoryManager->cudaAllocBladeGeometries(this);
    cudaMemoryManager->cudaAllocBladeOrientations(this);

    for(uint turbine=0; turbine<this->numberOfTurbines; turbine++)
    {
        for(uint node=0; node<this->numberOfBladeNodes; node++)
        {
            this->bladeRadiiH[calcNode(node, numberOfBladeNodes, 0, 1, turbine, numberOfTurbines)] = this->preInitBladeRadii[turbine][node];
        }

    }
    std::copy(preInitPosX.begin(), preInitPosX.end(), turbinePosXH);
    std::copy(preInitPosY.begin(), preInitPosY.end(), turbinePosYH);
    std::copy(preInitPosZ.begin(), preInitPosZ.end(), turbinePosZH);
    std::copy(preInitDiameters.begin(), preInitDiameters.end(), diametersH);

    cudaMemoryManager->cudaCopyBladeGeometriesHtoD(this);
    std::copy(preInitAzimuths.begin(), preInitAzimuths.end(), this->azimuthsH);
    std::copy(preInitOmegas.begin(), preInitOmegas.end(), this->omegasH);
    std::copy(preInitYaws.begin(), preInitYaws.end(), this->yawsH);

    cudaMemoryManager->cudaCopyBladeOrientationsHtoD(this);
    this->factorGaussian = pow(this->epsilon*sqrt(cPi),-c3o1)/this->forceRatio;
}

void ActuatorFarm::initBladeCoords(CudaMemoryManager* cudaMemoryManager)
{   
    cudaMemoryManager->cudaAllocBladeCoords(this);

    for(uint turbine=0; turbine<numberOfTurbines; turbine++)
    {
        for(uint blade=0; blade<this->numberOfBlades; blade++)
        {
            for(uint bladeNode=0; bladeNode<this->numberOfBladeNodes; bladeNode++)
            {
                uint node = calcNode(bladeNode, this->numberOfBladeNodes, blade, this->numberOfBlades, turbine, this->numberOfTurbines);

                this->bladeCoordsXH[node] = c0o1;
                this->bladeCoordsYH[node] = c0o1;
                this->bladeCoordsZH[node] = this->bladeRadiiH[calcNode(bladeNode, numberOfBladeNodes, 0, 1, turbine, numberOfTurbines)];
            }
        }
    }
    cudaMemoryManager->cudaCopyBladeCoordsHtoD(this);
    swapArrays(this->bladeCoordsXDCurrentTimestep, this->bladeCoordsXDPreviousTimestep);
    swapArrays(this->bladeCoordsYDCurrentTimestep, this->bladeCoordsYDPreviousTimestep);
    swapArrays(this->bladeCoordsZDCurrentTimestep, this->bladeCoordsZDPreviousTimestep);
    cudaMemoryManager->cudaCopyBladeCoordsHtoD(this);
}

void ActuatorFarm::initBladeVelocities(CudaMemoryManager* cudaMemoryManager)
{   
    cudaMemoryManager->cudaAllocBladeVelocities(this);

    std::fill_n(this->bladeVelocitiesXH, this->numberOfNodes, c0o1);
    std::fill_n(this->bladeVelocitiesYH, this->numberOfNodes, c0o1);
    std::fill_n(this->bladeVelocitiesZH, this->numberOfNodes, c0o1);

    cudaMemoryManager->cudaCopyBladeVelocitiesHtoD(this);
    swapArrays(this->bladeVelocitiesXDCurrentTimestep, this->bladeVelocitiesXDPreviousTimestep);
    swapArrays(this->bladeVelocitiesYDCurrentTimestep, this->bladeVelocitiesYDPreviousTimestep);
    swapArrays(this->bladeVelocitiesZDCurrentTimestep, this->bladeVelocitiesZDPreviousTimestep);
    cudaMemoryManager->cudaCopyBladeVelocitiesHtoD(this);
}

void ActuatorFarm::initBladeForces(CudaMemoryManager* cudaMemoryManager)
{   
    cudaMemoryManager->cudaAllocBladeForces(this);

    std::fill_n(this->bladeForcesXH, this->numberOfNodes, c0o1);
    std::fill_n(this->bladeForcesYH, this->numberOfNodes, c0o1);
    std::fill_n(this->bladeForcesZH, this->numberOfNodes, c0o1);

    cudaMemoryManager->cudaCopyBladeForcesHtoD(this);
    swapArrays(this->bladeForcesXDCurrentTimestep, this->bladeForcesXDPreviousTimestep);
    swapArrays(this->bladeForcesYDCurrentTimestep, this->bladeForcesYDPreviousTimestep);
    swapArrays(this->bladeForcesZDCurrentTimestep, this->bladeForcesZDPreviousTimestep);
    cudaMemoryManager->cudaCopyBladeForcesHtoD(this);
}

void ActuatorFarm::initBladeIndices(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{   
    cudaMemoryManager->cudaAllocBladeIndices(this);

    std::fill_n(this->bladeIndicesH, this->numberOfNodes, 1);

    cudaMemoryManager->cudaCopyBladeIndicesHtoD(this);
}

void ActuatorFarm::initBoundingSpheres(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    std::vector<int> nodesInSpheres;

    for(uint turbine=0; turbine<this->numberOfTurbines; turbine++)
    {
        real sphereRadius = c1o2*this->diametersH[turbine]+c4o1*this->epsilon;

        real posX = this->turbinePosXH[turbine];
        real posY = this->turbinePosYH[turbine];
        real posZ = this->turbinePosZH[turbine];

        real sphereRadiusSqrd = sphereRadius*sphereRadius;
            
        uint minimumNumberOfNodesPerSphere = (uint)(c4o3*cPi*pow(sphereRadius-this->deltaX, c3o1)/pow(this->deltaX, c3o1));
        uint nodesInThisSphere = 0;

        for (size_t pos = 1; pos <= para->getParH(this->level)->numberOfNodes; pos++)
        {
            const real distX = para->getParH(this->level)->coordinateX[pos]-posX;
            const real distY = para->getParH(this->level)->coordinateY[pos]-posY;
            const real distZ = para->getParH(this->level)->coordinateZ[pos]-posZ;
            if(distSqrd(distX,distY,distZ) < sphereRadiusSqrd) 
            {
                nodesInSpheres.push_back((int)pos);
                nodesInThisSphere++;
            }
        }

        if(nodesInThisSphere<minimumNumberOfNodesPerSphere)
        {
            VF_LOG_CRITICAL("Found only {} nodes in bounding sphere of turbine no. {}, expected at least {}!", nodesInThisSphere, turbine, minimumNumberOfNodesPerSphere);
            throw std::runtime_error("ActuatorFarm::initBoundingSpheres: Turbine bounding sphere partially out of domain.");
        }
    }

    this->numberOfIndices = uint(nodesInSpheres.size());

    cudaMemoryManager->cudaAllocSphereIndices(this);
    std::copy(nodesInSpheres.begin(), nodesInSpheres.end(), this->boundingSphereIndicesH);
    cudaMemoryManager->cudaCopySphereIndicesHtoD(this);
}

void ActuatorFarm::setAllAzimuths(real* _azimuths)
{ 
    std::copy_n(_azimuths, this->numberOfTurbines, this->azimuthsH);
}

void ActuatorFarm::setAllOmegas(real* _omegas)
{ 
    std::copy_n(_omegas, this->numberOfTurbines, this->omegasH);
}

void ActuatorFarm::setAllYaws(real* _yaws)
{ 
    std::copy_n(_yaws, this->numberOfTurbines, this->yawsH);
}

void ActuatorFarm::setAllBladeCoords(real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ)
{ 
    std::copy_n(_bladeCoordsX, this->numberOfNodes, this->bladeCoordsXH);
    std::copy_n(_bladeCoordsY, this->numberOfNodes, this->bladeCoordsYH);
    std::copy_n(_bladeCoordsZ, this->numberOfNodes, this->bladeCoordsZH);
}

void ActuatorFarm::setAllBladeVelocities(real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ)
{ 
    std::copy_n(_bladeVelocitiesX, this->numberOfNodes, this->bladeVelocitiesXH);
    std::copy_n(_bladeVelocitiesY, this->numberOfNodes, this->bladeVelocitiesYH);
    std::copy_n(_bladeVelocitiesZ, this->numberOfNodes, this->bladeVelocitiesZH);
}

void ActuatorFarm::setAllBladeForces(real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ)
{ 
    std::copy_n(_bladeForcesX, this->numberOfNodes, this->bladeForcesXH);
    std::copy_n(_bladeForcesY, this->numberOfNodes, this->bladeForcesYH);
    std::copy_n(_bladeForcesZ, this->numberOfNodes, this->bladeForcesZH);

}void ActuatorFarm::setTurbineBladeCoords(uint turbine, real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ)
{ 
    std::copy_n(_bladeCoordsX, numberOfBladeNodes*numberOfBlades, &this->bladeCoordsXH[turbine*numberOfBladeNodes*numberOfBlades]);
    std::copy_n(_bladeCoordsY, numberOfBladeNodes*numberOfBlades, &this->bladeCoordsYH[turbine*numberOfBladeNodes*numberOfBlades]);
    std::copy_n(_bladeCoordsZ, numberOfBladeNodes*numberOfBlades, &this->bladeCoordsZH[turbine*numberOfBladeNodes*numberOfBlades]);
}

void ActuatorFarm::setTurbineBladeVelocities(uint turbine, real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ)
{ 
    std::copy_n(_bladeVelocitiesX, numberOfBladeNodes*numberOfBlades, &this->bladeVelocitiesXH[turbine*numberOfBladeNodes*numberOfBlades]);
    std::copy_n(_bladeVelocitiesY, numberOfBladeNodes*numberOfBlades, &this->bladeVelocitiesYH[turbine*numberOfBladeNodes*numberOfBlades]);
    std::copy_n(_bladeVelocitiesZ, numberOfBladeNodes*numberOfBlades, &this->bladeVelocitiesZH[turbine*numberOfBladeNodes*numberOfBlades]);
}

void ActuatorFarm::setTurbineBladeForces(uint turbine, real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ)
{ 
    std::copy_n(_bladeForcesX, numberOfBladeNodes*numberOfBlades, &this->bladeForcesXH[turbine*numberOfBladeNodes*numberOfBlades]);
    std::copy_n(_bladeForcesY, numberOfBladeNodes*numberOfBlades, &this->bladeForcesYH[turbine*numberOfBladeNodes*numberOfBlades]);
    std::copy_n(_bladeForcesZ, numberOfBladeNodes*numberOfBlades, &this->bladeForcesZH[turbine*numberOfBladeNodes*numberOfBlades]);
}

void ActuatorFarm::swapDeviceArrays()
{
    swapArrays(this->bladeCoordsXDPreviousTimestep, this->bladeCoordsXDCurrentTimestep);
    swapArrays(this->bladeCoordsYDPreviousTimestep, this->bladeCoordsYDCurrentTimestep);
    swapArrays(this->bladeCoordsZDPreviousTimestep, this->bladeCoordsZDCurrentTimestep);

    swapArrays(this->bladeVelocitiesXDPreviousTimestep, this->bladeVelocitiesXDCurrentTimestep);
    swapArrays(this->bladeVelocitiesYDPreviousTimestep, this->bladeVelocitiesYDCurrentTimestep);
    swapArrays(this->bladeVelocitiesZDPreviousTimestep, this->bladeVelocitiesZDCurrentTimestep);

    swapArrays(this->bladeForcesXDPreviousTimestep, this->bladeForcesXDCurrentTimestep);
    swapArrays(this->bladeForcesYDPreviousTimestep, this->bladeForcesYDCurrentTimestep);
    swapArrays(this->bladeForcesZDPreviousTimestep, this->bladeForcesZDCurrentTimestep);
}
#include "ActuatorFarm.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cuda/CudaGrid.h>
#include "VirtualFluids_GPU/GPU/GeometryUtils.h"
#include "VirtualFluids_GPU/Kernel/Utilities/DistributionHelper.cuh"

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"
#include "lbm/constants/NumericConstants.h"

using namespace vf::lbm::constant;


__host__ __device__ __inline__ uint calcNode(uint bladeNode, uint numberOfBladeNodes, uint blade, uint numberOfBlades, uint turbine, uint numberOfTurbines)
{

    return bladeNode+numberOfBladeNodes*(blade+turbine*numberOfBlades);
}

__host__ __device__ __inline__ void calcTurbineBladeAndBladeNode(uint node, uint& bladeNode, uint numberOfBladeNodes, uint& blade, uint numberOfBlades, uint& turbine, uint numberOfTurbines)
{
    turbine = node/(numberOfBladeNodes*numberOfBlades);
    uint x_off = turbine*numberOfBladeNodes*numberOfBlades;
    blade = (node - x_off)/numberOfBlades;
    uint y_off = numberOfBladeNodes*blade+x_off;
    bladeNode = (node - y_off)/numberOfBladeNodes;
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
                                      uint numberOfTurbines, uint numberOfBlades, uint numberOfBladeNodes, 
                                      real* azimuths, real* yaws, real* omegas, 
                                      real* turbPosX, real* turbPosY, real* turbPosZ,
                                      uint* bladeIndices, real velocityRatio, real invDeltaX)
{

    const uint node =  vf::gpu::getNodeIndex();

    if(node>=numberOfBladeNodes*numberOfBlades*numberOfTurbines) return;

    uint turbine, bladeNode, blade;

    calcTurbineBladeAndBladeNode(node, bladeNode, numberOfBladeNodes, blade, numberOfBlades, turbine, numberOfTurbines);

    real bladeCoordX_BF = bladeCoordsX[node];
    real bladeCoordY_BF = bladeCoordsY[node];
    real bladeCoordZ_BF = bladeCoordsZ[node];

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

    k = findNearestCellBSW(bladeIndices[node], 
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


    rotateFromGlobalToBlade(bladeVelX_BF, bladeVelY_BF, bladeVelZ_BF, 
                            bladeVelX_GF, bladeVelY_GF, bladeVelZ_GF, 
                            localAzimuth, yaw);

    bladeVelocitiesX[node] = bladeVelX_BF;
    bladeVelocitiesY[node] = bladeVelY_BF+omegas[turbine]*bladeCoordZ_BF;
    bladeVelocitiesZ[node] = bladeVelZ_BF;
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
            uint nextNode = calcNode(0, numberOfBladeNodes, blade, numberOfBlades, turbine, numberOfTurbines);;

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

            // if(eta>1e-4) printf("gridForcs %f % f% f",gridForceX_RF, gridForceY_RF, gridForceZ_RF );
        }
    }
    // if(index==nIndices/2) printf("gridForcs %f % f% f",gridForceX_RF, gridForceY_RF, gridForceZ_RF );
    // if(index==nIndices/4) printf("gridForcs %f % f% f",gridForceX_RF, gridForceY_RF, gridForceZ_RF );
    // if(index==3*nIndices/4) printf("gridForcs %f % f% f",gridForceX_RF, gridForceY_RF, gridForceZ_RF );
    // if(index==nIndices/2+10) printf("gridForcs %f % f% f",gridForceX_RF, gridForceY_RF, gridForceZ_RF );
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
    this->initBoundingSphere(para, cudaMemoryManager);    
}


void ActuatorFarm::interact(Parameter* para, CudaMemoryManager* cudaMemoryManager, int level, unsigned int t)
{
    if (level != this->level) return;

    if(useHostArrays) cudaMemoryManager->cudaCopyBladeCoordsHtoD(this);

    vf::cuda::CudaGrid bladeGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->numberOfNodes);

    interpolateVelocities<<< bladeGrid.grid, bladeGrid.threads >>>(
        para->getParD(this->level)->coordinateX, para->getParD(this->level)->coordinateY, para->getParD(this->level)->coordinateZ,        
        para->getParD(this->level)->neighborX, para->getParD(this->level)->neighborY, para->getParD(this->level)->neighborZ, para->getParD(this->level)->neighborInverse,
        para->getParD(this->level)->velocityX, para->getParD(this->level)->velocityY, para->getParD(this->level)->velocityZ,
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeVelocitiesXD, this->bladeVelocitiesYD, this->bladeVelocitiesZD,  
        this->numberOfTurbines, this->numberOfBlades, this->numberOfBladeNodes,
        this->azimuthsD, this->yawsD, this->omegasD, 
        this->turbinePosXD, this->turbinePosYD, this->turbinePosZD,
        this->bladeIndicesD, para->getVelocityRatio(), this->invDeltaX);

    if(useHostArrays) cudaMemoryManager->cudaCopyBladeVelocitiesDtoH(this);

    this->calcBladeForces();

    if(useHostArrays) cudaMemoryManager->cudaCopyBladeForcesHtoD(this);

    vf::cuda::CudaGrid sphereGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->numberOfIndices);

    applyBodyForces<<<sphereGrid.grid, sphereGrid.threads>>>(
        para->getParD(this->level)->coordinateX, para->getParD(this->level)->coordinateY, para->getParD(this->level)->coordinateZ,        
        para->getParD(this->level)->forceX_SP, para->getParD(this->level)->forceY_SP, para->getParD(this->level)->forceZ_SP,        
        this->bladeCoordsXD, this->bladeCoordsYD, this->bladeCoordsZD,  
        this->bladeForcesXD, this->bladeForcesYD, this->bladeForcesZD,
        this->numberOfTurbines, this->numberOfBlades, this->numberOfBladeNodes,
        this->azimuthsD, this->yawsD, this->diametersD,
        this->turbinePosXD, this->turbinePosYD, this->turbinePosZD,
        this->boundingSphereIndicesD, this->numberOfIndices,
        this->invEpsilonSqrd, this->factorGaussian);

    for(uint turbine=0; turbine<this->numberOfTurbines; turbine++)
        this->azimuthsH[turbine] = fmod(this->azimuthsH[turbine]+this->omegasH[turbine]*this->deltaT, c2Pi);
    cudaMemoryManager->cudaCopyBladeOrientationsHtoD(this);    
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
    real c0 = c1o1;

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
            
                this->bladeForcesXH[node] = -c1o2*u_rel_sq*c*this->density*Cn;
                this->bladeForcesYH[node] = -c1o2*u_rel_sq*c*this->density*Ct;
                this->bladeForcesZH[node] = c0o1;
            }
        }
    }
}

void ActuatorFarm::calcBladeForces()
{
    // this->calcForcesEllipticWing();
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
}

void ActuatorFarm::initBladeVelocities(CudaMemoryManager* cudaMemoryManager)
{   
    cudaMemoryManager->cudaAllocBladeVelocities(this);

    std::fill_n(this->bladeVelocitiesXH, this->numberOfNodes, c0o1);
    std::fill_n(this->bladeVelocitiesYH, this->numberOfNodes, c0o1);
    std::fill_n(this->bladeVelocitiesZH, this->numberOfNodes, c0o1);

    cudaMemoryManager->cudaCopyBladeVelocitiesHtoD(this);
}

void ActuatorFarm::initBladeForces(CudaMemoryManager* cudaMemoryManager)
{   
    cudaMemoryManager->cudaAllocBladeForces(this);

    std::fill_n(this->bladeForcesXH, this->numberOfNodes, c0o1);
    std::fill_n(this->bladeForcesYH, this->numberOfNodes, c0o1);
    std::fill_n(this->bladeForcesZH, this->numberOfNodes, c0o1);

    cudaMemoryManager->cudaCopyBladeForcesHtoD(this);
}

void ActuatorFarm::initBladeIndices(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{   
    cudaMemoryManager->cudaAllocBladeIndices(this);

    std::fill_n(this->bladeIndicesH, this->numberOfNodes, 1);

    cudaMemoryManager->cudaCopyBladeIndicesHtoD(this);
}

void ActuatorFarm::initBoundingSphere(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    std::vector<int> nodesInSphere;

    for(uint turbine=0; turbine<this->numberOfTurbines; turbine++)
    {
        real sphereRadius = c1o2*this->diametersH[turbine]+c4o1*this->epsilon;

        real posX = this->turbinePosXH[turbine];
        real posY = this->turbinePosYH[turbine];
        real posZ = this->turbinePosZH[turbine];
            
        real sphereRadiusSqrd = sphereRadius*sphereRadius;

        for (uint j = 1; j <= para->getParH(this->level)->numberOfNodes; j++)
        {
            const real distX = para->getParH(this->level)->coordinateX[j]-posX;
            const real distY = para->getParH(this->level)->coordinateY[j]-posY;
            const real distZ = para->getParH(this->level)->coordinateZ[j]-posZ;
            if(distSqrd(distX,distY,distZ) < sphereRadiusSqrd) nodesInSphere.push_back(j);
        }
    }

    this->numberOfIndices = uint(nodesInSphere.size());

    cudaMemoryManager->cudaAllocSphereIndices(this);
    std::copy(nodesInSphere.begin(), nodesInSphere.end(), this->boundingSphereIndicesH);
    cudaMemoryManager->cudaCopySphereIndicesHtoD(this);
}

void ActuatorFarm::setAllBladeCoords(real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ)
{ 
    std::copy_n(this->bladeCoordsXH, this->numberOfNodes, _bladeCoordsX);
    std::copy_n(this->bladeCoordsYH, this->numberOfNodes, _bladeCoordsY);
    std::copy_n(this->bladeCoordsZH, this->numberOfNodes, _bladeCoordsZ);

}

void ActuatorFarm::setAllBladeVelocities(real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ)
{ 
    std::copy_n(this->bladeVelocitiesXH, this->numberOfNodes, _bladeVelocitiesX);
    std::copy_n(this->bladeVelocitiesYH, this->numberOfNodes, _bladeVelocitiesY);
    std::copy_n(this->bladeVelocitiesZH, this->numberOfNodes, _bladeVelocitiesZ);

}

void ActuatorFarm::setAllBladeForces(real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ)
{ 
    std::copy_n(this->bladeForcesXH, this->numberOfNodes, _bladeForcesX);
    std::copy_n(this->bladeForcesYH, this->numberOfNodes, _bladeForcesY);
    std::copy_n(this->bladeForcesZH, this->numberOfNodes, _bladeForcesZ);

}void ActuatorFarm::setTurbineBladeCoords(uint turbine, real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ)
{ 
    std::copy_n(&this->bladeCoordsXH[turbine*numberOfBladeNodes*numberOfBlades], numberOfBladeNodes*numberOfBlades, _bladeCoordsX);
    std::copy_n(&this->bladeCoordsYH[turbine*numberOfBladeNodes*numberOfBlades], numberOfBladeNodes*numberOfBlades, _bladeCoordsY);
    std::copy_n(&this->bladeCoordsZH[turbine*numberOfBladeNodes*numberOfBlades], numberOfBladeNodes*numberOfBlades, _bladeCoordsZ);
}

void ActuatorFarm::setTurbineBladeVelocities(uint turbine, real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ)
{ 
    std::copy_n(&this->bladeVelocitiesXH[turbine*numberOfBladeNodes*numberOfBlades], numberOfBladeNodes*numberOfBlades, _bladeVelocitiesX);
    std::copy_n(&this->bladeVelocitiesYH[turbine*numberOfBladeNodes*numberOfBlades], numberOfBladeNodes*numberOfBlades, _bladeVelocitiesY);
    std::copy_n(&this->bladeVelocitiesZH[turbine*numberOfBladeNodes*numberOfBlades], numberOfBladeNodes*numberOfBlades, _bladeVelocitiesZ);
}

void ActuatorFarm::setTurbineBladeForces(uint turbine, real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ)
{ 
    std::copy_n(&this->bladeForcesXH[turbine*numberOfBladeNodes*numberOfBlades], numberOfBladeNodes*numberOfBlades, _bladeForcesX);
    std::copy_n(&this->bladeForcesYH[turbine*numberOfBladeNodes*numberOfBlades], numberOfBladeNodes*numberOfBlades, _bladeForcesY);
    std::copy_n(&this->bladeForcesZH[turbine*numberOfBladeNodes*numberOfBlades], numberOfBladeNodes*numberOfBlades, _bladeForcesZ);
}

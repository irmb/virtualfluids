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
#include "ActuatorFarmInlines.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <basics/constants/NumericConstants.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <cuda_helper/CudaGrid.h>
#include <logger/Logger.h>

#include "DataStructureInitializer/GridProvider.h"
#include "Cuda/CudaMemoryManager.h"
#include "Utilities/GeometryUtils.h"
#include "Utilities/KernelUtilities.h"
#include "Cuda/CudaStreamManager.h"
#include "Parameter/Parameter.h"

using namespace vf::basics::constant;

struct GridData
{
    const uint* indices;
    const uint nIndices;
    const real *coordsX, *coordsY, *coordsZ;
    const uint *neighborsX, *neighborsY, *neighborsZ, *neighborsWSB;
    const real *vx, *vy, *vz;
    real *fx, *fy, *fz;
    const real inverseDeltaX, velocityRatio;
};

struct TurbineData
{
    const real *posX, *posY, *posZ;
    const uint numberOfTurbines;
    const real smearingWidth, factorGaussian;
};

struct ComponentData
{
    const real referenceLength;
    const uint numberOfNodesPerTurbine;
    const real *coordsX, *coordsY, *coordsZ;
    real *velocitiesX, *velocitiesY, *velocitiesZ;
    const real *forcesX, *forcesY, *forcesZ;
    uint* gridIndices;
};

__global__ void interpolateVelocities(const GridData gridData, const TurbineData turbineData, ComponentData componentData)
{
    const unsigned nodeIndex = vf::gpu::getNodeIndex();

    if (nodeIndex >= componentData.numberOfNodesPerTurbine * turbineData.numberOfTurbines)
        return;

    const real coordX = componentData.coordsX[nodeIndex];
    const real coordY = componentData.coordsY[nodeIndex];
    const real coordZ = componentData.coordsZ[nodeIndex];

    uint k, ke, kn, kt;
    uint kne, kte, ktn, ktne;

    k = findNearestCellBSW(componentData.gridIndices[nodeIndex],
                           gridData.coordsX, gridData.coordsY, gridData.coordsZ,
                           coordX, coordY, coordZ,
                           gridData.neighborsX, gridData.neighborsY, gridData.neighborsZ, gridData.neighborsWSB);

    componentData.gridIndices[nodeIndex] = k;

    getNeighborIndicesOfBSW(k, ke, kn, kt, kne, kte, ktn, ktne, gridData.neighborsX, gridData.neighborsY,
                            gridData.neighborsZ);

    real dW, dE, dN, dS, dT, dB;

    const real distX = gridData.inverseDeltaX * (coordX - gridData.coordsX[k]);
    const real distY = gridData.inverseDeltaX * (coordY - gridData.coordsY[k]);
    const real distZ = gridData.inverseDeltaX * (coordZ - gridData.coordsZ[k]);

    getInterpolationWeights(dW, dE, dN, dS, dT, dB, distX, distY, distZ);

    componentData.velocitiesX[nodeIndex] =
        trilinearInterpolation(dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, gridData.vx) *
        gridData.velocityRatio;
    componentData.velocitiesY[nodeIndex] =
        trilinearInterpolation(dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, gridData.vy) *
        gridData.velocityRatio;
    componentData.velocitiesZ[nodeIndex] =
        trilinearInterpolation(dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, gridData.vz) *
        gridData.velocityRatio;
}

__global__ void applyBodyForces(GridData gridData, const TurbineData turbineData, const ComponentData componentData)
{

    const uint index = vf::gpu::getNodeIndex();

    if (index >= gridData.nIndices)
        return;

    const uint gridIndex = gridData.indices[index];
    const real gridCoordX = gridData.coordsX[gridIndex];
    const real gridCoordY = gridData.coordsY[gridIndex];
    const real gridCoordZ = gridData.coordsZ[gridIndex];

    real gridForceX = c0o1;
    real gridForceY = c0o1;
    real gridForceZ = c0o1;

    for (uint turbine = 0; turbine < turbineData.numberOfTurbines; turbine++) {
        const real distToHubX = gridCoordX - turbineData.posX[turbine];
        const real distToHubY = gridCoordY - turbineData.posY[turbine];
        const real distToHubZ = gridCoordZ - turbineData.posZ[turbine];

        if (!inBoundingSphere(distToHubX, distToHubY, distToHubZ, componentData.referenceLength, turbineData.smearingWidth))
            continue;

        for (uint turbineNode = 0; turbineNode < componentData.numberOfNodesPerTurbine; turbineNode++) {
            const uint node = turbine * componentData.numberOfNodesPerTurbine + turbineNode;

            const real distX = componentData.coordsX[node] - gridCoordX;
            const real distY = componentData.coordsY[node] - gridCoordY;
            const real distZ = componentData.coordsZ[node] - gridCoordZ;

            const real eta = gaussianSmearing(distX, distY, distZ, turbineData.smearingWidth, turbineData.factorGaussian);
            gridForceX += componentData.forcesX[node] * eta;
            gridForceY += componentData.forcesY[node] * eta;
            gridForceZ += componentData.forcesZ[node] * eta;
        }
    }
    gridData.fx[gridIndex] += gridForceX;
    gridData.fy[gridIndex] += gridForceY;
    gridData.fz[gridIndex] += gridForceZ;
}

void ActuatorFarm::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{
    if (!para->getIsBodyForce())
        throw std::runtime_error("try to allocate ActuatorFarm but BodyForce is not set in Parameter.");
    this->forceRatio = para->getForceRatio();
    this->initTurbineGeometries(cudaManager);
    this->initBladeCoords(cudaManager);
    this->initBladeIndices(cudaManager);
    this->initBladeVelocities(cudaManager);
    this->initBladeForces(cudaManager);
    this->initBoundingSpheres(para, cudaManager);
    this->streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::ActuatorFarm);
}

void ActuatorFarm::interact(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t)
{
    if (level != this->level)
        return;

    cudaStream_t stream = para->getStreamManager()->getStream(CudaStreamIndex::ActuatorFarm, this->streamIndex);

    if (useHostArrays)
        cudaManager->cudaCopyBladeCoordsHtoD(this);

    if (this->writeOutput && ((t - this->tStartOut) % this->tOut == 0)) {
        if (!useHostArrays) {
            cudaManager->cudaCopyBladeCoordsDtoH(this);
            cudaManager->cudaCopyBladeVelocitiesDtoH(this);
            cudaManager->cudaCopyBladeForcesDtoH(this);
        }
        this->write(this->getFilename(para, t));
    }

    const GridData gridData {
        this->boundingSphereIndicesD, this->numberOfIndices,
        para->getParD(this->level)->coordinateX, para->getParD(this->level)->coordinateY, para->getParD(this->level)->coordinateZ,
        para->getParD(this->level)->neighborX, para->getParD(this->level)->neighborY, para->getParD(this->level)->neighborZ, para->getParD(this->level)->neighborInverse,
        para->getParD(this->level)->velocityX, para->getParD(this->level)->velocityY, para->getParD(this->level)->velocityZ,
        para->getParD(this->level)->forceX_SP,para->getParD(this->level)->forceY_SP,para->getParD(this->level)->forceZ_SP,
        this->invDeltaX, para->getVelocityRatio()};

    const TurbineData turbineData {
        this->turbinePosXD, this->turbinePosYD, this->turbinePosZD,
        this->numberOfTurbines,
        this->smearingWidth, this->factorGaussian};
    
    const ComponentData bladeData {
        this->diameter, this->numberOfNodesPerTurbine,
        this->bladeCoordsXDCurrentTimestep, this->bladeCoordsYDCurrentTimestep, this->bladeCoordsZDCurrentTimestep, 
        this->bladeVelocitiesXDCurrentTimestep, this->bladeVelocitiesYDCurrentTimestep, this->bladeVelocitiesZDCurrentTimestep, 
        this->bladeForcesXDCurrentTimestep, this->bladeForcesYDCurrentTimestep, this->bladeForcesZDCurrentTimestep,
        this->bladeIndicesD};

    vf::cuda::CudaGrid bladeGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->numberOfGridNodes);
    interpolateVelocities<<<bladeGrid.grid, bladeGrid.threads, 0, stream>>>(gridData, turbineData, bladeData);
    cudaStreamSynchronize(stream);

    if (useHostArrays)
        cudaManager->cudaCopyBladeVelocitiesDtoH(this);

    this->updateForcesAndCoordinates();
    this->swapDeviceArrays();

    if (useHostArrays)
        cudaManager->cudaCopyBladeForcesHtoD(this);

    vf::cuda::CudaGrid sphereGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->numberOfIndices);

    applyBodyForces<<<sphereGrid.grid, sphereGrid.threads, 0, stream>>>(gridData, turbineData, bladeData);
    cudaStreamSynchronize(stream);
}

void ActuatorFarm::free(Parameter* para, CudaMemoryManager* cudaManager)
{
    cudaManager->cudaFreeBladeGeometries(this);
    cudaManager->cudaFreeBladeCoords(this);
    cudaManager->cudaFreeBladeVelocities(this);
    cudaManager->cudaFreeBladeForces(this);
    cudaManager->cudaFreeBladeIndices(this);
    cudaManager->cudaFreeSphereIndices(this);
}

void ActuatorFarm::getTaggedFluidNodes(Parameter* para, GridProvider* gridProvider)
{
    std::vector<uint> indicesInSphere(this->boundingSphereIndicesH, this->boundingSphereIndicesH + this->numberOfIndices);
    gridProvider->tagFluidNodeIndices(indicesInSphere, CollisionTemplate::AllFeatures, this->level);
}

void ActuatorFarm::initTurbineGeometries(CudaMemoryManager* cudaManager)
{
    this->numberOfGridNodes = this->numberOfTurbines * this->numberOfNodesPerTurbine;

    cudaManager->cudaAllocBladeGeometries(this);

    std::copy(initialTurbinePositionsX.begin(), initialTurbinePositionsX.end(), turbinePosXH);
    std::copy(initialTurbinePositionsY.begin(), initialTurbinePositionsY.end(), turbinePosYH);
    std::copy(initialTurbinePositionsZ.begin(), initialTurbinePositionsZ.end(), turbinePosZH);

    cudaManager->cudaCopyBladeGeometriesHtoD(this);

    this->factorGaussian = pow(this->smearingWidth * sqrt(cPi), -c3o1) / this->forceRatio;
}

void ActuatorFarm::initBladeCoords(CudaMemoryManager* cudaManager)
{
    cudaManager->cudaAllocBladeCoords(this);

    for (uint turbine = 0; turbine < this->numberOfTurbines; turbine++) {
        for (uint blade = 0; blade < ActuatorFarm::numberOfBlades; blade++) {
            const real local_azimuth = this->azimuths[turbine] + blade * c2Pi / ActuatorFarm::numberOfBlades;

            for (uint bladeNode = 0; bladeNode < this->numberOfNodesPerBlade; bladeNode++) {
                const uint node = calcNodeIndexInBladeArrays({ turbine, blade, bladeNode }, this->numberOfNodesPerBlade,
                                                             ActuatorFarm::numberOfBlades);

                real x, y, z;
                rotateFromBladeToGlobal(c0o1, c0o1, this->bladeRadii[bladeNode], x, y, z, local_azimuth);
                bladeCoordsXH[node] = x + this->turbinePosXH[turbine];
                bladeCoordsYH[node] = y + this->turbinePosYH[turbine];
                bladeCoordsZH[node] = z + this->turbinePosZH[turbine];
            }
        }
    }
    cudaManager->cudaCopyBladeCoordsHtoD(this);
    swapArrays(this->bladeCoordsXDCurrentTimestep, this->bladeCoordsXDPreviousTimestep);
    swapArrays(this->bladeCoordsYDCurrentTimestep, this->bladeCoordsYDPreviousTimestep);
    swapArrays(this->bladeCoordsZDCurrentTimestep, this->bladeCoordsZDPreviousTimestep);
    cudaManager->cudaCopyBladeCoordsHtoD(this);
}

void ActuatorFarm::initBladeVelocities(CudaMemoryManager* cudaManager)
{
    cudaManager->cudaAllocBladeVelocities(this);

    std::fill_n(this->bladeVelocitiesXH, this->numberOfGridNodes, c0o1);
    std::fill_n(this->bladeVelocitiesYH, this->numberOfGridNodes, c0o1);
    std::fill_n(this->bladeVelocitiesZH, this->numberOfGridNodes, c0o1);

    cudaManager->cudaCopyBladeVelocitiesHtoD(this);
    swapArrays(this->bladeVelocitiesXDCurrentTimestep, this->bladeVelocitiesXDPreviousTimestep);
    swapArrays(this->bladeVelocitiesYDCurrentTimestep, this->bladeVelocitiesYDPreviousTimestep);
    swapArrays(this->bladeVelocitiesZDCurrentTimestep, this->bladeVelocitiesZDPreviousTimestep);
    cudaManager->cudaCopyBladeVelocitiesHtoD(this);
}

void ActuatorFarm::initBladeForces(CudaMemoryManager* cudaManager)
{
    cudaManager->cudaAllocBladeForces(this);

    std::fill_n(this->bladeForcesXH, this->numberOfGridNodes, c0o1);
    std::fill_n(this->bladeForcesYH, this->numberOfGridNodes, c0o1);
    std::fill_n(this->bladeForcesZH, this->numberOfGridNodes, c0o1);

    cudaManager->cudaCopyBladeForcesHtoD(this);
    swapArrays(this->bladeForcesXDCurrentTimestep, this->bladeForcesXDPreviousTimestep);
    swapArrays(this->bladeForcesYDCurrentTimestep, this->bladeForcesYDPreviousTimestep);
    swapArrays(this->bladeForcesZDCurrentTimestep, this->bladeForcesZDPreviousTimestep);
    cudaManager->cudaCopyBladeForcesHtoD(this);
}

void ActuatorFarm::initBladeIndices(CudaMemoryManager* cudaManager)
{
    cudaManager->cudaAllocBladeIndices(this);

    std::fill_n(this->bladeIndicesH, this->numberOfGridNodes, 1);

    cudaManager->cudaCopyBladeIndicesHtoD(this);
}

void ActuatorFarm::initBoundingSpheres(Parameter* para, CudaMemoryManager* cudaManager)
{
    std::vector<int> nodesInSpheres;
    const real sphereRadius = getBoundingSphereRadius(this->diameter, this->smearingWidth);
    const real sphereRadiusSqrd = sphereRadius * sphereRadius;
    const uint minimumNumberOfNodesPerSphere =
        (uint)(c4o3 * cPi * pow(sphereRadius - this->deltaX, c3o1) / pow(this->deltaX, c3o1));

    for (uint turbine = 0; turbine < this->numberOfTurbines; turbine++) {

        const real posX = this->turbinePosXH[turbine];
        const real posY = this->turbinePosYH[turbine];
        const real posZ = this->turbinePosZH[turbine];

        uint nodesInThisSphere = 0;

        for (size_t pos = 1; pos <= para->getParH(this->level)->numberOfNodes; pos++) {
            const real distX = para->getParH(this->level)->coordinateX[pos] - posX;
            const real distY = para->getParH(this->level)->coordinateY[pos] - posY;
            const real distZ = para->getParH(this->level)->coordinateZ[pos] - posZ;
            if (distSqrd(distX, distY, distZ) < sphereRadiusSqrd) {
                nodesInSpheres.push_back((int)pos);
                nodesInThisSphere++;
            }
        }

        if (nodesInThisSphere < minimumNumberOfNodesPerSphere) {
            VF_LOG_CRITICAL("Found only {} nodes in bounding sphere of turbine no. {}, expected at least {}!",
                            nodesInThisSphere, turbine, minimumNumberOfNodesPerSphere);
            throw std::runtime_error("ActuatorFarm::initBoundingSpheres: Turbine bounding sphere partially out of domain.");
        }
    }

    this->numberOfIndices = uint(nodesInSpheres.size());

    cudaManager->cudaAllocSphereIndices(this);
    std::copy(nodesInSpheres.begin(), nodesInSpheres.end(), this->boundingSphereIndicesH);
    cudaManager->cudaCopySphereIndicesHtoD(this);
}

void ActuatorFarm::setAllBladeCoords(const real* _bladeCoordsX, const real* _bladeCoordsY, const real* _bladeCoordsZ)
{
    std::copy_n(_bladeCoordsX, this->numberOfGridNodes, this->bladeCoordsXH);
    std::copy_n(_bladeCoordsY, this->numberOfGridNodes, this->bladeCoordsYH);
    std::copy_n(_bladeCoordsZ, this->numberOfGridNodes, this->bladeCoordsZH);
}

void ActuatorFarm::setAllBladeVelocities(const real* _bladeVelocitiesX, const real* _bladeVelocitiesY,
                                         const real* _bladeVelocitiesZ)
{
    std::copy_n(_bladeVelocitiesX, this->numberOfGridNodes, this->bladeVelocitiesXH);
    std::copy_n(_bladeVelocitiesY, this->numberOfGridNodes, this->bladeVelocitiesYH);
    std::copy_n(_bladeVelocitiesZ, this->numberOfGridNodes, this->bladeVelocitiesZH);
}

void ActuatorFarm::setAllBladeForces(const real* _bladeForcesX, const real* _bladeForcesY, const real* _bladeForcesZ)
{
    std::copy_n(_bladeForcesX, this->numberOfGridNodes, this->bladeForcesXH);
    std::copy_n(_bladeForcesY, this->numberOfGridNodes, this->bladeForcesYH);
    std::copy_n(_bladeForcesZ, this->numberOfGridNodes, this->bladeForcesZH);
}

void ActuatorFarm::setTurbineBladeCoords(uint turbine, const real* _bladeCoordsX, const real* _bladeCoordsY,
                                         const real* _bladeCoordsZ)
{
    std::copy_n(_bladeCoordsX, this->numberOfNodesPerTurbine, &this->bladeCoordsXH[turbine * this->numberOfNodesPerTurbine]);
    std::copy_n(_bladeCoordsY, this->numberOfNodesPerTurbine, &this->bladeCoordsYH[turbine * this->numberOfNodesPerTurbine]);
    std::copy_n(_bladeCoordsZ, this->numberOfNodesPerTurbine, &this->bladeCoordsZH[turbine * this->numberOfNodesPerTurbine]);
}

void ActuatorFarm::setTurbineBladeVelocities(uint turbine, const real* _bladeVelocitiesX, const real* _bladeVelocitiesY,
                                             const real* _bladeVelocitiesZ)
{
    std::copy_n(_bladeVelocitiesX, this->numberOfNodesPerTurbine, &this->bladeVelocitiesXH[turbine * this->numberOfNodesPerTurbine]);
    std::copy_n(_bladeVelocitiesY, this->numberOfNodesPerTurbine, &this->bladeVelocitiesYH[turbine * this->numberOfNodesPerTurbine]);
    std::copy_n(_bladeVelocitiesZ, this->numberOfNodesPerTurbine, &this->bladeVelocitiesZH[turbine * this->numberOfNodesPerTurbine]);
}

void ActuatorFarm::setTurbineBladeForces(uint turbine, const real* _bladeForcesX, const real* _bladeForcesY,
                                         const real* _bladeForcesZ)
{
    std::copy_n(_bladeForcesX, this->numberOfNodesPerTurbine, &this->bladeForcesXH[turbine * this->numberOfNodesPerTurbine]);
    std::copy_n(_bladeForcesY, this->numberOfNodesPerTurbine, &this->bladeForcesYH[turbine * this->numberOfNodesPerTurbine]);
    std::copy_n(_bladeForcesZ, this->numberOfNodesPerTurbine, &this->bladeForcesZH[turbine * this->numberOfNodesPerTurbine]);
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

std::string ActuatorFarm::getFilename(Parameter* para, uint t) const
{
    return para->getOutputPath() + this->outputName + "_ID_" + std::to_string(para->getMyProcessID()) + "_t_" + std::to_string(t);
}

void ActuatorFarm::write(const std::string& filename) const
{
    std::vector<std::string> dataNames = {
        "bladeVelocitiesX",
        "bladeVelocitiesY",
        "bladeVelocitiesZ",
        "bladeForcesX",
        "bladeForcesY",
        "bladeForcesZ",
    };
    std::vector<UbTupleFloat3> nodes(numberOfGridNodes);
    std::vector<std::vector<double>> nodeData(6);
    for (auto& data : nodeData)
        data.resize(numberOfGridNodes);
    for (uint i = 0; i < numberOfGridNodes; i++) {
        nodes[i] = UbTupleFloat3(this->bladeCoordsXH[i], this->bladeCoordsYH[i], this->bladeCoordsZH[i]);
        nodeData[0][i] = this->bladeVelocitiesXH[i];
        nodeData[1][i] = this->bladeVelocitiesYH[i];
        nodeData[2][i] = this->bladeVelocitiesZH[i];
        nodeData[3][i] = this->bladeForcesXH[i];
        nodeData[4][i] = this->bladeForcesYH[i];
        nodeData[5][i] = this->bladeForcesZH[i];
    }
    WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filename, nodes, dataNames, nodeData);
}

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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_Samplers Samplers
//! \ingroup gpu_core core
//! \{
//! \author Henrik Asmuth, Henry Korb
//======================================================================================
#include "PrecursorWriter.h"

#include <algorithm>
#include <memory>

#include <basics/constants/NumericConstants.h>
#include <basics/writer/WbWriterVtkXmlImageBinary.h>
#include <logger/Logger.h>

#include "Cuda/CudaMemoryManager.h"
#include "Cuda/CudaStreamManager.h"
#include "DataStructureInitializer/GridProvider.h"
#include "Output/FilePartCalculator.h"
#include "Parameter/Parameter.h"
#include "StringUtilities/StringUtil.h"
#include "Utilities/KernelUtilities.h"
#include "cuda_helper/CudaGrid.h"
#include "cuda_helper/CudaIndexCalculation.h"
#include "helper_cuda.h"

using namespace vf::lbm::dir;
using namespace vf::gpu;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO check everything for multiple level
int index1d(int y, int z, int ny, int nz)
{
    return y + ny * z;
}

constexpr uint linearIdx(const uint component, const uint node, const uint timeStep, const uint numberOfComponents,
                         const uint numberOfNodes)
{
    return node + numberOfNodes * (component + numberOfComponents * timeStep);
}

constexpr uint linearIdx(const uint component, const uint node, const uint numberOfNodes)
{
    return node + component * numberOfNodes;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void fillArrayVelocities(uint numberOfPrecursorNodes, const uint* indices, real* precursorData, const real* vx,
                                    const real* vy, const real* vz, real velocityRatio)

{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= numberOfPrecursorNodes)
        return;

    precursorData[linearIdx(0U, nodeIndex, numberOfPrecursorNodes)] = vx[indices[nodeIndex]] * velocityRatio;
    precursorData[linearIdx(1U, nodeIndex, numberOfPrecursorNodes)] = vy[indices[nodeIndex]] * velocityRatio;
    precursorData[linearIdx(2U, nodeIndex, numberOfPrecursorNodes)] = vz[indices[nodeIndex]] * velocityRatio;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void fillArrayDistributions(uint numberOfPrecursorNodes, const uint* indices, real* precursorData,
                                       real* distributions, const uint* neighborX, const uint* neighborY,
                                       const uint* neighborZ, bool isEvenTimeStep, unsigned long numberOfLBnodes)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= numberOfPrecursorNodes)
        return;

    Distributions27 dist;
    getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimeStep);

    ////////////////////////////////////////////////////////////////////////////////
    // ! - Set neighbor indices (necessary for indirect addressing)
    uint k_000 = indices[nodeIndex];
    // uint k_M00 = neighborX[k_000];
    uint k_0M0 = neighborY[k_000];
    uint k_00M = neighborZ[k_000];
    // uint k_MM0 = neighborY[k_M00];
    // uint k_M0M = neighborZ[k_M00];
    uint k_0MM = neighborZ[k_0M0];
    // uint k_MMM = neighborZ[k_MM0];

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Get local distributions in PX directions
    //!
    precursorData[linearIdx(PrecursorWriter::dP00, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dP00])[k_000];
    precursorData[linearIdx(PrecursorWriter::dPP0, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPP0])[k_000];
    precursorData[linearIdx(PrecursorWriter::dPM0, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPM0])[k_0M0];
    precursorData[linearIdx(PrecursorWriter::dP0P, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dP0P])[k_000];
    precursorData[linearIdx(PrecursorWriter::dP0M, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dP0M])[k_00M];
    precursorData[linearIdx(PrecursorWriter::dPPP, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPPP])[k_000];
    precursorData[linearIdx(PrecursorWriter::dPMP, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPMP])[k_0M0];
    precursorData[linearIdx(PrecursorWriter::dPPM, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPPM])[k_00M];
    precursorData[linearIdx(PrecursorWriter::dPMM, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPMM])[k_0MM];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PrecursorWriter::init()
{
    VF_LOG_INFO("PrecursorWriter: Start initializing...");
    VF_LOG_INFO("Writing yz-planes at x={}m every {}. timestep, starting at t={}", this->xPos, this->tSave,
                this->tStartWritingOutput);
    VF_LOG_INFO("Writing to {} with prefix {}", outputPath, probeName);

    precursorStructs.resize(para->getMaxLevel() + 1);
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        const real deltaX = abs(para->getParH(level)->coordinateX[para->getParH(level)->neighborX[1]] -
                                para->getParH(level)->coordinateX[1]);

        real lowestY, lowestZ, highestY, highestZ;

        lowestY = para->getParH(level)->coordinateY[para->getParH(level)->numberOfNodes - 1];
        highestY = para->getParH(level)->coordinateY[1];

        lowestZ = para->getParH(level)->coordinateZ[para->getParH(level)->numberOfNodes - 1];
        highestZ = para->getParH(level)->coordinateZ[1];

        std::vector<uint> indicesOnGrid;
        std::vector<int> indicesOnPlane;
        std::vector<real> coordY, coordZ;

        for (size_t pos = 1; pos < para->getParH(level)->numberOfNodes; pos++) {
            const real pointCoordX = para->getParH(level)->coordinateX[pos];
            const real pointCoordY = para->getParH(level)->coordinateY[pos];
            const real pointCoordZ = para->getParH(level)->coordinateZ[pos];
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID && pointCoordX < (deltaX + xPos) &&
                pointCoordX >= xPos && pointCoordY <= yMax && pointCoordY >= yMin && pointCoordZ <= zMax &&
                pointCoordZ >= zMin) {
                highestY = std::max(highestY, pointCoordY);
                highestZ = std::max(highestZ, pointCoordZ);

                lowestY = std::min(lowestY, pointCoordY);
                lowestZ = std::min(lowestZ, pointCoordZ);
                indicesOnGrid.push_back((uint)pos);
                coordY.push_back(pointCoordY);
                coordZ.push_back(pointCoordZ);
            }
        }
        if (indicesOnGrid.size() == 0)
            throw std::runtime_error("PrecursorWriter did not find any points on the grid");

        const int ny = int((highestY - lowestY) / deltaX) + 1;
        const int nz = int((highestZ - lowestZ) / deltaX) + 1;

        precursorStructs[level].indicesOnPlane.reserve(indicesOnGrid.size());

        for (uint i = 0; i < indicesOnGrid.size(); i++) {
            const int idxY = int((coordY[i] - lowestY) / deltaX);
            const int idxZ = int((coordZ[i] - lowestZ) / deltaX);
            precursorStructs[level].indicesOnPlane.push_back(index1d(idxY, idxZ, ny, nz));
        }

        precursorStructs[level].numberOfPointsInBC = (uint)indicesOnGrid.size();
        precursorStructs[level].spacing = makeUbTuple(deltaX, deltaX, tSave * para->getTimeRatio() * pow(2, -level));
        precursorStructs[level].origin = makeUbTuple(lowestY, lowestZ);
        precursorStructs[level].extent = makeUbTuple(0, ny - 1, 0, nz - 1);
        precursorStructs[level].numberOfPointsInData = ny * nz;
        precursorStructs[level].numberOfTimeStepsPerFile =
            std::min(FilePartCalculator::limitOfNodesForVTK / (ny * nz), maxTimestepsPerFile);
        precursorStructs[level].numberOfFilesWritten = 0;
        precursorStructs[level].numberOfTimeStepsBuffered = 0;
        precursorStructs[level].streamIndex =
            para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::PrecursorWriter);

        switch (outputVariable) {
            case OutputVariable::Velocities:
                precursorStructs[level].numberOfQuantities = 3;
                break;
            case OutputVariable::Distributions:
                precursorStructs[level].numberOfQuantities = 9;
                break;
            default:
                break;
        }

        cudaMemoryManager->cudaAllocPrecursorWriter(this, level);

        std::copy(indicesOnGrid.begin(), indicesOnGrid.end(), precursorStructs[level].indicesH);

        cudaMemoryManager->cudaCopyPrecursorWriterIndicesHtoD(this, level);

        VF_LOG_INFO("Found {} points in precursor plane on level {}", precursorStructs[level].numberOfPointsInBC, level);
    }
    VF_LOG_INFO("PrecursorWriter: Done initializing.");
}

void PrecursorWriter::sample(int level, uint t)
{
    const uint t_level = para->getTimeStep(level, t, true);
    const uint tStartOut_level = tStartWritingOutput * std::exp2(level);
    auto* levelStruct = getPrecursorStruct(level);

    if (t_level > tStartOut_level && ((t_level - tStartOut_level) % tSave) == 0) {
        vf::cuda::CudaGrid grid(para->getParH(level)->numberofthreads, levelStruct->numberOfPointsInBC);

        if (this->outputVariable == OutputVariable::Velocities) {
            fillArrayVelocities<<<grid.grid, grid.threads>>>(levelStruct->numberOfPointsInBC, levelStruct->indicesD,
                                                             levelStruct->bufferD, para->getParD(level)->velocityX,
                                                             para->getParD(level)->velocityY,
                                                             para->getParD(level)->velocityZ, para->getVelocityRatio());
            getLastCudaError("In PrecursorWriter::interact fillArrayVelocities execution failed");
        } else if (this->outputVariable == OutputVariable::Distributions) {
            fillArrayDistributions<<<grid.grid, grid.threads>>>(
                levelStruct->numberOfPointsInBC, levelStruct->indicesD, levelStruct->bufferD,
                para->getParD(level)->distributions.f[0], para->getParD(level)->neighborX, para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ, para->getEvenOrOdd(level), para->getParD(level)->numberOfNodes);
            getLastCudaError("In PrecursorWriter::interact fillArrayDistributions execution failed");
        }
        cudaMemoryManager->cudaCopyPrecursorWriterOutputVariablesDtoH(this, level);

        // switch device buffer and data pointer so precursor data is gathered in buffer and copied from bufferD to bufferH
        real* tmp = levelStruct->bufferD;
        levelStruct->bufferD = levelStruct->dataD;
        levelStruct->dataD = tmp;

        levelStruct->numberOfTimeStepsBuffered++;

        if (levelStruct->numberOfTimeStepsBuffered >= levelStruct->numberOfTimeStepsPerFile || t == para->getTimestepEnd()) {
            // switch host buffer and data pointer so precursor data is copied in buffer and written from data

            tmp = levelStruct->bufferH;
            levelStruct->bufferH = levelStruct->dataH;
            levelStruct->dataH = tmp;

            writeFuture.wait();
            writeFuture = std::async(
                std::launch::async, [this](uint level, uint timeSteps) { this->write(level, timeSteps); }, level,
                levelStruct->numberOfTimeStepsBuffered);
            levelStruct->numberOfTimeStepsBuffered = 0;
        }
    }
}

PrecursorWriter::~PrecursorWriter()
{
    writeFuture.wait();
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        if (getPrecursorStruct(level)->numberOfTimeStepsBuffered > 0)
            write(level, getPrecursorStruct(level)->numberOfTimeStepsBuffered);

        cudaMemoryManager->cudaFreePrecursorWriter(this, level);
    }
}

void PrecursorWriter::write(int level, uint numberOfTimeStepsBuffered)
{
    const std::string fileName = makeFileName(level, precursorStructs[level].numberOfFilesWritten);

    const std::size_t numberOfPointsInData = precursorStructs[level].numberOfPointsInData;

    const int startTime = precursorStructs[level].numberOfFilesWritten * precursorStructs[level].numberOfTimeStepsPerFile;

    UbTupleInt6 extent = makeUbTuple(val<1>(precursorStructs[level].extent), val<2>(precursorStructs[level].extent),
                                     val<3>(precursorStructs[level].extent), val<4>(precursorStructs[level].extent),
                                     startTime, startTime + (int)numberOfTimeStepsBuffered - 1);

    UbTupleFloat3 origin = makeUbTuple(static_cast<float>(val<1>(precursorStructs[level].origin)),
                                       static_cast<float>(val<2>(precursorStructs[level].origin)), 0.F);

    std::vector<std::vector<double>> nodeData;

    for (uint quant = 0; quant < precursorStructs[level].numberOfQuantities; quant++) {
        std::vector<double> doubleArr(numberOfPointsInData * numberOfTimeStepsBuffered, NAN);
        for (uint timeStep = 0; timeStep < numberOfTimeStepsBuffered; timeStep++) {
            for (uint pos = 0; pos < precursorStructs[level].numberOfPointsInBC; pos++) {
                const std::size_t indexOnPlane = precursorStructs[level].indicesOnPlane[pos] + timeStep * numberOfPointsInData;
                doubleArr[indexOnPlane] = double(
                    precursorStructs[level].dataH[linearIdx(quant, pos, timeStep, precursorStructs[level].numberOfQuantities,
                                                            precursorStructs[level].numberOfPointsInBC)]);
            }
        }
        nodeData.push_back(doubleArr);
    }

    std::vector<std::vector<double>> cellData;
    getWriter()->writeData(fileName, nodeDataNames, cellDataNames, nodeData, cellData, extent, origin,
                           precursorStructs[level].spacing, extent, this->writePrecision);
    VF_LOG_INFO("PrecursorWriter: Wrote file {}", fileName);
    precursorStructs[level].numberOfFilesWritten++;
}

std::string PrecursorWriter::makeFileName(int level, uint numberOfFilesWritten)
{
    return outputPath + probeName + "_lev_" + StringUtil::toString<int>(level) + "_ID_" +
           StringUtil::toString<int>(para->getMyProcessID()) + "_File_" + StringUtil::toString<int>(numberOfFilesWritten) +
           getWriter()->getFileExtension();
}

void PrecursorWriter::getTaggedFluidNodes(GridProvider* gridProvider)
{
    for (uint level = 0; level < (uint)para->getMaxLevel(); level++) {
        if (outputVariable == OutputVariable::Velocities) {
            std::vector<uint> indices(precursorStructs[level].indicesH,
                                      precursorStructs[level].indicesH + precursorStructs[level].numberOfPointsInBC);
            gridProvider->tagFluidNodeIndices(indices, CollisionTemplate::WriteMacroVars, level);
        }
    }
}
//! \}

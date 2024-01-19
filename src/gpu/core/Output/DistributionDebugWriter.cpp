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
//! \addtogroup gpu_Output Output
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "DistributionDebugWriter.h"

#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <logger/Logger.h>

#include <lbm/constants/D3Q27.h>

#include "Cuda/CudaMemoryManager.h"
#include "FilePartCalculator.h"
#include "Parameter/Parameter.h"
#include "WriterUtilities.h"

using namespace vf::lbm::dir;

void DistributionDebugWriter::writeDistributions(const Parameter& para, uint timestep)
{
    for (int level = para.getCoarse(); level <= para.getFine(); level++) {
        DistributionDebugWriter::writeDistributionsForLevel(para, level, timestep);
    }
}

void createFileNames(std::vector<std::string>& fileNames, uint numberOfParts, uint level, uint timestep, const Parameter& para)
{
    for (uint i = 1; i <= numberOfParts; i++) {
        fileNames.push_back(para.getFName() + "_bin_distributions" +
                            WriterUtilities::makePartFileNameEnding(level, para.getMyProcessID(), i, timestep));
    }
}

void createNodeDataNames(std::vector<std::string>& nodeDataNames)
{
    nodeDataNames.resize(NUMBER_Of_DIRECTIONS);

    for (uint dir = STARTDIR; dir <= ENDDIR; dir++) {
        const size_t minLenghtOfNumberString = 2; // the number is padded with zeros to this length
        const auto numberString = std::to_string(dir);
        nodeDataNames[dir] =
            "f_" + std::string(minLenghtOfNumberString - std::min(minLenghtOfNumberString, numberString.length()), '0') +
            numberString;
    }
}

void DistributionDebugWriter::writeDistributionsForLevel(const Parameter& para, uint level, uint timestep)
{
    const LBMSimulationParameter& parH = para.getParHostAsReference(level);
    const uint numberOfParts = FilePartCalculator::calculateNumberOfParts(parH.numberOfNodes);

    std::vector<std::string> fileNames;
    createFileNames(fileNames, numberOfParts, level, timestep, para);

    std::vector<std::string> nodeDataNames;
    createNodeDataNames(nodeDataNames);

    uint sizeOfNodes;
    uint startPosition;
    uint endPosition;
    std::array<uint, 8> indicesOfOct;
    std::array<uint, 8> relativePosInPart;
    uint relativePositionInPart;

    Distributions27 distributions = parH.distributions;

    if (distributions.f[0] == nullptr)
        throw std::runtime_error("Distributions (distributions.f[0]) at level " + std::to_string(level) +
                                 " are not allocated on the host. Can't write distributions.");

    for (unsigned int part = 0; part < numberOfParts; part++) {
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(parH.numberOfNodes, part);
        startPosition = FilePartCalculator::calculateStartingPostionOfPart(part);
        endPosition = startPosition + sizeOfNodes;

        std::vector<UbTupleFloat3> nodes(sizeOfNodes);
        std::vector<UbTupleUInt8> cells;
        std::vector<std::vector<double>> nodeData(nodeDataNames.size());
        for (uint i = 0; i < (uint)nodeDataNames.size(); i++)
            nodeData[i].resize(sizeOfNodes);

        for (unsigned int pos = startPosition; pos < endPosition; pos++) {

            if (parH.typeOfGridNode[pos] != GEO_FLUID)
                continue;

            relativePositionInPart = pos - startPosition;

            nodes[relativePositionInPart] =
                makeUbTuple((float)parH.coordinateX[pos], (float)parH.coordinateY[pos], (float)parH.coordinateZ[pos]);

            for (uint dir = STARTDIR; dir <= ENDDIR; dir++) {
                nodeData[dir][relativePositionInPart] = distributions.f[0][dir * parH.numberOfNodes + pos];
            }

            WriterUtilities::getIndicesOfAllNodesInOct(indicesOfOct, pos, parH);
            if (WriterUtilities::isPeriodicCell(parH, indicesOfOct[0], indicesOfOct[6])) {
                continue;
            }

            if (WriterUtilities::areAllNodesInOctValidForWriting(indicesOfOct, parH, endPosition)) {
                WriterUtilities::calculateRelativeNodeIndexInPart(relativePosInPart, indicesOfOct, startPosition);
                cells.push_back(makeUbTupleFromArray(relativePosInPart));
            }
        }

        std::string fileName = WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fileNames[part], nodes, cells,
                                                                                          nodeDataNames, nodeData);
        VF_LOG_DEBUG("DistributionDebugWriter wrote to {} ", fileName);
    }
}

void DistributionDebugWriter::allocateDistributionsOnHost(const CudaMemoryManager& cudaMemoryManager)
{
    cudaMemoryManager.cudaAllocFsForAllLevelsOnHost();
}

void DistributionDebugWriter::allocateDistributionsOnHost(const CudaMemoryManager& cudaMemoryManager, uint level)
{
    cudaMemoryManager.cudaAllocFsForCheckPointAndRestart(level);
}

void DistributionDebugWriter::copyDistributionsToHost(const Parameter& para, const CudaMemoryManager& cudaMemoryManager)
{
    for (int level = 0; level <= para.getMaxLevel(); level++)
        DistributionDebugWriter::copyDistributionsToHost(para, cudaMemoryManager, level);
}

void DistributionDebugWriter::copyDistributionsToHost(const Parameter& para, const CudaMemoryManager& cudaMemoryManager,
                                                      uint level)
{

    if (para.getParHostAsReference(level).distributions.f[0] == nullptr)
        throw std::runtime_error("Distributions (distributions.f[0]) at level " + std::to_string(level) +
                                 " are not allocated on the host. Can't copy distributions to host");
    cudaMemoryManager.cudaCopyFsForCheckPoint(level);
}

//! \}

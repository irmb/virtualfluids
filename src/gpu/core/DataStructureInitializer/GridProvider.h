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
//! \addtogroup gpu_DataStructureInitializer DataStructureInitializer
//! \ingroup gpu_core core
//! \{
//=======================================================================================
#ifndef GridReader_H
#define GridReader_H

#include <string>
#include <vector>
#include <memory>

#include "Calculation/Calculation.h"
#include <basics/PointerDefinitions.h>

#include "gpu/GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"
namespace vf::parallel
{
class Communicator;
}

class Parameter;
class GridBuilder;
class CudaMemoryManager;

class GridProvider
{
public:
    static std::shared_ptr<GridProvider> makeGridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager, vf::parallel::Communicator& communicator);
    static std::shared_ptr<GridProvider> makeGridReader(FILEFORMAT format, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager);

    virtual void allocArrays_CoordNeighborGeo() = 0;
    virtual void allocArrays_BoundaryValues() = 0;
    virtual void allocArrays_BoundaryQs() = 0;
    virtual void allocArrays_OffsetScale() = 0;
    virtual void allocArrays_taggedFluidNodes() = 0;

    virtual void tagFluidNodeIndices(const std::vector<uint>& taggedFluidNodeIndices, CollisionTemplate tag, uint level) = 0;
    virtual void sortFluidNodeTags() = 0;

    virtual void setDimensions() = 0;
    virtual void setBoundingBox() = 0;
    virtual void initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int> > > periodV, std::vector<std::vector<unsigned int> > periodIndex, std::string way) = 0;

    virtual void allocAndCopyForcing();
    virtual void allocAndCopyQuadricLimiters();
    virtual void freeMemoryOnHost();
    virtual void cudaCopyDataToHost(int level);

    virtual ~GridProvider() = default;
    virtual void initalGridInformations() = 0;

protected:
    void setNumberOfNodes(uint numberOfNodes, int level) const;
    void setNumberOfTaggedFluidNodes(uint numberOfNodes, CollisionTemplate tag, int level) const;
    virtual void setInitialNodeValues(uint numberOfNodes, int level) const;

    void setPressSizePerLevel(int level, int sizePerLevel) const;
    void setVelocitySizePerLevel(int level, int sizePerLevel) const;
    void setOutflowSizePerLevel(int level, int sizePerLevel) const;

    std::shared_ptr<Parameter> para;
    std::shared_ptr<CudaMemoryManager> cudaMemoryManager;
};

#endif

//! \}

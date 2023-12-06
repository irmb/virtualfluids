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
//=======================================================================================
#ifndef GridReaderFiles_H
#define GridReaderFiles_H

#include "../GridProvider.h"

#include <memory>
#include <string>
#include <vector>

#include "Calculation/Calculation.h"

class Parameter;
class BoundaryValues;
class BoundaryQs;
class CoordNeighborGeoV;

class GridReader : public GridProvider
{
private:
    bool binaer;
    std::vector<std::string> channelDirections;
    std::vector<std::string> channelBoundaryConditions;
    std::shared_ptr<CoordNeighborGeoV> neighX, neighY, neighZ, neighWSB;
    std::vector<std::shared_ptr<BoundaryValues>> BC_Values;

    std::vector<std::vector<real>> velocityX_BCvalues, velocityY_BCvalues, velocityZ_BCvalues;
    std::vector<std::vector<std::vector<real>>> velocityQs;
    std::vector<std::vector<int>> velocityIndex;

    std::vector<std::vector<real>> pressureBCvalues;
    std::vector<std::vector<real>> outflowBCvalues;

public:
    GridReader(FILEFORMAT format, std::shared_ptr<Parameter> para,
               std::shared_ptr<CudaMemoryManager> cudaMemoryManager);
     ~GridReader() = default;
    void allocArrays_CoordNeighborGeo() override;
    void allocArrays_BoundaryValues() override;
    void allocArrays_OffsetScale() override;
    void allocArrays_taggedFluidNodes() override;

    void tagFluidNodeIndices(const std::vector<uint> &taggedFluidNodeIndices, CollisionTemplate tag, uint level) override;

    void sortFluidNodeTags() override;

    void initalValuesDomainDecompostion(int level);

    void setChannelBoundaryCondition();

    void allocArrays_BoundaryQs() override;
    bool getBinaer();
    void setDimensions() override;
    void setBoundingBox() override;
    void initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int>>> periodV,
                           std::vector<std::vector<unsigned int>> periodIndex, std::string way) override;

private:
    void makeReader(std::shared_ptr<Parameter> para);
    void makeReader(std::vector<std::shared_ptr<BoundaryQs>> &BC_Qs, std::shared_ptr<Parameter> para);

    void setPressureValues(int channelSide) const;
    void setPressRhoBC(int sizePerLevel, int level, int channelSide) const;

    void fillVelocityVectors(int channelSide);
    void setVelocityValues();
    void setVelocity(int level, int sizePerLevel) const;

    void setOutflowValues(int channelSide) const;
    void setOutflow(int level, int sizePerLevel, int channelSide) const;

    // void fillVelocityQVectors(int channelSide);
    void setPressQs(std::shared_ptr<BoundaryQs> boundaryQ) const;
    void setVelocityQs(std::shared_ptr<BoundaryQs> boundaryQ);
    void setOutflowQs(std::shared_ptr<BoundaryQs> boundaryQ) const;
    void setNoSlipQs(std::shared_ptr<BoundaryQs> boundaryQ) const;
    void setGeoQs(std::shared_ptr<BoundaryQs> boundaryQ) const;
    void modifyQElement(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;

    void initalVectorForQStruct(std::vector<std::vector<std::vector<real>>> &Qs, std::vector<std::vector<int>> &index,
                                std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
    void copyVectorsToQStruct(std::vector<std::vector<real>> &Qs, std::vector<int> &index,
                              QforBoundaryConditions &Q) const;
    void initalQStruct(QforBoundaryConditions &Q, std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
    void printQSize(std::string bc, std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
    void setSizeNoSlip(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
    void setSizeGeoQs(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
    void setQ27Size(QforBoundaryConditions &Q, real *QQ, unsigned int sizeQ) const;
    bool hasQs(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;

public:
    void initalGridInformations() override;
};

#endif

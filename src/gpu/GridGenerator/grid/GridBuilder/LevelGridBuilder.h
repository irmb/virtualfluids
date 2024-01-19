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
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz, Martin Sch�nherr
//=======================================================================================
#ifndef LEVEL_GRID_BUILDER_H
#define LEVEL_GRID_BUILDER_H

#include <vector>
#include <string>
#include <memory>
#include <array>

#include <basics/constants/NumericConstants.h>

#include "gpu/GridGenerator/global.h"

#include "gpu/GridGenerator/grid/GridBuilder/GridBuilder.h"
#include "gpu/GridGenerator/grid/Grid.h"
#include "gpu/GridGenerator/grid/GridInterface.h"
#include "gpu/GridGenerator/grid/NodeValues.h"

using namespace vf::basics::constant;

struct Vertex;
class Grid;
class Transformator;
class ArrowTransformator;
class PolyDataWriterWrapper;
class BoundingBox;

class Side;
class VelocityBoundaryCondition;
class SlipBoundaryCondition;
class StressBoundaryCondition;
class PressureBoundaryCondition;
class GeometryBoundaryCondition;
class PrecursorBoundaryCondition;
enum class SideType;

class TransientBCInputFileReader;
class FileCollection;

class LevelGridBuilder : public GridBuilder
{
protected:
    LevelGridBuilder();

public:
    static std::shared_ptr<LevelGridBuilder> makeShared();

    SPtr<Grid> getGrid(uint level) override;

     ~LevelGridBuilder() override;

    virtual void setSlipBoundaryCondition(SideType sideType, real nomalX, real normalY, real normalZ);
    virtual void setStressBoundaryCondition(SideType sideType, real nomalX, real normalY, real normalZ,
                                                                 uint samplingOffset, real z0, real dx);
    virtual void setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz);
    virtual void setPressureBoundaryCondition(SideType sideType, real rho);
    virtual void setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z);
    void setPeriodicShiftOnXBoundaryInYDirection(real shift);
    void setPeriodicShiftOnXBoundaryInZDirection(real shift);
    void setPeriodicShiftOnYBoundaryInXDirection(real shift);
    void setPeriodicShiftOnYBoundaryInZDirection(real shift);
    void setPeriodicShiftOnZBoundaryInXDirection(real shift);
    void setPeriodicShiftOnZBoundaryInYDirection(real shift);
    virtual void setNoSlipBoundaryCondition(SideType sideType);
    virtual void setPrecursorBoundaryCondition(SideType sideType, SPtr<FileCollection> fileCollection,
                                                                    int timeStepsBetweenReads, real velocityX = c0o1,
                                                                    real velocityY = c0o1, real velocityZ = c0o1,
                                                                    std::vector<uint> fileLevelToGridLevelMap = {});

    void setEnableFixRefinementIntoTheWall(bool enableFixRefinementIntoTheWall);

    virtual void setCommunicationProcess(int direction, uint process);

    virtual uint getCommunicationProcess(int direction) override;

    std::shared_ptr<Grid> getGrid(int level, int box);

    virtual unsigned int getNumberOfNodes(unsigned int level) const override;

    virtual uint getNumberOfFluidNodes(unsigned int level) const override;
    virtual void getFluidNodeIndices(uint* fluidNodeIndices, const int level) const override;
    virtual uint getNumberOfFluidNodesBorder(unsigned int level) const override;
    virtual void getFluidNodeIndicesBorder(uint *fluidNodeIndices, const int level) const override;

    virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords,
                                         uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative,
                                         uint *geo, const int level) const override;
    virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const override;


    uint getSlipSize(int level) const override;
    virtual void getSlipValues(real* normalX, real* normalY, real* normalZ, int* indices, int level) const override;
    virtual void getSlipQs(real* qs[27], int level) const override;

    uint getStressSize(int level) const override;
    virtual void getStressValues(  real* normalX, real* normalY, real* normalZ,
                                                        real* vx,      real* vy,      real* vz,
                                                        real* vx1,     real* vy1,     real* vz1,
                                                        int* indices, int* samplingIndices, int* samplingOffsets, real* z0, int level) const override;
    virtual void getStressQs(real* qs[27], int level) const override;

    uint getVelocitySize(int level) const override;
    virtual void getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const override;
    virtual void getVelocityQs(real* qs[27], int level) const override;

    uint getPressureSize(int level) const override;
    void getPressureValues(real* rho, int* indices, int* neighborIndices, int level) const override;
    virtual void getPressureQs(real* qs[27], int level) const override;

    uint getPrecursorSize(int level) const override;
    void getPrecursorValues(   uint* neighbor0PP, uint* neighbor0PM, uint* neighbor0MP, uint* neighbor0MM,
                                                    real* weights0PP, real* weights0PM, real* weights0MP, real* weights0MM,
                                                    int* indices, std::vector<SPtr<TransientBCInputFileReader>>& reader,
                                                    int& numberOfPrecursorNodes, size_t& numberOfQuantities, uint& timeStepsBetweenReads,
                                                    real& velocityX, real& velocityY, real& velocityZ, int level) const override;
    virtual void getPrecursorQs(real* qs[27], int level) const override;

    virtual void getGeometryQs(real *qs[27], int level) const override;
    virtual uint getGeometrySize(int level) const override;
    virtual void getGeometryIndices(int *indices, int level) const override;
    virtual bool hasGeometryValues() const override;
    virtual void getGeometryValues(real *vx, real *vy, real *vz, int level) const override;

    void writeArrows(std::string fileName) const override;

    SPtr<gg::BoundaryCondition> getBoundaryCondition( SideType side, uint level ) const override;
    SPtr<GeometryBoundaryCondition> getGeometryBoundaryCondition(uint level) const override;

protected:


    struct BoundaryConditions
    {
        BoundaryConditions() = default;

        std::vector<SPtr<SlipBoundaryCondition>> slipBoundaryConditions;

        std::vector<SPtr<StressBoundaryCondition>> stressBoundaryConditions;

        std::vector<SPtr<VelocityBoundaryCondition>> velocityBoundaryConditions;

        std::vector<SPtr<PressureBoundaryCondition>> pressureBoundaryConditions;

        std::vector<SPtr<VelocityBoundaryCondition>> noSlipBoundaryConditions;

        std::vector<SPtr<PrecursorBoundaryCondition>> precursorBoundaryConditions;

        SPtr<GeometryBoundaryCondition> geometryBoundaryCondition;
    };
    bool geometryHasValues = false;

    std::vector<std::shared_ptr<Grid> > grids;
    std::vector<SPtr<BoundaryConditions> > boundaryConditions;

    std::array<uint, 6> communicationProcesses;

    void checkLevel(int level);

protected:
    void setVelocityGeometryBoundaryCondition(real vx, real vy, real vz);
    void setNoSlipGeometryBoundaryCondition();
    void setSlipGeometryBoundaryCondition(real normalX, real normalY, real normalZ);

    void createBCVectors();
    void addShortQsToVector(int index);
    void addQsToVector(int index);
    void fillRBForNode(int index, int direction, int directionSign, int rb);

    Vertex getVertex(const int matrixIndex) const;

public:
    void getGridInformations(std::vector<int>& gridX, std::vector<int>& gridY,
                                       std::vector<int>& gridZ, std::vector<int>& distX, std::vector<int>& distY,
                                       std::vector<int>& distZ) override;
    virtual uint getNumberOfGridLevels() const override;

    uint getNumberOfNodesCF(int level) override;
    uint getNumberOfNodesFC(int level) override;

    void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const override;

    void getOffsetFC(real* xOffCf, real* yOffCf, real* zOffCf, int level) override;
    void getOffsetCF(real* xOffFc, real* yOffFc, real* zOffFc, int level) override;

    virtual uint getNumberOfSendIndices(int direction, uint level) override;
    virtual uint getNumberOfReceiveIndices(int direction, uint level) override;
    virtual void getSendIndices(int *sendIndices, int direction, int level) override;
    virtual void getReceiveIndices(int *sendIndices, int direction, int level) override;


    // needed for CUDA Streams MultiGPU (Communication Hiding)
    void findFluidNodes(bool splitDomain) override;

    void addFluidNodeIndicesMacroVars(const std::vector<uint>& fluidNodeIndicesMacroVars, uint level) override;
    void addFluidNodeIndicesApplyBodyForce(const std::vector<uint>& fluidNodeIndicesApplyBodyForce, uint level) override;
    void addFluidNodeIndicesAllFeatures(const std::vector<uint>& fluidNodeIndicesAllFeatures, uint level) override;

    void sortFluidNodeIndicesMacroVars(uint level) override;
    void sortFluidNodeIndicesApplyBodyForce(uint level) override;
    void sortFluidNodeIndicesAllFeatures(uint level) override;

    uint getNumberOfFluidNodesMacroVars(unsigned int level) const override;
    void getFluidNodeIndicesMacroVars(uint *fluidNodeIndicesMacroVars, const int level) const override;
    uint getNumberOfFluidNodesApplyBodyForce(unsigned int level) const override;
    void getFluidNodeIndicesApplyBodyForce(uint *fluidNodeIndicesApplyBodyForce, const int level) const override;
    uint getNumberOfFluidNodesAllFeatures(unsigned int level) const override;
    void getFluidNodeIndicesAllFeatures(uint *fluidNodeIndicesAllFeatures, const int level) const override;
};

#endif

//! \}

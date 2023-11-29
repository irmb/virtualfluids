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
//! \file LevelGridBuilder.h
//! \ingroup grid
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
    GRIDGENERATOR_EXPORT LevelGridBuilder();

public:
    GRIDGENERATOR_EXPORT static std::shared_ptr<LevelGridBuilder> makeShared();

    GRIDGENERATOR_EXPORT SPtr<Grid> getGrid(uint level) override;

    GRIDGENERATOR_EXPORT  ~LevelGridBuilder() override;

    GRIDGENERATOR_EXPORT virtual void setSlipBoundaryCondition(SideType sideType, real nomalX, real normalY, real normalZ);
    GRIDGENERATOR_EXPORT virtual void setStressBoundaryCondition(SideType sideType, real nomalX, real normalY, real normalZ,
                                                                 uint samplingOffset, real z0, real dx);
    GRIDGENERATOR_EXPORT virtual void setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz);
    GRIDGENERATOR_EXPORT virtual void setPressureBoundaryCondition(SideType sideType, real rho);
    GRIDGENERATOR_EXPORT virtual void setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z);
    GRIDGENERATOR_EXPORT void setPeriodicShiftOnXBoundaryInYDirection(real shift);
    GRIDGENERATOR_EXPORT void setPeriodicShiftOnXBoundaryInZDirection(real shift);
    GRIDGENERATOR_EXPORT void setPeriodicShiftOnYBoundaryInXDirection(real shift);
    GRIDGENERATOR_EXPORT void setPeriodicShiftOnYBoundaryInZDirection(real shift);
    GRIDGENERATOR_EXPORT void setPeriodicShiftOnZBoundaryInXDirection(real shift);
    GRIDGENERATOR_EXPORT void setPeriodicShiftOnZBoundaryInYDirection(real shift);
    GRIDGENERATOR_EXPORT virtual void setNoSlipBoundaryCondition(SideType sideType);
    GRIDGENERATOR_EXPORT virtual void setPrecursorBoundaryCondition(SideType sideType, SPtr<FileCollection> fileCollection,
                                                                    int timeStepsBetweenReads, real velocityX = c0o1,
                                                                    real velocityY = c0o1, real velocityZ = c0o1,
                                                                    std::vector<uint> fileLevelToGridLevelMap = {});

    GRIDGENERATOR_EXPORT void setEnableFixRefinementIntoTheWall(bool enableFixRefinementIntoTheWall);

    GRIDGENERATOR_EXPORT virtual void setCommunicationProcess(int direction, uint process);

    GRIDGENERATOR_EXPORT virtual uint getCommunicationProcess(int direction) override;

    GRIDGENERATOR_EXPORT std::shared_ptr<Grid> getGrid(int level, int box);

    GRIDGENERATOR_EXPORT virtual unsigned int getNumberOfNodes(unsigned int level) const override;

    GRIDGENERATOR_EXPORT virtual uint getNumberOfFluidNodes(unsigned int level) const override;
    GRIDGENERATOR_EXPORT virtual void getFluidNodeIndices(uint* fluidNodeIndices, const int level) const override;
    GRIDGENERATOR_EXPORT virtual uint getNumberOfFluidNodesBorder(unsigned int level) const override;
    GRIDGENERATOR_EXPORT virtual void getFluidNodeIndicesBorder(uint *fluidNodeIndices, const int level) const override;

    GRIDGENERATOR_EXPORT virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords,
                                         uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative,
                                         uint *geo, const int level) const override;
    GRIDGENERATOR_EXPORT virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const override;


    GRIDGENERATOR_EXPORT uint getSlipSize(int level) const override;
    GRIDGENERATOR_EXPORT virtual void getSlipValues(real* normalX, real* normalY, real* normalZ, int* indices, int level) const override;
    GRIDGENERATOR_EXPORT virtual void getSlipQs(real* qs[27], int level) const override;

    GRIDGENERATOR_EXPORT uint getStressSize(int level) const override;
    GRIDGENERATOR_EXPORT virtual void getStressValues(  real* normalX, real* normalY, real* normalZ,
                                                        real* vx,      real* vy,      real* vz,
                                                        real* vx1,     real* vy1,     real* vz1,
                                                        int* indices, int* samplingIndices, int* samplingOffsets, real* z0, int level) const override;
    GRIDGENERATOR_EXPORT virtual void getStressQs(real* qs[27], int level) const override;

    GRIDGENERATOR_EXPORT uint getVelocitySize(int level) const override;
    GRIDGENERATOR_EXPORT virtual void getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const override;
    GRIDGENERATOR_EXPORT virtual void getVelocityQs(real* qs[27], int level) const override;

    GRIDGENERATOR_EXPORT uint getPressureSize(int level) const override;
    GRIDGENERATOR_EXPORT void getPressureValues(real* rho, int* indices, int* neighborIndices, int level) const override;
    GRIDGENERATOR_EXPORT virtual void getPressureQs(real* qs[27], int level) const override;

    GRIDGENERATOR_EXPORT uint getPrecursorSize(int level) const override;
    GRIDGENERATOR_EXPORT void getPrecursorValues(   uint* neighbor0PP, uint* neighbor0PM, uint* neighbor0MP, uint* neighbor0MM,
                                                    real* weights0PP, real* weights0PM, real* weights0MP, real* weights0MM,
                                                    int* indices, std::vector<SPtr<TransientBCInputFileReader>>& reader,
                                                    int& numberOfPrecursorNodes, size_t& numberOfQuantities, uint& timeStepsBetweenReads,
                                                    real& velocityX, real& velocityY, real& velocityZ, int level) const override;
    GRIDGENERATOR_EXPORT virtual void getPrecursorQs(real* qs[27], int level) const override;

    GRIDGENERATOR_EXPORT virtual void getGeometryQs(real *qs[27], int level) const override;
    GRIDGENERATOR_EXPORT virtual uint getGeometrySize(int level) const override;
    GRIDGENERATOR_EXPORT virtual void getGeometryIndices(int *indices, int level) const override;
    GRIDGENERATOR_EXPORT virtual bool hasGeometryValues() const override;
    GRIDGENERATOR_EXPORT virtual void getGeometryValues(real *vx, real *vy, real *vz, int level) const override;

    GRIDGENERATOR_EXPORT void writeArrows(std::string fileName) const override;

    GRIDGENERATOR_EXPORT SPtr<gg::BoundaryCondition> getBoundaryCondition( SideType side, uint level ) const override;
    GRIDGENERATOR_EXPORT SPtr<GeometryBoundaryCondition> getGeometryBoundaryCondition(uint level) const override;

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
    GRIDGENERATOR_EXPORT void getGridInformations(std::vector<int>& gridX, std::vector<int>& gridY,
                                       std::vector<int>& gridZ, std::vector<int>& distX, std::vector<int>& distY,
                                       std::vector<int>& distZ) override;
    GRIDGENERATOR_EXPORT virtual uint getNumberOfGridLevels() const override;

    GRIDGENERATOR_EXPORT uint getNumberOfNodesCF(int level) override;
    GRIDGENERATOR_EXPORT uint getNumberOfNodesFC(int level) override;

    GRIDGENERATOR_EXPORT void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const override;

    GRIDGENERATOR_EXPORT void getOffsetFC(real* xOffCf, real* yOffCf, real* zOffCf, int level) override;
    GRIDGENERATOR_EXPORT void getOffsetCF(real* xOffFc, real* yOffFc, real* zOffFc, int level) override;

    GRIDGENERATOR_EXPORT virtual uint getNumberOfSendIndices(int direction, uint level) override;
    GRIDGENERATOR_EXPORT virtual uint getNumberOfReceiveIndices(int direction, uint level) override;
    GRIDGENERATOR_EXPORT virtual void getSendIndices(int *sendIndices, int direction, int level) override;
    GRIDGENERATOR_EXPORT virtual void getReceiveIndices(int *sendIndices, int direction, int level) override;


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

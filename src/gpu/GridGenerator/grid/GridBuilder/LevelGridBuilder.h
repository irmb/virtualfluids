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

#include "gpu/GridGenerator/global.h"

#include "gpu/GridGenerator/grid/GridBuilder/GridBuilder.h"
#include "gpu/GridGenerator/grid/Grid.h"
#include "gpu/GridGenerator/grid/GridInterface.h"
#include "gpu/GridGenerator/grid/NodeValues.h"

struct Vertex;
class  Grid;
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
enum class SideType;



class LevelGridBuilder : public GridBuilder
{
protected:
    GRIDGENERATOR_EXPORT LevelGridBuilder();

public:
    GRIDGENERATOR_EXPORT static std::shared_ptr<LevelGridBuilder> makeShared();

    GRIDGENERATOR_EXPORT SPtr<Grid> getGrid(uint level) override;

    GRIDGENERATOR_EXPORT virtual ~LevelGridBuilder();

    GRIDGENERATOR_EXPORT void setSlipBoundaryCondition(SideType sideType, real nomalX, real normalY, real normalZ);
    GRIDGENERATOR_EXPORT void setStressBoundaryCondition(SideType sideType, real nomalX, real normalY, real normalZ, uint samplingOffset, real z0);
    GRIDGENERATOR_EXPORT void setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz);
    GRIDGENERATOR_EXPORT void setPressureBoundaryCondition(SideType sideType, real rho);
    GRIDGENERATOR_EXPORT void setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z);
    GRIDGENERATOR_EXPORT void setNoSlipBoundaryCondition(SideType sideType);

    GRIDGENERATOR_EXPORT void setEnableFixRefinementIntoTheWall(bool enableFixRefinementIntoTheWall);

    GRIDGENERATOR_EXPORT void setCommunicationProcess(int direction, uint process);

    GRIDGENERATOR_EXPORT uint getCommunicationProcess(int direction) override;

    GRIDGENERATOR_EXPORT virtual std::shared_ptr<Grid> getGrid(int level, int box);

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
		BoundaryConditions() {}

        std::vector<SPtr<SlipBoundaryCondition>> slipBoundaryConditions;

        std::vector<SPtr<StressBoundaryCondition>> stressBoundaryConditions;

        std::vector<SPtr<VelocityBoundaryCondition>> velocityBoundaryConditions;

        std::vector<SPtr<PressureBoundaryCondition>> pressureBoundaryConditions;

        std::vector<SPtr<VelocityBoundaryCondition>> noSlipBoundaryConditions;

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
    GRIDGENERATOR_EXPORT uint getNumberOfGridLevels() const override;

    GRIDGENERATOR_EXPORT uint getNumberOfNodesCF(int level) override;
    GRIDGENERATOR_EXPORT uint getNumberOfNodesFC(int level) override;

    GRIDGENERATOR_EXPORT void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const override;

    GRIDGENERATOR_EXPORT void getOffsetFC(real* xOffCf, real* yOffCf, real* zOffCf, int level) override;
    GRIDGENERATOR_EXPORT void getOffsetCF(real* xOffFc, real* yOffFc, real* zOffFc, int level) override;

    GRIDGENERATOR_EXPORT uint getNumberOfSendIndices(int direction, uint level) override;
    GRIDGENERATOR_EXPORT uint getNumberOfReceiveIndices(int direction, uint level) override;
    GRIDGENERATOR_EXPORT void getSendIndices(int *sendIndices, int direction, int level) override;
    GRIDGENERATOR_EXPORT void getReceiveIndices(int *sendIndices, int direction, int level) override;
};

#endif


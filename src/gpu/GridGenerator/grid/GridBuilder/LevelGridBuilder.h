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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef LEVEL_GRID_BUILDER_H
#define LEVEL_GRID_BUILDER_H

#include <vector>
#include <string>
#include <memory>
#include <array>

#include "global.h"

#include "grid/GridBuilder/GridBuilder.h"
#include "grid/Grid.h"
#include "grid/NodeValues.h"

struct Vertex;
class  Grid;
class PolyDataWriterWrapper;
class BoundingBox;
enum class Device;

class Side;
class VelocityBoundaryCondition;
enum class SideType;



class LevelGridBuilder : public GridBuilder
{
protected:
    GRIDGENERATOR_EXPORT LevelGridBuilder(Device device, const std::string& d3qxx);

public:
    GRIDGENERATOR_EXPORT static std::shared_ptr<LevelGridBuilder> makeShared(Device device, const std::string& d3qxx);

    GRIDGENERATOR_EXPORT SPtr<Grid> getGrid(uint level) override;

    GRIDGENERATOR_EXPORT virtual ~LevelGridBuilder();

    GRIDGENERATOR_EXPORT void setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz);
    GRIDGENERATOR_EXPORT void setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z);
    GRIDGENERATOR_EXPORT void setNoSlipBoundaryCondition(SideType sideType);

    GRIDGENERATOR_EXPORT virtual std::shared_ptr<Grid> getGrid(int level, int box);

    GRIDGENERATOR_EXPORT virtual unsigned int getNumberOfNodes(unsigned int level) const;


    GRIDGENERATOR_EXPORT virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords,
                                         uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, 
                                         uint *geo, const int level) const override;
    GRIDGENERATOR_EXPORT virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const;


    GRIDGENERATOR_EXPORT uint getVelocitySize(int level) const;
    GRIDGENERATOR_EXPORT virtual void getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const;
    GRIDGENERATOR_EXPORT virtual void getVelocityQs(real* qs[27], int level) const;

    GRIDGENERATOR_EXPORT SPtr<BoundaryCondition> getBoundaryCondition( SideType side, uint level ) const override;

protected:
    

    struct BoundaryConditions
    {
		BoundaryConditions() {}

        std::vector<SPtr<VelocityBoundaryCondition> > velocityBoundaryConditions;

        std::vector<SPtr<VelocityBoundaryCondition> > noSlipBoundaryConditions;
    };
    bool geometryHasValues = false;

    std::vector<std::shared_ptr<Grid> > grids;
    std::vector<SPtr<BoundaryConditions> > boundaryConditions;
    
    void checkLevel(int level);

private:
    Device device;
    std::string d3qxx;

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

};

#endif


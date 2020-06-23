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
    VF_PUBLIC LevelGridBuilder(Device device, const std::string& d3qxx);

public:
    VF_PUBLIC static std::shared_ptr<LevelGridBuilder> makeShared(Device device, const std::string& d3qxx);

    VF_PUBLIC SPtr<Grid> getGrid(uint level) override;

    VF_PUBLIC virtual ~LevelGridBuilder();

    VF_PUBLIC void setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz);
    VF_PUBLIC void setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z);
    VF_PUBLIC void setNoSlipBoundaryCondition(SideType sideType);

    VF_PUBLIC virtual std::shared_ptr<Grid> getGrid(int level, int box);

    VF_PUBLIC virtual unsigned int getNumberOfNodes(unsigned int level) const;


    VF_PUBLIC virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, 
                                         uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, 
                                         uint *geo, const int level) const override;
    VF_PUBLIC virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const;


    VF_PUBLIC uint getVelocitySize(int level) const;
    VF_PUBLIC virtual void getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const;
    VF_PUBLIC virtual void getVelocityQs(real* qs[27], int level) const;

    VF_PUBLIC SPtr<BoundaryCondition> getBoundaryCondition( SideType side, uint level ) const override;

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
    VF_PUBLIC void getGridInformations(std::vector<int>& gridX, std::vector<int>& gridY,
                                       std::vector<int>& gridZ, std::vector<int>& distX, std::vector<int>& distY,
                                       std::vector<int>& distZ) override;
    VF_PUBLIC uint getNumberOfGridLevels() const override;

    VF_PUBLIC uint getNumberOfNodesCF(int level) override;
    VF_PUBLIC uint getNumberOfNodesFC(int level) override;

    VF_PUBLIC void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const override;

    VF_PUBLIC void getOffsetFC(real* xOffCf, real* yOffCf, real* zOffCf, int level) override;
    VF_PUBLIC void getOffsetCF(real* xOffFc, real* yOffFc, real* zOffFc, int level) override;

};

#endif


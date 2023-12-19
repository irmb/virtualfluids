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
//! \file PlaneProbe.h
//! \author Henry Korb, Henrik Asmuth
//! \date 13/05/2022
//! \brief Probe computing point-wise statistics for a set of points across a plane
//!
//! The set of points can be defined by providing a list or on an x-normal plane.
//! All statistics are temporal.
//!
//=======================================================================================

#ifndef PlaneProbe_H
#define PlaneProbe_H

#include <logger/Logger.h>

#include "Probe.h"

class PlaneProbe : public Probe
{
public: 
    PlaneProbe(
        const std::string probeName,
        const std::string outputPath,
        uint tStartAvg,
        uint tAvg,
        uint tStartOut,
        uint tOut
    ): Probe(probeName, 
             outputPath,
             tStartAvg, 
             tStartAvg+1,
             tAvg,
             tStartOut, 
             tOut,
             true,
             false)
    {}

    ~PlaneProbe() = default;

    void setProbePlane(real posX, real posY, real posZ, real deltaX, real deltaY, real deltaZ)
    {
        this->posX = posX; 
        this->posY = posY; 
        this->posZ = posZ;         
        this->deltaX = deltaX; 
        this->deltaY = deltaY; 
        this->deltaZ = deltaZ; 
    }

    void getTaggedFluidNodes(GridProvider* gridProvider) override;

private:
    bool isAvailableStatistic(Statistic _variable) override;

    std::vector<PostProcessingVariable> getPostProcessingVariables(Statistic variable) override;

    void findPoints(std::vector<int>& probeIndices_level,
                    std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                    std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                    int level) override;
    void calculateQuantities(SPtr<ProbeStruct> probeStruct, uint t, int level) override;

private:
    real posX, posY, posZ;
    real deltaX, deltaY, deltaZ;
};

#endif
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

#include "Probe.h"

class PlaneProbe : public Probe
{
public: 
    PlaneProbe(
        const std::string _probeName,
        const std::string _outputPath,
        uint _tStartAvg,
        uint _tAvg,
        uint _tStartOut,
        uint _tOut
    ): Probe(_probeName, 
             _outputPath,
             _tStartAvg, 
             0,
             _tAvg,
             _tStartOut, 
             _tOut,
             true,
             false)
    {}

    void setProbePlane(real _posX, real _posY, real _posZ, real _deltaX, real _deltaY, real _deltaZ)
    {
        this->posX = _posX; 
        this->posY = _posY; 
        this->posZ = _posZ;         
        this->deltaX = _deltaX; 
        this->deltaY = _deltaY; 
        this->deltaZ = _deltaZ; 
    }

    void getInteractorFluidNodes(Parameter *para, GridProvider* gridProvider) override;

private:
    bool isAvailableStatistic(Statistic _variable) override;

    std::vector<PostProcessingVariable> getPostProcessingVariables(Statistic variable) override;

    void findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                    std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                    std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                    int level) override;
    void calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, uint t, int level) override;

private:
    real posX, posY, posZ;
    real deltaX, deltaY, deltaZ;
};

#endif
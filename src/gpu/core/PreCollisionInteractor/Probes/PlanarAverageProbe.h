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
//! \file PlanarAverageProbe.h
//! \author Henrik Asmuth
//! \date 13/05/2022
//! \brief Probe computing statistics across planes spanning the entire domain
//!
//! Computes spatial statistics across x, y or z-normal planes defined by planeNormal. 
//! The planes include all points of the domain at each respective position along that normal direction.
//! The spatial statistics can additionally be averaged in time.
//!
//=======================================================================================

#ifndef PlanarAverageProbe_H
#define PlanarAverageProbe_H

#include "Probe.h"

__global__ void moveIndicesInNegNormalDir( uint* pointIndices, uint nPoints, uint* neighborWSB, uint* neighborInplane1, uint* neighborInplane2, real* coordsX, real* coordsY, real* coordsZ ); 

__global__ void moveIndicesInPosNormalDir( uint* pointIndices, uint nPoints, uint* neighborNormal, real* coordsX, real* coordsY, real* coordsZ );

///////////////////////////////////////////////////////////////////////////////////

class PlanarAverageProbe : public Probe
{
public:
    PlanarAverageProbe(const std::string _probeName, const std::string _outputPath, uint _tStartAvg, uint _tStartTmpAvg,
                       uint _tAvg, uint _tStartOut, uint _tOut, char _planeNormal)
        : Probe(_probeName, _outputPath, _tStartAvg, _tStartTmpAvg, _tAvg, _tStartOut, _tOut, false, false),
          planeNormal(_planeNormal)
    {
        if (_tStartTmpAvg<_tStartAvg)   throw std::runtime_error("Probe: tStartTmpAvg must be larger than tStartAvg!");
        if(!(_planeNormal == 'x' || _planeNormal == 'y' || _planeNormal == 'z')) 
            throw std::runtime_error("PlanarAverageProbe: planeNormal must be 'x', 'y' or 'z'!");
    }
    ~PlanarAverageProbe() = default;

private:
    bool isAvailableStatistic(Statistic _variable) override;

    std::vector<PostProcessingVariable> getPostProcessingVariables(Statistic variable) override;

    void findPoints(std::vector<int>& probeIndices_level,
                    std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                    std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                    int level) override;
    void calculateQuantities(SPtr<ProbeStruct> probeStruct, uint t, int level) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override {};


private:
    real posX, posY, posZ;
    real deltaX, deltaY, deltaZ;
    char planeNormal;
};

#endif
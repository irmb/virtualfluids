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
//! \file WallModelProbe.h
//! \author Henrik Asmuth
//! \date 13/05/2022
//! \brief Probe computing statistics of all relevant wall model quantities used in the StressBC kernels
//!
//! Computes spatial statistics for all grid points of the StressBC 
//! The spatial statistics can additionally be averaged in time.
//!
//=======================================================================================

#ifndef WallModelProbe_H
#define WallModelProbe_H

#include "Probe.h"

///////////////////////////////////////////////////////////////////////////////////

class WallModelProbe : public Probe
{
public: 
    WallModelProbe(
        const std::string _probeName,
        const std::string _outputPath,
        uint _tStartAvg,
        uint _tStartTmpAvg,
        uint _tAvg,
        uint _tStartOut,
        uint _tOut
    ):  Probe(_probeName, 
            _outputPath,
            _tStartAvg,
            _tStartTmpAvg,
            _tAvg,
            _tStartOut, 
            _tOut,
            false,
            true)
    {
        if (_tStartTmpAvg<_tStartAvg)   throw std::runtime_error("Probe: tStartTmpAvg must be larger than tStartAvg!");
    }


    void setForceOutputToStress(bool _outputStress){ this->outputStress = _outputStress; }
    void setEvaluatePressureGradient(bool _evalPressGrad){ this->evaluatePressureGradient = _evalPressGrad; }

private:
    bool isAvailableStatistic(Statistic _variable) override;

    std::vector<PostProcessingVariable> getPostProcessingVariables(Statistic variable) override;

    void findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                    std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                    std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                    int level) override;
    void calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, uint t, int level) override;

private:
    bool outputStress = false; //!> if true, output wall force is converted to a stress 
    bool evaluatePressureGradient = false; //!> if true, mean global pressure gradient will also be evaluated
};

#endif
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
//! \file Parameters.h
//! \ingroup Parameters
//! \author Stephan Lenz
//=======================================================================================
#ifndef Parameters_H
#define Parameters_H

#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include <VirtualFluidsDefinitions.h>

enum class VF_PUBLIC ViscosityModel{
    constant,
    sutherlandsLaw
};

//! Comprises all simulation parameters for passing them to kernels
struct  VF_PUBLIC Parameters
{

    real mu = real(0.01);
    real K  = real(2.0);
    real Pr = real(1.0);
    real D  = real(0.01);
    real D1 = real(0.01);
    real D2 = real(0.01);

    real dt = real(0.01);
    real dx = real(0.01);

    Vec3 force;

    real lambdaRef = real(1.0);

    real rhoRef = real(1.0);

    ViscosityModel viscosityModel = ViscosityModel::constant;

    real boussinesqT0   = real(1.0);
    real boussinesqBeta = real(1.0);

    //////////////////////////////////////////////////////////////////////////

    bool useSmagorinsky = false;
    real smagorinskyConstant = real(0.2);

    //////////////////////////////////////////////////////////////////////////

    bool useSpongeLayer = false;
    uint spongeLayerIdx = 0;

    //////////////////////////////////////////////////////////////////////////

    uint forcingSchemeIdx = 0;

    //////////////////////////////////////////////////////////////////////////

    bool useTemperatureLimiter     = false;
    bool usePassiveScalarLimiter   = false;

    real temperatureLimiter     = real(1.0e-3);
    real passiveScalarLimiter   = real(0.1);
};

#endif

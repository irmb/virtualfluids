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
//! \file LBMSystem.h
//! \ingroup LBM
//! \author Sebastian Geller
//=======================================================================================
#ifndef LBMSYSTEM_H
#define LBMSYSTEM_H

#include <cmath>
#include <string>
#include <iostream>


//! \brief namespace for global system-functions

namespace LBMSystem
{

//#define SINGLEPRECISION

#ifdef SINGLEPRECISION
   typedef float real;
   #define REAL_CAST(x) ( (LBMSystem::real)(x) )
#else
   typedef double real;
   #define REAL_CAST(x) ( x )
#endif

   extern real SMAG_CONST;

   //////////////////////////////////////////////////////////////////////////
   //!get LBM deltaT is equal LBM DeltaX
   //!deltaT is dependent from grid level 
   //!for first grid level is deltaT = 1.0
   //!for next grid level 1/2 etc.
   static real getDeltaT(int level)
   {
      return REAL_CAST(1.0/REAL_CAST(1<<level)); 
   }

   //////////////////////////////////////////////////////////////////////////
   //!calculate collision factor omega = 1.0/(3.0*viscosity/deltaT+0.5)
   //!deltaT is dependent from grid level 
   //!for first grid level is deltaT = 1.0
   //!for next grid level 1/2 etc.
   static real calcCollisionFactor(real viscosity, int level)
   {
      //return REAL_CAST(1.0/(3.0*viscosity/deltaT+0.5));
      return REAL_CAST(1.0/(3.0*viscosity/(1.0/REAL_CAST(1<<level))+0.5));
   }

   //!bulk viscosity
   static real calcOmega2(real viscosity, int level)
   {
      return REAL_CAST(1.0/(4.5*viscosity/(1.0/REAL_CAST(1<<level))+0.5));
   }
   //!bulk viscosity
   static real calcOmega2(real viscosity, real deltaT)
   {
      return REAL_CAST(1.0/(4.5*viscosity/deltaT+0.5));
   }
}

//some typedefs for global namespace
typedef LBMSystem::real LBMReal;

#endif


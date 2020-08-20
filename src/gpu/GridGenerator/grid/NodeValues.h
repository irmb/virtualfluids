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
//! \file NodeValues.h
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef NodeValues_H
#define NodeValues_H

#define FLUID 0

#define FLUID_CFC 1
#define FLUID_CFF 2

#define FLUID_FCC 3
#define FLUID_FCF 4

#define MULTI_GPU_SEND 10
#define MULTI_GPU_RECIEVE 11

#define BC_PRESSURE 20
#define BC_VELOCITY 21
#define BC_SOLID 22

#define BC_SLIP 23
#define BC_NOSLIP 24
#define BC_OUTFLOW 25

#define STOPPER_OUT_OF_GRID 30
#define STOPPER_COARSE_UNDER_FINE 31
#define STOPPER_SOLID 32
#define STOPPER_OUT_OF_GRID_BOUNDARY 33

#define INVALID_OUT_OF_GRID 40
#define INVALID_COARSE_UNDER_FINE 41
#define INVALID_SOLID 42

#define OVERLAP_TMP 60

#endif

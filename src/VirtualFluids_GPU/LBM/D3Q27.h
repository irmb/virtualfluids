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
//! \file D3Q27.h
//! \ingroup LBM
//! \author Martin Schoenherr
//=======================================================================================
#ifndef _LB_D3Q27_H_
#define _LB_D3Q27_H_

//! \brief definitions of the 27 speeds (D3Q27)
#define dirE    0
#define dirW    1
#define dirN    2
#define dirS    3
#define dirT    4
#define dirB    5
#define dirNE   6
#define dirSW   7
#define dirSE   8
#define dirNW   9
#define dirTE   10
#define dirBW   11
#define dirBE   12
#define dirTW   13
#define dirTN   14
#define dirBS   15
#define dirBN   16
#define dirTS   17
#define dirREST 18

#define dirTNE  19
#define dirBNE  20
#define dirTSE  21
#define dirBSE  22

#define dirTNW  23
#define dirBNW  24
#define dirTSW  25
#define dirBSW  26

//! \brief definitions of start and end value (useful for loops)
#define dirSTART  0
#define dirEND   26

#endif



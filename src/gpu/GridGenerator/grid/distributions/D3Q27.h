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
//! \file D3Q27.h
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef D3Q27_H_
#define D3Q27_H_

#define DIR_27_E    0
#define DIR_27_W    1
#define DIR_27_N    2
#define DIR_27_S    3
#define DIR_27_T    4
#define DIR_27_B    5

#define DIR_27_NE   6
#define DIR_27_SW   7
#define DIR_27_SE   8
#define DIR_27_NW   9
#define DIR_27_TE   10
#define DIR_27_BW   11
#define DIR_27_BE   12
#define DIR_27_TW   13
#define DIR_27_TN   14
#define DIR_27_BS   15
#define DIR_27_BN   16
#define DIR_27_TS   17
#define DIR_27_ZERO 18

#define DIR_27_TNE  19
#define DIR_27_BNE  20
#define DIR_27_TSE  21
#define DIR_27_BSE  22

#define DIR_27_TNW  23
#define DIR_27_BNW  24
#define DIR_27_TSW  25
#define DIR_27_BSW  26

#define DIR_27_START  0
#define DIR_27_END   26


#define DIR_27_E_X  1
#define DIR_27_E_Y  0
#define DIR_27_E_Z  0

#define DIR_27_W_X  -1
#define DIR_27_W_Y  0
#define DIR_27_W_Z  0

#define DIR_27_N_X  0
#define DIR_27_N_Y  1
#define DIR_27_N_Z  0

#define DIR_27_S_X  0
#define DIR_27_S_Y  -1
#define DIR_27_S_Z  0

#define DIR_27_T_X  0
#define DIR_27_T_Y  0
#define DIR_27_T_Z  1

#define DIR_27_B_X  0
#define DIR_27_B_Y  0
#define DIR_27_B_Z  -1


#define DIR_27_NE_X  1
#define DIR_27_NE_Y  1
#define DIR_27_NE_Z  0

#define DIR_27_SW_X  -1
#define DIR_27_SW_Y  -1
#define DIR_27_SW_Z  0

#define DIR_27_SE_X  1
#define DIR_27_SE_Y  -1
#define DIR_27_SE_Z  0

#define DIR_27_NW_X  -1
#define DIR_27_NW_Y  1
#define DIR_27_NW_Z  0

#define DIR_27_TE_X  1
#define DIR_27_TE_Y  0
#define DIR_27_TE_Z  1

#define DIR_27_BW_X  -1
#define DIR_27_BW_Y  0
#define DIR_27_BW_Z  -1

#define DIR_27_BE_X  1
#define DIR_27_BE_Y  0
#define DIR_27_BE_Z  -1

#define DIR_27_TW_X  -1
#define DIR_27_TW_Y  0
#define DIR_27_TW_Z  1

#define DIR_27_TN_X  0
#define DIR_27_TN_Y  1
#define DIR_27_TN_Z  1

#define DIR_27_BS_X  0
#define DIR_27_BS_Y  -1
#define DIR_27_BS_Z  -1

#define DIR_27_BN_X  0
#define DIR_27_BN_Y  1
#define DIR_27_BN_Z  -1

#define DIR_27_TS_X  0
#define DIR_27_TS_Y  -1
#define DIR_27_TS_Z  1



#define DIR_27_TNE_X  1
#define DIR_27_TNE_Y  1
#define DIR_27_TNE_Z  1

#define DIR_27_BNE_X  1
#define DIR_27_BNE_Y  1
#define DIR_27_BNE_Z  -1

#define DIR_27_TSE_X  1
#define DIR_27_TSE_Y  -1
#define DIR_27_TSE_Z  1

#define DIR_27_BSE_X  1
#define DIR_27_BSE_Y  -1
#define DIR_27_BSE_Z  -1


#define DIR_27_TNW_X  -1
#define DIR_27_TNW_Y  1
#define DIR_27_TNW_Z  1

#define DIR_27_BNW_X  -1
#define DIR_27_BNW_Y  1
#define DIR_27_BNW_Z  -1

#define DIR_27_TSW_X  -1
#define DIR_27_TSW_Y  -1
#define DIR_27_TSW_Z  1

#define DIR_27_BSW_X  -1
#define DIR_27_BSW_Y  -1
#define DIR_27_BSW_Z  -1

#endif

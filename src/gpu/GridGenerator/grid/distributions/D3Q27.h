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

static constexpr int STARTDIR = 0;
static constexpr int ENDDIR = 26;

static constexpr int E = 0;
static constexpr int W = 1;
static constexpr int N = 2;
static constexpr int S = 3;
static constexpr int T = 4;
static constexpr int B = 5;
								  
static constexpr int NE = 6;
static constexpr int SW = 7;
static constexpr int SE = 8;
static constexpr int NW = 9;
static constexpr int TE = 10;
static constexpr int BW = 11;
static constexpr int BE = 12;
static constexpr int TW = 13;
static constexpr int TN = 14;
static constexpr int BS = 15;
static constexpr int BN = 16;
static constexpr int TS = 17;

static constexpr int  TNE  = 18;
static constexpr int  TNW  = 19;
static constexpr int  TSE  = 20;
static constexpr int  TSW  = 21;
static constexpr int  BNE  = 22;
static constexpr int  BNW  = 23;
static constexpr int  BSE  = 24;
static constexpr int  BSW  = 25;
								  
static constexpr int REST = 26;




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

#define DIR_27_REST_X  0
#define DIR_27_REST_Y  0
#define DIR_27_REST_Z  0

#endif

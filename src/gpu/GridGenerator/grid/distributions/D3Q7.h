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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef D3Q7_H
#define D3Q7_H

#define DIR_7_E    0
#define DIR_7_W    1
#define DIR_7_N    2
#define DIR_7_S    3
#define DIR_7_T    4
#define DIR_7_B    5
#define DIR_7_ZERO 6

#define DIR_7_START  0
#define DIR_7_END   6


#define DIR_7_E_X  1
#define DIR_7_E_Y  0
#define DIR_7_E_Z  0

#define DIR_7_W_X  -1
#define DIR_7_W_Y  0
#define DIR_7_W_Z  0

#define DIR_7_N_X  0
#define DIR_7_N_Y  1
#define DIR_7_N_Z  0

#define DIR_7_S_X  0
#define DIR_7_S_Y  -1
#define DIR_7_S_Z  0

#define DIR_7_T_X  0
#define DIR_7_T_Y  0
#define DIR_7_T_Z  1

#define DIR_7_B_X  0
#define DIR_7_B_Y  0
#define DIR_7_B_Z  -1


#endif

//! \}

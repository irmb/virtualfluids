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
//! \author Soeren Peters, Stephan Lenz, Martin Schoenherr
//=======================================================================================
#ifndef NodeValues_H
#define NodeValues_H

namespace vf::gpu
{

static constexpr char FLUID = 0;

static constexpr char FLUID_CFC = 1;
static constexpr char FLUID_CFF = 2;

static constexpr char FLUID_FCC = 3;
static constexpr char FLUID_FCF = 4;

static constexpr char MULTI_GPU_SEND    = 10;
static constexpr char MULTI_GPU_RECIEVE = 11;

static constexpr char BC_PRESSURE = 20;
static constexpr char BC_VELOCITY = 21;
static constexpr char BC_SOLID    = 22;

static constexpr char BC_SLIP    = 23;
static constexpr char BC_NOSLIP  = 24;
static constexpr char BC_OUTFLOW = 25;
static constexpr char BC_STRESS   = 26;

static constexpr char STOPPER_OUT_OF_GRID          = 30;
static constexpr char STOPPER_COARSE_UNDER_FINE    = 31;
static constexpr char STOPPER_SOLID                = 32;
static constexpr char STOPPER_OUT_OF_GRID_BOUNDARY = 33;

static constexpr char INVALID_OUT_OF_GRID       = 40;
static constexpr char INVALID_COARSE_UNDER_FINE = 41;
static constexpr char INVALID_SOLID             = 42;

static constexpr char INSIDE                    = 50;
static constexpr char NEGATIVE_DIRECTION_BORDER = 51;
static constexpr char Q_DEPRECATED              = 52;

static constexpr char OVERLAP_TMP = 60;

}

#endif

//! \}

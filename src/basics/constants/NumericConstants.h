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
//! \addtogroup constants
//! \ingroup basics
//! \{
//=======================================================================================
#ifndef BASICS_NUMERIC_CONSTANT_H
#define BASICS_NUMERIC_CONSTANT_H

#ifndef __CUDACC__
#include <cmath>
#endif

namespace vf::basics::constant
{

#ifdef VF_DOUBLE_ACCURACY
static constexpr double c1o2 = 1. / 2.;
static constexpr double c3o2 = 3. / 2.;
static constexpr double c1o3 = 1. / 3.;
static constexpr double c2o3 = 2. / 3.;
static constexpr double c1o4 = 1. / 4.;
static constexpr double c3o4 = 3. / 4.;
static constexpr double c1o6 = 1. / 6.;
static constexpr double c1o7 = 1. / 7.;
static constexpr double c1o8 = 1. / 8.;
static constexpr double c1o9 = 1. / 9.;
static constexpr double c2o9 = 2. / 9.;
static constexpr double c4o9 = 4. / 9.;
static constexpr double c4o10 = 4. / 10.;
static constexpr double c1o10 = 1. / 10.;
static constexpr double c1o12 = 1. / 12.;
static constexpr double c1o16 = 1. / 16.;
static constexpr double c3o16 = 3. / 16.;
static constexpr double c9o16 = 9. / 16.;
static constexpr double c1o18 = 1. / 18.;
static constexpr double c1o20 = 1. / 20.;
static constexpr double c19o20 = 19. / 20.;
static constexpr double c21o20 = 21. / 20.;
static constexpr double c1o24 = 1. / 24.;
static constexpr double c1o27 = 1. / 27.;
static constexpr double c3o32 = 3. / 32.;
static constexpr double c4o32 = 4. / 32.;
static constexpr double c1o36 = 1. / 36.;
static constexpr double c1o48 = 1. / 48.;
static constexpr double c1o64 = 1. / 64.;
static constexpr double c3o64 = 3. / 64.;
static constexpr double c9o64 = 9. / 64.;
static constexpr double c27o64 = 27. / 64.;
static constexpr double c1o66 = 1. / 66.;
static constexpr double c1o72 = 1. / 72.;
static constexpr double c1o264 = 1. / 264.;
static constexpr double c8o27 = 8. / 27.;
static constexpr double c2o27 = 2. / 27.;
static constexpr double c1o54 = 1. / 54.;
static constexpr double c1o100 = 1. / 100.;
static constexpr double c99o100 = 99. / 100;
static constexpr double c1o126 = 1. / 126.;
static constexpr double c1o216 = 1. / 216.;
static constexpr double c5o4 = 5. / 4.;
static constexpr double c4o3 = 4. / 3.;
static constexpr double c9o4 = 9. / 4.;
static constexpr double c5o2 = 5. / 2.;
static constexpr double c9o2 = 9. / 2.;

static constexpr double c0o1 = 0.;
static constexpr double c1o1 = 1.;
static constexpr double c2o1 = 2.;
static constexpr double c3o1 = 3.;
static constexpr double c4o1 = 4.;
static constexpr double c5o1 = 5.;
static constexpr double c6o1 = 6.;
static constexpr double c7o1 = 7.;
static constexpr double c8o1 = 8.;
static constexpr double c9o1 = 9.;
static constexpr double c10o1 = 10.;
static constexpr double c11o1 = 11.;
static constexpr double c12o1 = 12.;
static constexpr double c13o1 = 13.;
static constexpr double c14o1 = 14.;
static constexpr double c15o1 = 15.;
static constexpr double c16o1 = 16.;
static constexpr double c17o1 = 17.;
static constexpr double c18o1 = 18.;
static constexpr double c21o1 = 21.;
static constexpr double c24o1 = 24.;
static constexpr double c25o1 = 25.;
static constexpr double c26o1 = 26.;
static constexpr double c27o1 = 27.;
static constexpr double c28o1 = 28.;
static constexpr double c29o1 = 29.;
static constexpr double c30o1 = 30.;
static constexpr double c32o1 = 32.;
static constexpr double c33o1 = 33.;
static constexpr double c34o1 = 34.;
static constexpr double c36o1 = 36.;
static constexpr double c40o1 = 40.;
static constexpr double c42o1 = 42.;
static constexpr double c46o1 = 46.;
static constexpr double c48o1 = 48.;
static constexpr double c50o1 = 50.;
static constexpr double c52o1 = 52.;
static constexpr double c54o1 = 54.;
static constexpr double c56o1 = 56.;
static constexpr double c64o1 = 64.;
static constexpr double c66o1 = 66.;
static constexpr double c68o1 = 68.;
static constexpr double c69o1 = 69.;
static constexpr double c72o1 = 72.;
static constexpr double c84o1 = 84.;
static constexpr double c88o1 = 88.;
static constexpr double c96o1 = 96.;
static constexpr double c100o1 = 100.;
static constexpr double c130o1 = 130.;
static constexpr double c152o1 = 152.;
static constexpr double c166o1 = 166.;
static constexpr double c195o1 = 195.;
static constexpr double c216o1 = 216.;
static constexpr double c264o1 = 264.;
static constexpr double c290o1 = 290.;
static constexpr double c367o1 = 367.;

static constexpr double c0p0000002 = 0.0000002;
static constexpr double c10eM30 = 1e-30;
static constexpr double c10eM10 = 1e-10;
static constexpr double cSmallSingle = 0.0000000002;

#ifndef __CUDACC__
static const double cPi = 4.0 * std::atan(1.0);               // 3.1415926535
static const double c2Pi = 8.0 * std::atan(1.0);              // 6.2831853071
static const double cPio180 = 4.0 * std::atan(1.0) / 180.0;   // 1.74532925199e-2
static const double c180oPi = 180.0 / (4.0 * std::atan(1.0)); // 57.2957795131
#else
static constexpr double cPi = 3.1415926535;
static constexpr double c2Pi = 6.28318530717;
static constexpr double cPio180 = 1.74532925199e-2;
static constexpr double c180oPi = 57.2957795131;
#endif

static const double c1oSqrt2 = 1.0 / sqrt(2.0); // 0.707106781
static const double c1oSqrt3 = 1.0 / sqrt(3.0); // 0.577350269
static const double cSqrt2 = sqrt(2.0);         // 1.4142135
static const double cSqrt3 = sqrt(3.0);         // 1.7320508

#else
static constexpr float c1o2 = 1.0f / 2.0f;
static constexpr float c3o2 = 3.0f / 2.0f;
static constexpr float c1o3 = 1.0f / 3.0f;
static constexpr float c2o3 = 2.0f / 3.0f;
static constexpr float c1o4 = 1.0f / 4.0f;
static constexpr float c3o4 = 3.0f / 4.0f;
static constexpr float c1o6 = 1.0f / 6.0f;
static constexpr float c1o7 = 1.0f / 7.0f;
static constexpr float c1o8 = 1.0f / 8.0f;
static constexpr float c1o9 = 1.0f / 9.0f;
static constexpr float c2o9 = 2.0f / 9.0f;
static constexpr float c4o9 = 4.0f / 9.0f;
static constexpr float c4o10 = 4.0f / 10.0f;
static constexpr float c1o10 = 1.0f / 10.0f;
static constexpr float c1o12 = 1.0f / 12.0f;
static constexpr float c1o16 = 1.0f / 16.0f;
static constexpr float c3o16 = 3.0f / 16.0f;
static constexpr float c9o16 = 9.0f / 16.0f;
static constexpr float c1o18 = 1.0f / 18.0f;
static constexpr float c1o20 = 1.0f / 20.0f;
static constexpr float c19o20 = 19.0f / 20.0f;
static constexpr float c21o20 = 21.0f / 20.0f;
static constexpr float c1o24 = 1.0f / 24.0f;
static constexpr float c1o27 = 1.0f / 27.0f;
static constexpr float c3o32 = 3.0f / 32.0f;
static constexpr float c4o32 = 4.0f / 32.0f;
static constexpr float c1o36 = 1.0f / 36.0f;
static constexpr float c1o48 = 1.0f / 48.0f;
static constexpr float c1o64 = 1.0f / 64.0f;
static constexpr float c3o64 = 3.0f / 64.0f;
static constexpr float c9o64 = 9.0f / 64.0f;
static constexpr float c27o64 = 27.0f / 64.0f;
static constexpr float c1o66 = 1.0f / 66.0f;
static constexpr float c1o72 = 1.0f / 72.0f;
static constexpr float c1o264 = 1.0f / 264.0f;
static constexpr float c8o27 = 8.0f / 27.0f;
static constexpr float c2o27 = 2.0f / 27.0f;
static constexpr float c1o54 = 1.0f / 54.0f;
static constexpr float c1o100 = 1.0f / 100.0f;
static constexpr float c99o100 = 99.0f / 100.0f;
static constexpr float c1o126 = 1.0f / 126.0f;
static constexpr float c1o216 = 1.0f / 216.0f;
static constexpr float c5o4 = 5.0f / 4.0f;
static constexpr float c4o3 = 4.0f / 3.0f;
static constexpr float c9o4 = 9.0f / 4.0f;
static constexpr float c5o2 = 5.0f / 2.0f;
static constexpr float c9o2 = 9.0f / 2.0f;

static constexpr float c0o1 = 0.f;
static constexpr float c1o1 = 1.f;
static constexpr float c2o1 = 2.f;
static constexpr float c3o1 = 3.f;
static constexpr float c4o1 = 4.f;
static constexpr float c5o1 = 5.f;
static constexpr float c6o1 = 6.f;
static constexpr float c7o1 = 7.f;
static constexpr float c8o1 = 8.f;
static constexpr float c9o1 = 9.f;
static constexpr float c10o1 = 10.f;
static constexpr float c11o1 = 11.f;
static constexpr float c12o1 = 12.f;
static constexpr float c13o1 = 13.f;
static constexpr float c14o1 = 14.f;
static constexpr float c15o1 = 15.f;
static constexpr float c16o1 = 16.f;
static constexpr float c17o1 = 17.f;
static constexpr float c18o1 = 18.f;
static constexpr float c21o1 = 21.f;
static constexpr float c24o1 = 24.f;
static constexpr float c25o1 = 25.f;
static constexpr float c26o1 = 26.f;
static constexpr float c27o1 = 27.f;
static constexpr float c28o1 = 28.f;
static constexpr float c29o1 = 29.f;
static constexpr float c30o1 = 30.f;
static constexpr float c32o1 = 32.f;
static constexpr float c33o1 = 33.f;
static constexpr float c34o1 = 34.f;
static constexpr float c36o1 = 36.f;
static constexpr float c40o1 = 40.f;
static constexpr float c42o1 = 42.f;
static constexpr float c46o1 = 46.f;
static constexpr float c48o1 = 48.f;
static constexpr float c50o1 = 50.f;
static constexpr float c52o1 = 52.f;
static constexpr float c54o1 = 54.f;
static constexpr float c56o1 = 56.f;
static constexpr float c64o1 = 64.f;
static constexpr float c66o1 = 66.f;
static constexpr float c68o1 = 68.f;
static constexpr float c69o1 = 69.f;
static constexpr float c72o1 = 72.f;
static constexpr float c84o1 = 84.f;
static constexpr float c88o1 = 88.f;
static constexpr float c96o1 = 96.f;
static constexpr float c100o1 = 100.0f;
static constexpr float c130o1 = 130.0f;
static constexpr float c152o1 = 152.0f;
static constexpr float c166o1 = 166.0f;
static constexpr float c195o1 = 195.0f;
static constexpr float c216o1 = 216.0f;
static constexpr float c264o1 = 264.0f;
static constexpr float c290o1 = 290.0f;
static constexpr float c367o1 = 367.0f;

static constexpr float c0p0000002 = 0.0000002f;
static constexpr float c10eM30 = 1e-30f;
static constexpr float c10eM10 = 1e-10f;
static constexpr float cSmallSingle = 0.0000000002f;

#ifndef __CUDACC__
static const float cPi = 4.0f * std::atan(1.0f);                // 3.1415926535
static const float c2Pi = 8.0f * std::atan(1.0f);               // 6.2831853071
static const float cPio180 = 4.0f * std::atan(1.0f) / 180.0f;   // 1.74532925199e-2
static const float c180oPi = 180.0f / (4.0f * std::atan(1.0f)); // 57.2957795131
#else
static constexpr float cPi = 3.1415926535f;
static constexpr float c2Pi = 6.28318530717f;
static constexpr float cPio180 = 1.74532925199e-2f;
static constexpr float c180oPi = 57.2957795131f;
#endif

static const float c1oSqrt2 = 1.0 / sqrtf(2.0); // 0.707106781
static const float c1oSqrt3 = 1.0 / sqrtf(3.0); // 0.577350269
static const float cSqrt2 = sqrtf(2.0);         // 1.4142135
static const float cSqrt3 = sqrtf(3.0);         // 1.7320508

#endif

}

#endif

//! \}

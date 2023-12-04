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
//! \file D3Q27System.cpp
//! \ingroup LBM
//! \author Konstantin Kutscher, Sebastian Geller, Soeren Freudiger
//=======================================================================================
#include "D3Q27System.h"

#include "lbm/MacroscopicQuantities.h"

namespace D3Q27System
{
//using namespace UbMath;
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;

// index            0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26
const int DX1[] = { 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1 };
const int DX2[] = { 0,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1, -1 };
const int DX3[] = { 0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1,  1, -1, -1, -1, -1 };

const real WEIGTH[] = { c8o27,  
                          c2o27,  c2o27,  c2o27,  c2o27,  c2o27,  c2o27,  
                          c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,
                          c1o216, c1o216, c1o216, c1o216, c1o216, c1o216, c1o216, c1o216 };

const int INVDIR[] = { d000, iP00, iM00, i0P0, i0M0, i00P, i00M, 
                                iPP0, iMM0, iPM0, iMP0, iP0P, iM0M, iP0M, iM0P, i0PP, i0MM, i0PM, i0MP,
                                iPPP, iMPP, iPMP, iMMP, iPPM, iMPM, iPMM, iMMM };

//////////////////////////////////////////////////////////////////////////



real getDensity(const real *const &f /*[27]*/)
{
    return vf::lbm::getDensity(f);
}

real getIncompVelocityX1(const real *const &f /*[27]*/)
{
    return vf::lbm::getIncompressibleVelocityX1(f);
}

real getIncompVelocityX2(const real *const &f /*[27]*/)
{
    return vf::lbm::getIncompressibleVelocityX2(f);
}

real getIncompVelocityX3(const real *const &f /*[27]*/)
{
    return vf::lbm::getIncompressibleVelocityX3(f);
}





} // namespace D3Q27System
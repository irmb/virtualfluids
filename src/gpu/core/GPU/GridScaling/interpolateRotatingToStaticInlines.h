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
//! \file scaleCF_compressible.cu
//! \ingroup GPU/GridScaling
//! \author Anna Wellmann, Martin Schoenherr
//=======================================================================================

#ifndef INTERPOLATE_ROTATING_TO_STATIC_INLINES
#define INTERPOLATE_ROTATING_TO_STATIC_INLINES

#include "DataTypes.h"
#include "basics/constants/NumericConstants.h"

using namespace vf::basics::constant;

__inline__ __host__ __device__ void rotateSecondOrderMomentsGlobalToRotating(real &m011, real &m101, real &m110,
                                                                             real &mxxMyy, real &mxxMzz, real angleX,
                                                                             real angleY, real angleZ)
{
    real mxxMyyRotated = mxxMyy;
    real mxxMzzRotated = mxxMzz;
    real m011Rotated = m011;
    real m101Rotated = m101;
    real m110Rotated = m110;

    if (angleX != c0o1) {
        mxxMyyRotated = (mxxMyy + mxxMzz + (mxxMyy - mxxMzz) * cos(angleX * c2o1) - c2o1 * m011 * sin(angleX * c2o1)) * c1o2;
        mxxMzzRotated =
            (mxxMyy + mxxMzz + (-mxxMyy + mxxMzz) * cos(angleX * c2o1) + c2o1 * m011 * sin(angleX * c2o1)) * c1o2;

        m011Rotated = m011 * cos(angleX * c2o1) + (mxxMyy - mxxMzz) * cos(angleX) * sin(angleX);
        m101Rotated = m101 * cos(angleX) - m110 * sin(angleX);
        m110Rotated = m110 * cos(angleX) + m101 * sin(angleX);
    } else if (angleY != c0o1) {
        mxxMyyRotated = mxxMyy - mxxMzz * c1o2 + c1o2 * mxxMzz * cos(angleY * c2o1) - m101 * sin(angleY * c2o1);
        mxxMzzRotated = mxxMzz * cos(angleY * c2o1) - 2 * m101 * sin(angleY * c2o1);

        m011Rotated = m011 * cos(angleY) + m110 * sin(angleY);
        m101Rotated = m101 * cos(angleY * c2o1) + mxxMzz * cos(angleY) * sin(angleY);
        m110Rotated = m110 * cos(angleY) - m011 * sin(angleY);
    } else if (angleZ != c0o1) {
        mxxMyyRotated = mxxMyy * cos(angleZ * c2o1) + c2o1 * m110 * sin(angleZ * c2o1);
        mxxMzzRotated = -c1o2 * mxxMyy + mxxMzz + c1o2 * mxxMyy * cos(angleZ * c2o1) + m110 * sin(angleZ * c2o1);

        m011Rotated = m011 * cos(angleZ) - m101 * sin(angleZ);
        m101Rotated = m101 * cos(angleZ) + m011 * sin(angleZ);
        m110Rotated = m110 * cos(angleZ * c2o1) - mxxMyy * cos(angleZ) * sin(angleZ);
    }

    mxxMyy = mxxMyyRotated;
    mxxMzz = mxxMzzRotated;
    m011 = m011Rotated;
    m101 = m101Rotated;
    m110 = m110Rotated;
}

#endif

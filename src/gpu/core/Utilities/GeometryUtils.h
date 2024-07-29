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
//! \addtogroup gpu_Utilities Utilities
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================

#ifndef _GEOMETRYUTILS_H
#define _GEOMETRYUTILS_H

#include <cuda_runtime.h>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

__inline__ __host__ __device__ void getNeighborIndicesOfBSW(uint k_MMM, uint& k_PMM, uint& k_MPM, uint& k_MMP, uint& k_PPM,
                                                            uint& k_PMP, uint& k_MPP, uint& k_PPP, const uint* neighborX,
                                                            const uint* neighborY, const uint* neighborZ)
{
    k_PMM = neighborX[k_MMM];
    k_MPM = neighborY[k_MMM];
    k_MMP = neighborZ[k_MMM];
    k_PPM = neighborY[k_PMM];
    k_PMP = neighborZ[k_PMM];
    k_MPP = neighborZ[k_MPM];
    k_PPP = neighborX[k_MPP];
}

__inline__ __host__ __device__ uint findNearestCellBSW(const uint index, const real* coordsX, const real* coordsY,
                                                       const real* coordsZ, const real posX, const real posY,
                                                       const real posZ, const uint* neighborsX, const uint* neighborsY,
                                                       const uint* neighborsZ, const uint* neighborsWSB)
{
    uint new_index = index;

    while (coordsX[new_index] > posX && coordsY[new_index] > posY && coordsZ[new_index] > posZ) {
        new_index = max(1U, neighborsWSB[new_index]);
    }

    while (coordsX[new_index] > posX && coordsY[new_index] > posY) {
        new_index = max(1, neighborsZ[neighborsWSB[new_index]]);
    }
    while (coordsX[new_index] > posX && coordsZ[new_index] > posZ) {
        new_index = max(1, neighborsY[neighborsWSB[new_index]]);
    }
    while (coordsY[new_index] > posY && coordsZ[new_index] > posZ) {
        new_index = max(1, neighborsX[neighborsWSB[new_index]]);
    }

    while (coordsX[new_index] > posX) {
        new_index = max(1U, neighborsY[neighborsZ[neighborsWSB[new_index]]]);
    }
    while (coordsY[new_index] > posY) {
        new_index = max(1U, neighborsX[neighborsZ[neighborsWSB[new_index]]]);
    }
    while (coordsZ[new_index] > posZ) {
        new_index = max(1U, neighborsX[neighborsY[neighborsWSB[new_index]]]);
    }

    while (coordsX[new_index] < posX) {
        new_index = max(1U, neighborsX[new_index]);
    }
    while (coordsY[new_index] < posY) {
        new_index = max(1U, neighborsY[new_index]);
    }
    while (coordsZ[new_index] < posZ) {
        new_index = max(1U, neighborsZ[new_index]);
    }

    return neighborsWSB[new_index];
}

__inline__ __host__ __device__ real trilinearInterpolation(real dXM, real dYM, real dZM, uint kMMM, uint kPMM, uint kMPM,
                                                           uint kMMP, uint kPPM, uint kPMP, uint kMPP, uint kPPP,
                                                           const real* quantity)
{
    const real dXP = vf::basics::constant::c1o1 - dXM;
    const real dYP = vf::basics::constant::c1o1 - dYM;
    const real dZP = vf::basics::constant::c1o1 - dZM;
    return (dXP * dYP * dZP * quantity[kMMM] + dXM * dYP * dZP * quantity[kPMM] + dXP * dYM * dZP * quantity[kMPM] +
            dXM * dYM * dZP * quantity[kPPM] + dXP * dYP * dZM * quantity[kMMP] + dXM * dYP * dZM * quantity[kPMP] +
            dXP * dYM * dZM * quantity[kMPP] + dXM * dYM * dZM * quantity[kPPP]);
}

__inline__ __host__ __device__ void translate2D(real posX, real posY, real& newPosX, real& newPosY, real translationX,
                                                real translationY)
{
    newPosX = posX + translationX;
    newPosY = posY + translationY;
}

__inline__ __host__ __device__ void invTranslate2D(real posX, real posY, real& newPosX, real& newPosY, real translationX,
                                                   real translationY)
{
    newPosX = posX - translationX;
    newPosY = posY - translationY;
}

__inline__ __host__ __device__ void translate3D(real posX, real posY, real posZ, real& newPosX, real& newPosY, real& newPosZ,
                                                real translationX, real translationY, real translationZ)
{
    newPosX = posX + translationX;
    newPosY = posY + translationY;
    newPosZ = posZ + translationZ;
}

__inline__ __host__ __device__ void invTranslate3D(real posX, real posY, real posZ, real& newPosX, real& newPosY,
                                                   real& newPosZ, real translationX, real translationY, real translationZ)
{
    newPosX = posX - translationX;
    newPosY = posY - translationY;
    newPosZ = posZ - translationZ;
}

__inline__ __host__ __device__ void rotate2D(real angle, real posX, real posY, real& newPosX, real& newPosY)
{
    newPosX = posX * cos(angle) - posY * sin(angle);
    newPosY = posX * sin(angle) + posY * cos(angle);
}

__inline__ __host__ __device__ void rotate2D(real angle, real posX, real posY, real& newPosX, real& newPosY, real originX,
                                             real originY)
{
    real tmpX, tmpY;
    invTranslate2D(posX, posY, newPosX, newPosY, originX, originY);
    rotate2D(angle, newPosX, newPosY, tmpX, tmpY);
    translate2D(tmpX, tmpY, newPosX, newPosY, originX, originY);
}

__inline__ __host__ __device__ void invRotate2D(real angle, real posX, real posY, real& newPosX, real& newPosY)
{
    newPosX = posX * cos(angle) + posY * sin(angle);
    newPosY = -posX * sin(angle) + posY * cos(angle);
}

__inline__ __host__ __device__ void invRotate2D(real angle, real posX, real posY, real& newPosX, real& newPosY, real originX,
                                                real originY)
{
    real tmpX, tmpY;
    invTranslate2D(posX, posY, newPosX, newPosY, originX, originY);
    invRotate2D(angle, newPosX, newPosY, tmpX, tmpY);
    translate2D(tmpX, tmpY, newPosX, newPosY, originX, originY);
}

__inline__ __host__ __device__ void rotateAboutX3D(real angle, real posX, real posY, real posZ, real& newPosX, real& newPosY,
                                                   real& newPosZ)
{
    newPosX = posX;
    rotate2D(angle, posY, posZ, newPosY, newPosZ);
}

__inline__ __host__ __device__ void rotateAboutX3D(real angle, real posX, real posY, real posZ, real& newPosX, real& newPosY,
                                                   real& newPosZ, real originX, real originY, real originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    rotateAboutX3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

__inline__ __host__ __device__ void invRotateAboutX3D(real angle, real posX, real posY, real posZ, real& newPosX,
                                                      real& newPosY, real& newPosZ)
{
    newPosX = posX;
    invRotate2D(angle, posY, posZ, newPosY, newPosZ);
}

__inline__ __host__ __device__ void invRotateAboutX3D(real angle, real posX, real posY, real posZ, real& newPosX,
                                                      real& newPosY, real& newPosZ, real originX, real originY, real originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    invRotateAboutX3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

__inline__ __host__ __device__ void rotateAboutY3D(real angle, real posX, real posY, real posZ, real& newPosX, real& newPosY,
                                                   real& newPosZ)
{
    newPosY = posY;
    rotate2D(angle, posX, posZ, newPosX, newPosZ);
}

__inline__ __host__ __device__ void rotateAboutY3D(real angle, real posX, real posY, real posZ, real& newPosX, real& newPosY,
                                                   real& newPosZ, real originX, real originY, real originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    rotateAboutY3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

__inline__ __host__ __device__ void invRotateAboutY3D(real angle, real posX, real posY, real posZ, real& newPosX,
                                                      real& newPosY, real& newPosZ)
{
    newPosY = posY;
    invRotate2D(angle, posX, posZ, newPosX, newPosZ);
}

__inline__ __host__ __device__ void invRotateAboutY3D(real angle, real posX, real posY, real posZ, real& newPosX,
                                                      real& newPosY, real& newPosZ, real originX, real originY, real originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    invRotateAboutY3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

__inline__ __host__ __device__ void rotateAboutZ3D(real angle, real posX, real posY, real posZ, real& newPosX, real& newPosY,
                                                   real& newPosZ)
{
    newPosZ = posZ;
    rotate2D(angle, posX, posY, newPosX, newPosY);
}

__inline__ __host__ __device__ void rotateAboutZ3D(real angle, real posX, real posY, real posZ, real& newPosX, real& newPosY,
                                                   real& newPosZ, real originX, real originY, real originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    rotateAboutZ3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

__inline__ __host__ __device__ void invRotateAboutZ3D(real angle, real posX, real posY, real posZ, real& newPosX,
                                                      real& newPosY, real& newPosZ)
{
    newPosZ = posZ;
    invRotate2D(angle, posX, posY, newPosX, newPosY);
}

__inline__ __host__ __device__ void invRotateAboutZ3D(real angle, real posX, real posY, real posZ, real& newPosX,
                                                      real& newPosY, real& newPosZ, real originX, real originY, real originZ)
{
    real tmpX, tmpY, tmpZ;
    invTranslate3D(posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    invRotateAboutZ3D(angle, newPosX, newPosY, newPosZ, tmpX, tmpY, tmpZ);
    translate3D(tmpX, tmpY, tmpZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
}

#endif

//! \}

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
//! \addtogroup interpolation
//! \ingroup lbm
//! \{
//=======================================================================================
#ifndef LBM_SCALING_HELPER_FUNCTIONS_H
#define LBM_SCALING_HELPER_FUNCTIONS_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/constants/NumericConstants.h>

#include "lbm/constants/D3Q27.h"

#include "lbm/MacroscopicQuantities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

namespace vf::lbm
{

// The Coefficients struct needs be created like this:
// MomentsOnSourceNodeSet momentsSet;
// momentsSet.calculatePPP(f, omega);
// ... and so on
// InterpolationCoefficients coeffs;
// momentsSet.calculateCoefficients(coeffs);

// Coefficients of the interpolation polynomial
struct InterpolationCoefficients
{
    real a000, a100, a010, a001, a200, a020, a002, a110, a101, a011;
    real b000, b100, b010, b001, b200, b020, b002, b110, b101, b011;
    real c000, c100, c010, c001, c200, c020, c002, c110, c101, c011;
    real d000, d100, d010, d001, d110, d101, d011;
    real a111, b111, c111, d111;
    real LaplaceRho;
};

// Private struct - is only used within the MomentsOnSourceNodeSet
struct MomentsOnSourceNode
{
    real drho;
    real velocityX;
    real velocityY;
    real velocityZ;
    real kxyFromfcNEQ;
    real kyzFromfcNEQ;
    real kxzFromfcNEQ;
    real kxxMyyFromfcNEQ;
    real kxxMzzFromfcNEQ;

    __host__ __device__ void calculate(const real *const f, const real omega)
    {
        // const real f_000 = f[dir::d000];
        const real fP00 = f[dir::dP00];
        const real fM00 = f[dir::dM00];
        const real f0P0 = f[dir::d0P0];
        const real f0M0 = f[dir::d0M0];
        const real f00P = f[dir::d00P];
        const real f00M = f[dir::d00M];
        const real fPP0 = f[dir::dPP0];
        const real fMM0 = f[dir::dMM0];
        const real fPM0 = f[dir::dPM0];
        const real fMP0 = f[dir::dMP0];
        const real fP0P = f[dir::dP0P];
        const real fM0M = f[dir::dM0M];
        const real fP0M = f[dir::dP0M];
        const real fM0P = f[dir::dM0P];
        const real f0PP = f[dir::d0PP];
        const real f0MM = f[dir::d0MM];
        const real f0PM = f[dir::d0PM];
        const real f0MP = f[dir::d0MP];
        const real fPPP = f[dir::dPPP];
        const real fMPP = f[dir::dMPP];
        const real fPMP = f[dir::dPMP];
        const real fMMP = f[dir::dMMP];
        const real fPPM = f[dir::dPPM];
        const real fMPM = f[dir::dMPM];
        const real fPMM = f[dir::dPMM];
        const real fMMM = f[dir::dMMM];

        real oneOverRho;
        getCompressibleMacroscopicValues(f, this->drho, oneOverRho, this->velocityX, this->velocityY, this->velocityZ);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Calculate second order moments for interpolation
        //!
        // example: kxxMzz: moment, second derivative in x direction minus the second derivative in z direction

        this->kxyFromfcNEQ = -c3o1 * omega *
                    ((fMM0 + fMMM + fMMP - fMP0 - fMPM - fMPP - fPM0 - fPMM - fPMP + fPP0 + fPPM + fPPP) *
                        oneOverRho -
                        ((this->velocityX * this->velocityY)));
        this->kyzFromfcNEQ = -c3o1 * omega *
                    ((f0MM + fPMM + fMMM - f0MP - fPMP - fMMP - f0PM - fPPM - fMPM + f0PP + fPPP + fMPP) *
                        oneOverRho -
                        ((this->velocityY * this->velocityZ)));
        this->kxzFromfcNEQ = -c3o1 * omega *
                    ((fM0M + fMMM + fMPM - fM0P - fMMP - fMPP - fP0M - fPMM - fPPM + fP0P + fPMP + fPPP) *
                        oneOverRho -
                        ((this->velocityX * this->velocityZ)));
        this->kxxMyyFromfcNEQ = -c3o2 * omega *
                        ((fM0M + fM00 + fM0P - f0MM - f0M0 - f0MP - f0PM - f0P0 - f0PP + fP0M + fP00 + fP0P) *
                        oneOverRho -
                        ((this->velocityX * this->velocityX - this->velocityY * this->velocityY)));
        this->kxxMzzFromfcNEQ = -c3o2 * omega *
                        ((fMM0 + fM00 + fMP0 - f0MM - f0MP - f00M - f00P - f0PM - f0PP + fPM0 + fP00 + fPP0) *
                        oneOverRho -
                        ((this->velocityX * this->velocityX - this->velocityZ * this->velocityZ)));
    }

};

class MomentsOnSourceNodeSet
{
private:
    vf::lbm::MomentsOnSourceNode momentsPPP;
    vf::lbm::MomentsOnSourceNode momentsMPP;
    vf::lbm::MomentsOnSourceNode momentsPMP;
    vf::lbm::MomentsOnSourceNode momentsMMP;
    vf::lbm::MomentsOnSourceNode momentsPPM;
    vf::lbm::MomentsOnSourceNode momentsMPM;
    vf::lbm::MomentsOnSourceNode momentsPMM;
    vf::lbm::MomentsOnSourceNode momentsMMM;

public:
    __host__ __device__ void calculatePPP(const real *const f, const real omega)
    {
        momentsPPP.calculate(f, omega);
    }

    __host__ __device__ void calculateMPP(const real *const f, const real omega)
    {
        momentsMPP.calculate(f, omega);
    }

    __host__ __device__ void calculatePMP(const real *const f, const real omega)
    {
        momentsPMP.calculate(f, omega);
    }

    __host__ __device__ void calculateMMP(const real *const f, const real omega)
    {
        momentsMMP.calculate(f, omega);
    }

    __host__ __device__ void calculatePPM(const real *const f, const real omega)
    {
        momentsPPM.calculate(f, omega);
    }

    __host__ __device__ void calculateMPM(const real *const f, const real omega)
    {
        momentsMPM.calculate(f, omega);
    }

    __host__ __device__ void calculatePMM(const real *const f, const real omega)
    {
        momentsPMM.calculate(f, omega);
    }

    __host__ __device__ void calculateMMM(const real *const f, const real omega)
    {
        momentsMMM.calculate(f, omega);
    }

    __host__ __device__ void calculateCoefficients(InterpolationCoefficients &coefficients, real xoff, real yoff,
                                                   real zoff) const
    {
        real& a000 = coefficients.a000;
        real& b000 = coefficients.b000;
        real& c000 = coefficients.c000;
        real& d000 = coefficients.d000;

        real& a100 = coefficients.a100;
        real& b100 = coefficients.b100;
        real& c100 = coefficients.c100;
        real& d100 = coefficients.d100;

        real& a010 = coefficients.a010;
        real& b010 = coefficients.b010;
        real& c010 = coefficients.c010;
        real& d010 = coefficients.d010;

        real& a001 = coefficients.a001;
        real& b001 = coefficients.b001;
        real& c001 = coefficients.c001;
        real& d001 = coefficients.d001;

        real& d110 = coefficients.d110, &d101 = coefficients.d101, &d011 = coefficients.d011;
        
        real& a200 = coefficients.a200, &a020 = coefficients.a020, &a002 = coefficients.a002;
        real& b200 = coefficients.b200, &b020 = coefficients.b020, &b002 = coefficients.b002;
        real& c200 = coefficients.c200, &c020 = coefficients.c020, &c002 = coefficients.c002;

        real& a110 = coefficients.a110, &a101 = coefficients.a101, &a011 = coefficients.a011;
        real& b110 = coefficients.b110, &b101 = coefficients.b101, &b011 = coefficients.b011;
        real& c110 = coefficients.c110, &c101 = coefficients.c101, &c011 = coefficients.c011;

        real &a111 = coefficients.a111, &b111 = coefficients.b111, &c111 = coefficients.c111, &d111 = coefficients.d111;

        real &LaplaceRho = coefficients.LaplaceRho;

        const real xoffsq = xoff * xoff;
        const real yoffsq = yoff * yoff;
        const real zoffsq = zoff * zoff;

        const real drhoPPP = momentsPPP.drho, vx1PPP = momentsPPP.velocityX, vx2PPP = momentsPPP.velocityY, vx3PPP = momentsPPP.velocityZ;
        const real drhoMPP = momentsMPP.drho, vx1MPP = momentsMPP.velocityX, vx2MPP = momentsMPP.velocityY, vx3MPP = momentsMPP.velocityZ;
        const real drhoPMP = momentsPMP.drho, vx1PMP = momentsPMP.velocityX, vx2PMP = momentsPMP.velocityY, vx3PMP = momentsPMP.velocityZ;
        const real drhoMMP = momentsMMP.drho, vx1MMP = momentsMMP.velocityX, vx2MMP = momentsMMP.velocityY, vx3MMP = momentsMMP.velocityZ;
        const real drhoPPM = momentsPPM.drho, vx1PPM = momentsPPM.velocityX, vx2PPM = momentsPPM.velocityY, vx3PPM = momentsPPM.velocityZ;
        const real drhoMPM = momentsMPM.drho, vx1MPM = momentsMPM.velocityX, vx2MPM = momentsMPM.velocityY, vx3MPM = momentsMPM.velocityZ;
        const real drhoPMM = momentsPMM.drho, vx1PMM = momentsPMM.velocityX, vx2PMM = momentsPMM.velocityY, vx3PMM = momentsPMM.velocityZ;
        const real drhoMMM = momentsMMM.drho, vx1MMM = momentsMMM.velocityX, vx2MMM = momentsMMM.velocityY, vx3MMM = momentsMMM.velocityZ;

        // second order moments at the source nodes
        const real kxyFromfcNEQPPP = momentsPPP.kxyFromfcNEQ, kyzFromfcNEQPPP = momentsPPP.kyzFromfcNEQ, kxzFromfcNEQPPP = momentsPPP.kxzFromfcNEQ, kxxMyyFromfcNEQPPP = momentsPPP.kxxMyyFromfcNEQ, kxxMzzFromfcNEQPPP = momentsPPP.kxxMzzFromfcNEQ;
        const real kxyFromfcNEQMPP = momentsMPP.kxyFromfcNEQ, kyzFromfcNEQMPP = momentsMPP.kyzFromfcNEQ, kxzFromfcNEQMPP = momentsMPP.kxzFromfcNEQ, kxxMyyFromfcNEQMPP = momentsMPP.kxxMyyFromfcNEQ, kxxMzzFromfcNEQMPP = momentsMPP.kxxMzzFromfcNEQ;
        const real kxyFromfcNEQPMP = momentsPMP.kxyFromfcNEQ, kyzFromfcNEQPMP = momentsPMP.kyzFromfcNEQ, kxzFromfcNEQPMP = momentsPMP.kxzFromfcNEQ, kxxMyyFromfcNEQPMP = momentsPMP.kxxMyyFromfcNEQ, kxxMzzFromfcNEQPMP = momentsPMP.kxxMzzFromfcNEQ;
        const real kxyFromfcNEQMMP = momentsMMP.kxyFromfcNEQ, kyzFromfcNEQMMP = momentsMMP.kyzFromfcNEQ, kxzFromfcNEQMMP = momentsMMP.kxzFromfcNEQ, kxxMyyFromfcNEQMMP = momentsMMP.kxxMyyFromfcNEQ, kxxMzzFromfcNEQMMP = momentsMMP.kxxMzzFromfcNEQ;
        const real kxyFromfcNEQPPM = momentsPPM.kxyFromfcNEQ, kyzFromfcNEQPPM = momentsPPM.kyzFromfcNEQ, kxzFromfcNEQPPM = momentsPPM.kxzFromfcNEQ, kxxMyyFromfcNEQPPM = momentsPPM.kxxMyyFromfcNEQ, kxxMzzFromfcNEQPPM = momentsPPM.kxxMzzFromfcNEQ;
        const real kxyFromfcNEQMPM = momentsMPM.kxyFromfcNEQ, kyzFromfcNEQMPM = momentsMPM.kyzFromfcNEQ, kxzFromfcNEQMPM = momentsMPM.kxzFromfcNEQ, kxxMyyFromfcNEQMPM = momentsMPM.kxxMyyFromfcNEQ, kxxMzzFromfcNEQMPM = momentsMPM.kxxMzzFromfcNEQ;
        const real kxyFromfcNEQPMM = momentsPMM.kxyFromfcNEQ, kyzFromfcNEQPMM = momentsPMM.kyzFromfcNEQ, kxzFromfcNEQPMM = momentsPMM.kxzFromfcNEQ, kxxMyyFromfcNEQPMM = momentsPMM.kxxMyyFromfcNEQ, kxxMzzFromfcNEQPMM = momentsPMM.kxxMzzFromfcNEQ;
        const real kxyFromfcNEQMMM = momentsMMM.kxyFromfcNEQ, kyzFromfcNEQMMM = momentsMMM.kyzFromfcNEQ, kxzFromfcNEQMMM = momentsMMM.kxzFromfcNEQ, kxxMyyFromfcNEQMMM = momentsMMM.kxxMyyFromfcNEQ, kxxMzzFromfcNEQMMM = momentsMMM.kxxMzzFromfcNEQ;

        a000 = c1o64 * (
                c2o1 * (
                ((kxyFromfcNEQMMM - kxyFromfcNEQPPP) + (kxyFromfcNEQMMP - kxyFromfcNEQPPM)) + ((kxyFromfcNEQPMM - kxyFromfcNEQMPP) + (kxyFromfcNEQPMP - kxyFromfcNEQMPM)) + 
                ((kxzFromfcNEQMMM - kxzFromfcNEQPPP) + (kxzFromfcNEQPPM - kxzFromfcNEQMMP)) + ((kxzFromfcNEQPMM - kxzFromfcNEQMPP) + (kxzFromfcNEQMPM - kxzFromfcNEQPMP)) + 
                ((vx2PPP + vx2MMM) + (vx2PPM + vx2MMP)) - ((vx2MPP + vx2PMM) + (vx2MPM + vx2PMP)) + 
                ((vx3PPP + vx3MMM) - (vx3PPM + vx3MMP)) + ((vx3PMP + vx3MPM) - (vx3MPP + vx3PMM))) + 
                c8o1 * (((vx1PPP + vx1MMM) + (vx1PPM + vx1MMP)) + ((vx1MPP + vx1PMM) + (vx1PMP + vx1MPM))) +
                ((kxxMyyFromfcNEQMMM - kxxMyyFromfcNEQPPP) + (kxxMyyFromfcNEQMMP - kxxMyyFromfcNEQPPM)) + 
                ((kxxMyyFromfcNEQMPP - kxxMyyFromfcNEQPMM) + (kxxMyyFromfcNEQMPM - kxxMyyFromfcNEQPMP)) +
                ((kxxMzzFromfcNEQMMM - kxxMzzFromfcNEQPPP) + (kxxMzzFromfcNEQMMP - kxxMzzFromfcNEQPPM)) + 
                ((kxxMzzFromfcNEQMPP - kxxMzzFromfcNEQPMM) + (kxxMzzFromfcNEQMPM - kxxMzzFromfcNEQPMP)));
        b000 = c1o64 * (
                c2o1 * (
                ((kxxMyyFromfcNEQPPP - kxxMyyFromfcNEQMMM) + (kxxMyyFromfcNEQPPM - kxxMyyFromfcNEQMMP)) + 
                ((kxxMyyFromfcNEQMPP - kxxMyyFromfcNEQPMM) + (kxxMyyFromfcNEQMPM - kxxMyyFromfcNEQPMP)) + 
                ((kxyFromfcNEQMMM - kxyFromfcNEQPPP) + (kxyFromfcNEQMMP - kxyFromfcNEQPPM)) + 
                ((kxyFromfcNEQMPP - kxyFromfcNEQPMM) + (kxyFromfcNEQMPM - kxyFromfcNEQPMP)) + 
                ((kyzFromfcNEQMMM - kyzFromfcNEQPPP) + (kyzFromfcNEQPPM - kyzFromfcNEQMMP)) + 
                ((kyzFromfcNEQPMM - kyzFromfcNEQMPP) + (kyzFromfcNEQMPM - kyzFromfcNEQPMP)) + 
                ((vx1PPP + vx1MMM) + (vx1PPM + vx1MMP)) - ((vx1MPM + vx1MPP) + (vx1PMM + vx1PMP)) + 
                ((vx3PPP + vx3MMM) - (vx3PPM + vx3MMP)) + ((vx3MPP + vx3PMM) - (vx3MPM + vx3PMP))) + 
                c8o1 * (((vx2PPP + vx2MMM) + (vx2PPM + vx2MMP)) + ((vx2MPP + vx2PMM) + (vx2MPM + vx2PMP))) + 
                ((kxxMzzFromfcNEQMMM - kxxMzzFromfcNEQPPP) + (kxxMzzFromfcNEQMMP - kxxMzzFromfcNEQPPM)) +
                ((kxxMzzFromfcNEQPMM - kxxMzzFromfcNEQMPP) + (kxxMzzFromfcNEQPMP - kxxMzzFromfcNEQMPM)));
        c000 = c1o64 * ( 
                c2o1 * (
                ((kxxMzzFromfcNEQPPP - kxxMzzFromfcNEQMMM) + (kxxMzzFromfcNEQMMP - kxxMzzFromfcNEQPPM)) + 
                ((kxxMzzFromfcNEQMPP - kxxMzzFromfcNEQPMM) + (kxxMzzFromfcNEQPMP - kxxMzzFromfcNEQMPM)) + 
                ((kxzFromfcNEQMMM - kxzFromfcNEQPPP) + (kxzFromfcNEQMMP - kxzFromfcNEQPPM)) + 
                ((kxzFromfcNEQMPP - kxzFromfcNEQPMM) + (kxzFromfcNEQMPM - kxzFromfcNEQPMP)) + 
                ((kyzFromfcNEQMMM - kyzFromfcNEQPPP) + (kyzFromfcNEQMMP - kyzFromfcNEQPPM)) + 
                ((kyzFromfcNEQPMM - kyzFromfcNEQMPP) + (kyzFromfcNEQPMP - kyzFromfcNEQMPM)) + 
                ((vx1PPP + vx1MMM) - (vx1MMP + vx1PPM)) + ((vx1MPM + vx1PMP) - (vx1MPP + vx1PMM)) + 
                ((vx2PPP + vx2MMM) - (vx2MMP + vx2PPM)) + ((vx2MPP + vx2PMM) - (vx2MPM + vx2PMP))) + 
                c8o1 * (((vx3PPP + vx3MMM) + (vx3PPM + vx3MMP)) + ((vx3PMM + vx3MPP) + (vx3PMP + vx3MPM))) +
                ((kxxMyyFromfcNEQMMM - kxxMyyFromfcNEQPPP) + (kxxMyyFromfcNEQPPM - kxxMyyFromfcNEQMMP)) + 
                ((kxxMyyFromfcNEQPMM - kxxMyyFromfcNEQMPP) + (kxxMyyFromfcNEQMPM - kxxMyyFromfcNEQPMP)));

        a100 = c1o4 * (((vx1PPP - vx1MMM) + (vx1PPM - vx1MMP)) + ((vx1PMM - vx1MPP) + (vx1PMP - vx1MPM)));
        b100 = c1o4 * (((vx2PPP - vx2MMM) + (vx2PPM - vx2MMP)) + ((vx2PMM - vx2MPP) + (vx2PMP - vx2MPM)));
        c100 = c1o4 * (((vx3PPP - vx3MMM) + (vx3PPM - vx3MMP)) + ((vx3PMM - vx3MPP) + (vx3PMP - vx3MPM)));

        a200 = c1o16 * ( 
                c2o1 * (
                ((vx2PPP + vx2MMM) + (vx2PPM - vx2MPP)) + ((vx2MMP - vx2PMM) - (vx2MPM + vx2PMP)) + 
                ((vx3PPP + vx3MMM) - (vx3PPM + vx3MPP)) + ((vx3MPM + vx3PMP) - (vx3MMP + vx3PMM))) + 
                ((kxxMyyFromfcNEQPPP - kxxMyyFromfcNEQMMM) + (kxxMyyFromfcNEQPPM - kxxMyyFromfcNEQMMP)) + 
                ((kxxMyyFromfcNEQPMM - kxxMyyFromfcNEQMPP) + (kxxMyyFromfcNEQPMP - kxxMyyFromfcNEQMPM)) + 
                ((kxxMzzFromfcNEQPPP - kxxMzzFromfcNEQMMM) + (kxxMzzFromfcNEQPPM - kxxMzzFromfcNEQMMP)) + 
                ((kxxMzzFromfcNEQPMM - kxxMzzFromfcNEQMPP) + (kxxMzzFromfcNEQPMP - kxxMzzFromfcNEQMPM)));
        b200 = c1o8 * (
                c2o1 * (
                -((vx1PPP + vx1MMM) + (vx1PPM + vx1MMP)) + ((vx1MPP + vx1PMM) + (vx1MPM + vx1PMP))) +
                ((kxyFromfcNEQPPP - kxyFromfcNEQMMM) + (kxyFromfcNEQPPM - kxyFromfcNEQMMP)) + 
                ((kxyFromfcNEQPMM - kxyFromfcNEQMPP) + (kxyFromfcNEQPMP - kxyFromfcNEQMPM)));
        c200 = c1o8 * (
                c2o1 * (
                ((vx1PPM + vx1MMP) - (vx1PPP + vx1MMM)) + ((vx1MPP + vx1PMM) - (vx1MPM + vx1PMP))) +
                ((kxzFromfcNEQPPP - kxzFromfcNEQMMM) + (kxzFromfcNEQPPM - kxzFromfcNEQMMP)) + 
                ((kxzFromfcNEQPMM - kxzFromfcNEQMPP) + (kxzFromfcNEQPMP - kxzFromfcNEQMPM)));

        a010 = c1o4 * (((vx1PPP - vx1MMM) + (vx1PPM - vx1MMP)) + ((vx1MPP - vx1PMM) + (vx1MPM - vx1PMP)));
        b010 = c1o4 * (((vx2PPP - vx2MMM) + (vx2PPM - vx2MMP)) + ((vx2MPP - vx2PMM) + (vx2MPM - vx2PMP)));
        c010 = c1o4 * (((vx3PPP - vx3MMM) + (vx3PPM - vx3MMP)) + ((vx3MPP - vx3PMM) + (vx3MPM - vx3PMP)));

        a020 = c1o8 * (
                c2o1 * (-((vx2PPP + vx2MMM) + (vx2MMP + vx2PPM)) + ((vx2MPP + vx2PMM) + (vx2MPM + vx2PMP))) +
                ((kxyFromfcNEQPPP - kxyFromfcNEQMMM) + (kxyFromfcNEQPPM - kxyFromfcNEQMMP)) + 
                ((kxyFromfcNEQMPP - kxyFromfcNEQPMM) + (kxyFromfcNEQMPM - kxyFromfcNEQPMP)));
        b020 = c1o16 * (
                c2o1 * (
                ((kxxMyyFromfcNEQMMM - kxxMyyFromfcNEQPPP) + (kxxMyyFromfcNEQMMP - kxxMyyFromfcNEQPPM)) +
                ((kxxMyyFromfcNEQPMM - kxxMyyFromfcNEQMPP) + (kxxMyyFromfcNEQPMP - kxxMyyFromfcNEQMPM)) +
                ((vx1PPP + vx1MMM) + (vx1PPM + vx1MMP)) - ((vx1MPP + vx1PMM) + (vx1PMP + vx1MPM)) + 
                ((vx3PPP + vx3MMM) - (vx3PPM + vx3MMP)) + ((vx3MPP + vx3PMM) - (vx3MPM + vx3PMP))) +
                ((kxxMzzFromfcNEQPPP - kxxMzzFromfcNEQMMM) + (kxxMzzFromfcNEQPPM - kxxMzzFromfcNEQMMP)) + 
                ((kxxMzzFromfcNEQMPP - kxxMzzFromfcNEQPMM) + (kxxMzzFromfcNEQMPM - kxxMzzFromfcNEQPMP)));
        c020 = c1o8 * (
                c2o1 * (((vx2MMP + vx2PPM) - (vx2PPP + vx2MMM)) + ((vx2PMP + vx2MPM) - (vx2MPP + vx2PMM))) +
                ((kyzFromfcNEQPPP - kyzFromfcNEQMMM) + (kyzFromfcNEQPPM - kyzFromfcNEQMMP)) +
                ((kyzFromfcNEQMPP - kyzFromfcNEQPMM) + (kyzFromfcNEQMPM - kyzFromfcNEQPMP)));

        a001 = c1o4 * (((vx1PPP - vx1MMM) + (vx1MMP - vx1PPM)) + ((vx1MPP - vx1PMM) + (vx1PMP - vx1MPM)));
        b001 = c1o4 * (((vx2PPP - vx2MMM) + (vx2MMP - vx2PPM)) + ((vx2MPP - vx2PMM) + (vx2PMP - vx2MPM)));
        c001 = c1o4 * (((vx3PPP - vx3MMM) + (vx3MMP - vx3PPM)) + ((vx3MPP - vx3PMM) + (vx3PMP - vx3MPM)));

        a002 = c1o8 * (
                c2o1 * (((vx3PPM + vx3MMP) - (vx3PPP + vx3MMM)) + ((vx3MPP + vx3PMM) - (vx3PMP + vx3MPM))) +
                        ((kxzFromfcNEQPPP - kxzFromfcNEQMMM) + (kxzFromfcNEQMMP - kxzFromfcNEQPPM)) +
                        ((kxzFromfcNEQPMP - kxzFromfcNEQMPM) + (kxzFromfcNEQMPP - kxzFromfcNEQPMM)));
        b002 = c1o8 * (
                c2o1 * (((vx3PPM + vx3MMP) - (vx3PPP + vx3MMM)) + ((vx3MPM + vx3PMP) - (vx3PMM + vx3MPP))) + 
                        ((kyzFromfcNEQPPP - kyzFromfcNEQMMM) + (kyzFromfcNEQMMP - kyzFromfcNEQPPM)) + 
                        ((kyzFromfcNEQPMP - kyzFromfcNEQMPM) + (kyzFromfcNEQMPP - kyzFromfcNEQPMM)));
        c002 = c1o16 * (
                c2o1 * (
                ((kxxMzzFromfcNEQMMM - kxxMzzFromfcNEQPPP) + (kxxMzzFromfcNEQPPM - kxxMzzFromfcNEQMMP)) + 
                ((kxxMzzFromfcNEQMPM - kxxMzzFromfcNEQPMP) + (kxxMzzFromfcNEQPMM - kxxMzzFromfcNEQMPP)) + 
                ((vx1PPP + vx1MMM) - (vx1MMP + vx1PPM)) + ((vx1MPM + vx1PMP) - (vx1PMM + vx1MPP)) + 
                ((vx2PPP + vx2MMM) - (vx2MMP + vx2PPM)) + ((vx2PMM + vx2MPP) - (vx2MPM + vx2PMP))) + 
                ((kxxMyyFromfcNEQPPP - kxxMyyFromfcNEQMMM) + (kxxMyyFromfcNEQMMP - kxxMyyFromfcNEQPPM)) +
                ((kxxMyyFromfcNEQPMP - kxxMyyFromfcNEQMPM) + (kxxMyyFromfcNEQMPP - kxxMyyFromfcNEQPMM)));

        a110 = c1o2 * (((vx1PPP + vx1MMM) + (vx1MMP + vx1PPM)) - ((vx1MPM + vx1PMP) + (vx1PMM + vx1MPP)));
        b110 = c1o2 * (((vx2PPP + vx2MMM) + (vx2MMP + vx2PPM)) - ((vx2MPM + vx2PMP) + (vx2PMM + vx2MPP)));
        c110 = c1o2 * (((vx3PPP + vx3MMM) + (vx3MMP + vx3PPM)) - ((vx3MPM + vx3PMP) + (vx3PMM + vx3MPP)));

        a101 = c1o2 * (((vx1PPP + vx1MMM) - (vx1MMP + vx1PPM)) + ((vx1MPM + vx1PMP) - (vx1PMM + vx1MPP)));
        b101 = c1o2 * (((vx2PPP + vx2MMM) - (vx2MMP + vx2PPM)) + ((vx2MPM + vx2PMP) - (vx2PMM + vx2MPP)));
        c101 = c1o2 * (((vx3PPP + vx3MMM) - (vx3MMP + vx3PPM)) + ((vx3MPM + vx3PMP) - (vx3PMM + vx3MPP)));
        
        a011 = c1o2 * (((vx1PPP + vx1MMM) - (vx1MMP + vx1PPM)) + ((vx1PMM + vx1MPP) - (vx1MPM + vx1PMP)));
        b011 = c1o2 * (((vx2PPP + vx2MMM) - (vx2MMP + vx2PPM)) + ((vx2PMM + vx2MPP) - (vx2MPM + vx2PMP)));
        c011 = c1o2 * (((vx3PPP + vx3MMM) - (vx3MMP + vx3PPM)) + ((vx3PMM + vx3MPP) - (vx3MPM + vx3PMP)));

        a111 = ((vx1PPP - vx1MMM) + (vx1MMP - vx1PPM)) + ((vx1MPM - vx1PMP) + (vx1PMM - vx1MPP));
        b111 = ((vx2PPP - vx2MMM) + (vx2MMP - vx2PPM)) + ((vx2MPM - vx2PMP) + (vx2PMM - vx2MPP));
        c111 = ((vx3PPP - vx3MMM) + (vx3MMP - vx3PPM)) + ((vx3MPM - vx3PMP) + (vx3PMM - vx3MPP));

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //!- Calculate coefficients for the polynomial interpolation of the pressure
        //! 
        LaplaceRho = 
            ((xoff != c0o1) || (yoff != c0o1) || (zoff != c0o1))
            ? c0o1 : -c3o1 * (a100 * a100 + b010 * b010 + c001 * c001) - c6o1 * (b100 * a010 + c100 * a001 + c010 * b001);
        d000 = c1o8 * (((drhoPPP + drhoMMM) + (drhoPPM + drhoMMP)) + ((drhoPMM + drhoMPP) + (drhoPMP + drhoMPM)));
        d100 = c1o4 * (((drhoPPP - drhoMMM) + (drhoPPM - drhoMMP)) + ((drhoPMM - drhoMPP) + (drhoPMP - drhoMPM)));
        d010 = c1o4 * (((drhoPPP - drhoMMM) + (drhoPPM - drhoMMP)) + ((drhoMPP - drhoPMM) + (drhoMPM - drhoPMP)));
        d001 = c1o4 * (((drhoPPP - drhoMMM) + (drhoMMP - drhoPPM)) + ((drhoMPP - drhoPMM) + (drhoPMP - drhoMPM)));
        d110 = c1o2 * (((drhoPPP + drhoMMM) + (drhoPPM + drhoMMP)) - ((drhoPMM + drhoMPP) + (drhoPMP + drhoMPM)));
        d101 = c1o2 * (((drhoPPP + drhoMMM) - (drhoPPM + drhoMMP)) + ((drhoPMP + drhoMPM) - (drhoPMM + drhoMPP)));
        d011 = c1o2 * (((drhoPPP + drhoMMM) - (drhoPPM + drhoMMP)) + ((drhoPMM + drhoMPP) - (drhoPMP + drhoMPM)));

        d111 = (((drhoPPP - drhoMMM) + (drhoMMP - drhoPPM)) + ((drhoPMM - drhoMPP) + (drhoMPM - drhoPMP)));

        //////////////////////////////////////////////////////////////////////////
        //! - Extrapolation for refinement in to the wall (polynomial coefficients)
        //!
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // x------x
        // |      |
        // |   ---+--->X
        // |      |  |
        // x------x  |
        //          offset-vector
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        a000 = a000 + xoff * a100 + yoff * a010 + zoff * a001 + xoffsq * a200 + yoffsq * a020 + zoffsq * a002 +
            xoff * yoff * a110 + xoff * zoff * a101 + yoff * zoff * a011;
        a100 = a100 + c2o1 * xoff * a200 + yoff * a110 + zoff * a101;
        a010 = a010 + c2o1 * yoff * a020 + xoff * a110 + zoff * a011;
        a001 = a001 + c2o1 * zoff * a002 + xoff * a101 + yoff * a011;
        b000 = b000 + xoff * b100 + yoff * b010 + zoff * b001 + xoffsq * b200 + yoffsq * b020 + zoffsq * b002 +
                xoff * yoff * b110 + xoff * zoff * b101 + yoff * zoff * b011;
        b100 = b100 + c2o1 * xoff * b200 + yoff * b110 + zoff * b101;
        b010 = b010 + c2o1 * yoff * b020 + xoff * b110 + zoff * b011;
        b001 = b001 + c2o1 * zoff * b002 + xoff * b101 + yoff * b011;
        c000 = c000 + xoff * c100 + yoff * c010 + zoff * c001 + xoffsq * c200 + yoffsq * c020 + zoffsq * c002 +
                xoff * yoff * c110 + xoff * zoff * c101 + yoff * zoff * c011;
        c100 = c100 + c2o1 * xoff * c200 + yoff * c110 + zoff * c101;
        c010 = c010 + c2o1 * yoff * c020 + xoff * c110 + zoff * c011;
        c001 = c001 + c2o1 * zoff * c002 + xoff * c101 + yoff * c011;
        d000 = d000 + xoff * d100 + yoff * d010 + zoff * d001 + 
                xoff * yoff * d110 + xoff * zoff * d101 + yoff * zoff * d011;

        d100 = d100 + yoff * d110 + zoff * d101;
        d010 = d010 + xoff * d110 + zoff * d011;
        d001 = d001 + xoff * d101 + yoff * d011;
    }

};

} // namespace vf::lbm

#endif

//! \}

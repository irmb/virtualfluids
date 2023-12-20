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
//! \addtogroup cpu_BoundaryConditions BoundaryConditions
//! \ingroup cpu_core core
//! \{
//! \author Sören Freudiger
//=======================================================================================
#ifndef BoundaryConditions_H
#define BoundaryConditions_H

#include <string>
#include <vector>

#include "D3Q27System.h"
#include "UbException.h"
#include "UbSystem.h"
#include "UbTuple.h"
#include "Vector3D.h"
#include <PointerDefinitions.h>
#include "basics/constants/NumericConstants.h"

//! Difenition of baundary conditions in grid generation
class BoundaryConditions
{
public:
    BoundaryConditions()
    {
        UB_STATIC_ASSERT(sizeof(long long) >= 8);
        UB_STATIC_ASSERT((sizeof(long long) * 8) >= (D3Q27System::FENDDIR + 1) * BoundaryConditions::optionDigits);

        for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++)
            q[fdir] = -999.;
    }
    virtual ~BoundaryConditions() = default;

    virtual bool isEmpty()
    {
        return (noslipBoundaryFlags & slipBoundaryFlags & velocityBoundaryFlags & densityBoundaryFlags) == 0;
    }
    virtual bool hasBoundaryCondition()
    {
        return (hasNoSlipBoundary() || hasSlipBoundary() || hasDensityBoundary() || hasVelocityBoundary() ||
                hasWallModelBoundary());
    }

    virtual bool hasBoundaryConditionFlag(const int &direction)
    {
        assert(direction >= D3Q27System::FSTARTDIR && direction <= D3Q27System::FENDDIR);

        return (hasNoSlipBoundaryFlag(direction) || hasSlipBoundaryFlag(direction) ||
                hasDensityBoundaryFlag(direction) || hasVelocityBoundaryFlag(direction) ||
                hasWallModelBoundaryFlag(direction));
    }

protected:
    void setFlagBits(long long &flag, const int &direction, const short &secOpt)
    {
        if ((secOpt + 1) > maxOptionVal)
            throw UbException(UB_EXARGS, "error: option > " + UbSystem::toString(maxOptionVal - 1));

        // all digits at the respective positions to "0"
        flag &= ~(maxOptionVal << (direction * optionDigits));
        // set all digits according to the flag at the respective positions
        flag |= ((long long)(secOpt + 1) << (direction * optionDigits));
    }

public:
    /*===================== NoSlip Boundary ==================================================*/
    void setNoSlipBoundaryFlag(const int &direction, const short &secOpt = 0)
    {
        this->setFlagBits(noslipBoundaryFlags, direction, secOpt);
    }
    void unsetNoSlipBoundaryFlag(const int &direction)
    {
        this->noslipBoundaryFlags &= ~(maxOptionVal << (direction * optionDigits));
    }
    void unsetNoSlipBoundary() { this->noslipBoundaryFlags = 0; }
    long long getNoSlipBoundary() { return this->noslipBoundaryFlags; }
    bool hasNoSlipBoundary() { return (noslipBoundaryFlags != 0); }
    bool hasNoSlipBoundaryFlag(const int &direction)
    {
        return (((noslipBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) != 0);
    }
    short getNoSlipSecondaryOption(const int &direction)
    {
        return (short)(((noslipBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) - 1);
    }
    /*===================== WallModel Boundary ==================================================*/
    void setWallModelBoundaryFlag(const int &direction, const short &secOpt = 0)
    {
        this->setFlagBits(wallModelBoundaryFlags, direction, secOpt);
    }
    void unsetWallModelBoundaryFlag(const int &direction)
    {
        this->wallModelBoundaryFlags &= ~(maxOptionVal << (direction * optionDigits));
    }
    void unsetWallModelBoundary() { this->wallModelBoundaryFlags = 0; }
    long long getWallModelBoundary() { return this->wallModelBoundaryFlags; }
    bool hasWallModelBoundary() { return (wallModelBoundaryFlags != 0); }
    bool hasWallModelBoundaryFlag(const int &direction)
    {
        return (((wallModelBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) != 0);
    }
    short getWallModelSecondaryOption(const int &direction)
    {
        return (short)(((wallModelBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) - 1);
    }
    /*===================== Slip-Solid Boundary ==================================================*/
    void setSlipBoundaryFlag(const int &direction, const short &secOpt = 0)
    {
        this->setFlagBits(slipBoundaryFlags, direction, secOpt);
    }
    void unsetSlipBoundaryFlag(const int &direction)
    {
        this->slipBoundaryFlags &= ~(maxOptionVal << (direction * optionDigits));
    }
    void unsetSlipBoundary() { this->slipBoundaryFlags = 0; }
    long long getSlipBoundary() { return this->slipBoundaryFlags; }
    bool hasSlipBoundary() { return (slipBoundaryFlags != 0); }
    bool hasSlipBoundaryFlag(const int &direction)
    {
        return (((slipBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) != 0);
    }
    short getSlipSecondaryOption(const int &direction)
    {
        return (short)(((slipBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) - 1);
    }
    void setNormalVector(const float &nx1, const float &nx2, const float &nx3)
    {
        this->nx1 = nx1;
        this->nx2 = nx2;
        this->nx3 = nx3;
    }
    UbTupleFloat3 getNormalVector() { return makeUbTuple(nx1, nx2, nx3); }

    /*============== Velocity Boundary ========================*/
    void setVelocityBoundaryFlag(const int &direction, const short &secOpt = 0)
    {
        this->setFlagBits(velocityBoundaryFlags, direction, secOpt);
    }
    void unsetVelocityBoundaryFlag(const int &direction)
    {
        this->velocityBoundaryFlags &= ~(maxOptionVal << (direction * optionDigits));
    }
    void unsetVelocityBoundary() { this->velocityBoundaryFlags = 0; }
    long long getVelocityBoundary() { return this->velocityBoundaryFlags; }
    bool hasVelocityBoundary() { return this->velocityBoundaryFlags != 0; }
    bool hasVelocityBoundaryFlag(const int &direction)
    {
        return (((velocityBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) != 0);
    }
    short getVelocitySecondaryOption(const int &direction)
    {
        return (short)(((velocityBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) - 1);
    }

    void setBoundaryVelocity(const Vector3D &vx)
    {
        setBoundaryVelocityX1((float)vx[0]);
        setBoundaryVelocityX2((float)vx[1]);
        setBoundaryVelocityX3((float)vx[2]);
    }
    void setBoundaryVelocityX1(const float &vx1) { this->bcVelocityX1 = vx1; }
    void setBoundaryVelocityX2(const float &vx2) { this->bcVelocityX2 = vx2; }
    void setBoundaryVelocityX3(const float &vx3) { this->bcVelocityX3 = vx3; }
    float getBoundaryVelocityX1() { return this->bcVelocityX1; }
    float getBoundaryVelocityX2() { return this->bcVelocityX2; }
    float getBoundaryVelocityX3() { return this->bcVelocityX3; }
    float getBoundaryVelocity(const int &direction)
    {
        using namespace vf::lbm::dir;

        switch (direction) {
            case dP00:
                return (float)(vf::basics::constant::c4o9 *
                               (+bcVelocityX1)); //(2/cs^2)(=6)*rho_0(=1 bei inkompr)*wi*u*ei mit cs=1/sqrt(3)
            case dM00:
                return (float)(vf::basics::constant::c4o9 *
                               (-bcVelocityX1)); // z.B. aus paper manfred MRT LB models in three dimensions (2002)
            case d0P0:
                return (float)(vf::basics::constant::c4o9 * (+bcVelocityX2));
            case d0M0:
                return (float)(vf::basics::constant::c4o9 * (-bcVelocityX2));
            case d00P:
                return (float)(vf::basics::constant::c4o9 * (+bcVelocityX3));
            case d00M:
                return (float)(vf::basics::constant::c4o9 * (-bcVelocityX3));
            case dPP0:
                return (float)(vf::basics::constant::c1o9 * (+bcVelocityX1 + bcVelocityX2));
            case dMM0:
                return (float)(vf::basics::constant::c1o9 * (-bcVelocityX1 - bcVelocityX2));
            case dPM0:
                return (float)(vf::basics::constant::c1o9 * (+bcVelocityX1 - bcVelocityX2));
            case dMP0:
                return (float)(vf::basics::constant::c1o9 * (-bcVelocityX1 + bcVelocityX2));
            case dP0P:
                return (float)(vf::basics::constant::c1o9 * (+bcVelocityX1 + bcVelocityX3));
            case dM0M:
                return (float)(vf::basics::constant::c1o9 * (-bcVelocityX1 - bcVelocityX3));
            case dP0M:
                return (float)(vf::basics::constant::c1o9 * (+bcVelocityX1 - bcVelocityX3));
            case dM0P:
                return (float)(vf::basics::constant::c1o9 * (-bcVelocityX1 + bcVelocityX3));
            case d0PP:
                return (float)(vf::basics::constant::c1o9 * (+bcVelocityX2 + bcVelocityX3));
            case d0MM:
                return (float)(vf::basics::constant::c1o9 * (-bcVelocityX2 - bcVelocityX3));
            case d0PM:
                return (float)(vf::basics::constant::c1o9 * (+bcVelocityX2 - bcVelocityX3));
            case d0MP:
                return (float)(vf::basics::constant::c1o9 * (-bcVelocityX2 + bcVelocityX3));
            case dPPP:
                return (float)(vf::basics::constant::c1o36 * (+bcVelocityX1 + bcVelocityX2 + bcVelocityX3));
            case dMMM:
                return (float)(vf::basics::constant::c1o36 * (-bcVelocityX1 - bcVelocityX2 - bcVelocityX3));
            case dPPM:
                return (float)(vf::basics::constant::c1o36 * (+bcVelocityX1 + bcVelocityX2 - bcVelocityX3));
            case dMMP:
                return (float)(vf::basics::constant::c1o36 * (-bcVelocityX1 - bcVelocityX2 + bcVelocityX3));
            case dPMP:
                return (float)(vf::basics::constant::c1o36 * (+bcVelocityX1 - bcVelocityX2 + bcVelocityX3));
            case dMPM:
                return (float)(vf::basics::constant::c1o36 * (-bcVelocityX1 + bcVelocityX2 - bcVelocityX3));
            case dPMM:
                return (float)(vf::basics::constant::c1o36 * (+bcVelocityX1 - bcVelocityX2 - bcVelocityX3));
            case dMPP:
                return (float)(vf::basics::constant::c1o36 * (-bcVelocityX1 + bcVelocityX2 + bcVelocityX3));
            default:
                throw UbException(UB_EXARGS, "unknown error");
        }
    }
    
    //Multiphase
    void setBoundaryPhaseField(const float &phiBC) { this->bcPhaseField = phiBC; }
    float getBoundaryPhaseField() { return this->bcPhaseField; }

    /*============== Density Boundary ========================*/
    void setDensityBoundaryFlag(const int &direction, const short &secOpt = 0)
    {
        this->setFlagBits(densityBoundaryFlags, direction, secOpt);
    }
    void unsetDensityBoundaryFlag(const int &direction)
    {
        this->densityBoundaryFlags &= ~(maxOptionVal << (direction * optionDigits));
    }
    void unsetDensityBoundary() { this->densityBoundaryFlags = 0; }
    long long getDensityBoundary() { return this->densityBoundaryFlags; }
    bool hasDensityBoundary() { return (this->densityBoundaryFlags != 0); }
    bool hasDensityBoundaryFlag(const int &direction)
    {
        return (((densityBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) != 0);
    }
    short getDensitySecondaryOption(const int &direction)
    {
        return (short)(((densityBoundaryFlags >> (optionDigits * direction)) & maxOptionVal) - 1);
    }

    void setBoundaryDensity(float density) { this->bcDensity = density; }
    float getBoundaryDensity() { return this->bcDensity; }

    /*======================= Qs =============================*/
    void setQ(const float &val, const int &direction) { q[direction] = val; }
    float getQ(const int &direction) { return q[direction]; }

    virtual std::vector<std::string> getBCNames()
    {
        std::vector<std::string> tmp;
        tmp.push_back("NoSlipBC");
        tmp.push_back("SlipBC");
        tmp.push_back("VelocityBC");
        tmp.push_back("DensityBC");
        return tmp;
    }
    virtual std::vector<long long> getBCFlags()
    {
        std::vector<long long> tmp;
        tmp.push_back(noslipBoundaryFlags);
        tmp.push_back(slipBoundaryFlags);
        tmp.push_back(velocityBoundaryFlags);
        tmp.push_back(densityBoundaryFlags);
        return tmp;
    }

    static bool hasFlagForDirection(const long long &flag, const int &direction)
    {
        return (((flag >> (optionDigits * direction)) & maxOptionVal) != 0);
    }

    void setBCStrategyKey(char key) { bcStrategyKey = key; }
    char getBCStrategyKey() { return bcStrategyKey; }

public:
    static const int optionDigits = 2;   //--> 2 bits for secondary Option --> maxOptionVal = 7
    static const long long maxOptionVal; // = ( 1<<optionDigits ) - 1; //2^3-1 -> 7

protected:
    float q[D3Q27System::FENDDIR + 1];

    long long noslipBoundaryFlags{ 0 };
    long long slipBoundaryFlags{ 0 };
    long long velocityBoundaryFlags{ 0 };
    long long densityBoundaryFlags{ 0 };
    long long wallModelBoundaryFlags{ 0 };

    float bcVelocityX1{ 0.0 };
    float bcVelocityX2{ 0.0 };
    float bcVelocityX3{ 0.0 };
    float bcDensity{ 0.0 };
    float bcPhaseField{ 0.0 };

    float nx1{ vf::basics::constant::c0o1 }, nx2{ vf::basics::constant::c0o1 }, nx3{ vf::basics::constant::c0o1 };

    char bcStrategyKey { -1 };

private:
    friend class MPIIORestartSimulationObserver;
    friend class MPIIOMigrationSimulationObserver;
    friend class MPIIOMigrationBESimulationObserver;
};

#endif

//! \}

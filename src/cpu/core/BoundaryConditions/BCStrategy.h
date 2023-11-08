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
//! \file BCStrategy.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <PointerDefinitions.h>

#include "D3Q27System.h"

class DistributionArray3D;
class BCArray3D;
class BoundaryConditions;
class Block3D;

//! \brief Abstract class of baundary conditions algorithm
//! \details  BCStrategy provides interface for implementation of diferent boundary conditions
class BCStrategy
{
public:
    static const char VelocityBCStrategy                           = 0;
    static const char EqDensityBCStrategy                          = 1;
    static const char NonEqDensityBCStrategy                       = 2;
    static const char NoSlipBCStrategy                             = 3;
    static const char SlipBCStrategy                               = 4;
    static const char HighViscosityNoSlipBCStrategy                = 5;
    static const char ThinWallNoSlipBCStrategy                     = 6;
    static const char VelocityWithDensityBCStrategy                = 7;
    static const char NonReflectingOutflowBCStrategy               = 8;
    static const char ThixotropyVelocityBCStrategy                 = 9;
    static const char ThixotropyDensityBCStrategy                  = 10;
    static const char ThixotropyNoSlipBCStrategy                   = 11;
    static const char ThixotropyNonReflectingOutflowBCStrategy     = 12;
    static const char ThixotropyVelocityWithDensityBCStrategy      = 13;
    static const char RheologyBinghamModelNoSlipBCStrategy         = 14;
    static const char RheologyHerschelBulkleyModelNoSlipBCStrategy = 15;
    static const char SimpleVelocityBCStrategy                     = 16;
    static const char SimpleSlipBCStrategy                         = 17;
    static const char RheologyPowellEyringModelNoSlipBCStrategy    = 18;
    static const char RheologyBinghamModelVelocityBCStrategy       = 19;
    static const char MultiphaseNoSlipBCStrategy                   = 20;
    static const char MultiphaseVelocityBCStrategy                 = 21;
    static const char NonReflectingInflowBCStrategy                = 22;
    static const char NonReflectingOutflowWithRelaxationBCStrategy = 23;
    static const char MultiphasePressureBCStrategy                 = 24;

public:
    BCStrategy() = default;
    virtual ~BCStrategy() = default;

    virtual void addDistributions(SPtr<DistributionArray3D> distributions)   = 0;
    virtual void addDistributionsH(SPtr<DistributionArray3D> distributionsH) {}
    virtual void addDistributionsH2(SPtr<DistributionArray3D> distributionsH2) {}
    void setBlock(SPtr<Block3D> block);
    void setNodeIndex(int x1, int x2, int x3);
    void setBcPointer(SPtr<BoundaryConditions> bcPtr);
    void setCompressible(bool c);
    void setCollFactor(real cf);

    void setCollFactorL(real cf);
    void setCollFactorG(real cf);
    void setCollFactorPh(real cf);
    void setDensityRatio(real dr);
    void setPhiBound(real phiL, real phiH);

    char getType();
    bool isPreCollision();
    virtual SPtr<BCStrategy> clone() = 0;
    SPtr<BCArray3D> getBcArray();
    void setBcArray(SPtr<BCArray3D> bcarray);
    virtual void applyBC() = 0;
    bool getThixotropy(){ return thixotropy; };

protected:
    bool compressible { false };
    char type;
    bool preCollision;
    bool thixotropy { false };

    SPtr<BoundaryConditions> bcPtr;
    SPtr<DistributionArray3D> distributions;
    SPtr<DistributionArray3D> distributionsH;
    SPtr<DistributionArray3D> distributionsH2;
    SPtr<BCArray3D> bcArray;
    SPtr<Block3D> block;

    real collFactor;
    real collFactorL, collFactorG, collFactorPh;
    real densityRatio;
    real phiL, phiH;
    int x1, x2, x3;

    real compressibleFactor;

    using CalcMacrosFct    = void (*)(const real *const &, real &, real &, real &, real &);
    using CalcFeqForDirFct = real (*)(const int &, const real &, const real &, const real &,
                                         const real &);
    using CalcFeqFct = void (*)(real *const &, const real &, const real &, const real &, const real &);

    CalcFeqForDirFct calcFeqsForDirFct;
    CalcMacrosFct calcMacrosFct;
    CalcFeqFct calcFeqFct;
};

#endif

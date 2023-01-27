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
//! \file BCAlgorithm.h
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
//! \details  BCAlgorithm provides interface for implementation of diferent boundary conditions
class BCAlgorithm
{
public:
    static const char VelocityBCAlgorithm                          = 0;
    static const char EqDensityBCAlgorithm                         = 1;
    static const char NonEqDensityBCAlgorithm                      = 2;
    static const char NoSlipBCAlgorithm                            = 3;
    static const char SlipBCAlgorithm                              = 4;
    static const char HighViscosityNoSlipBCAlgorithm               = 5;
    static const char ThinWallNoSlipBCAlgorithm                    = 6;
    static const char VelocityWithDensityBCAlgorithm               = 7;
    static const char NonReflectingOutflowBCAlgorithm              = 8;
    static const char ThixotropyVelocityBCAlgorithm             = 9;
    static const char ThixotropyDensityBCAlgorithm              = 10;
    static const char ThixotropyNoSlipBCAlgorithm               = 11;
    static const char ThixotropyNonReflectingOutflowBCAlgorithm = 12;
    static const char ThixotropyVelocityWithDensityBCAlgorithm  = 13;
    static const char RheologyBinghamModelNoSlipBCAlgorithm                = 14;
    static const char RheologyHerschelBulkleyModelNoSlipBCAlgorithm        = 15;
    static const char SimpleVelocityBCAlgorithm                    = 16;
    static const char SimpleSlipBCAlgorithm                        = 17;
    static const char RheologyPowellEyringModelNoSlipBCAlgorithm           = 18;
    static const char RheologyBinghamModelVelocityBCAlgorithm              = 19;
    static const char MultiphaseNoSlipBCAlgorithm                  = 20;
    static const char MultiphaseVelocityBCAlgorithm = 21;



public:
    BCAlgorithm() = default;
    virtual ~BCAlgorithm() = default;

    virtual void addDistributions(SPtr<DistributionArray3D> distributions)   = 0;
    virtual void addDistributionsH(SPtr<DistributionArray3D> distributionsH) {}
    virtual void addDistributionsH2(SPtr<DistributionArray3D> distributionsH2) {}
    void setBlock(SPtr<Block3D> block);
    void setNodeIndex(int x1, int x2, int x3);
    void setBcPointer(SPtr<BoundaryConditions> bcPtr);
    void setCompressible(bool c);
    void setCollFactor(LBMReal cf);

    void setCollFactorL(LBMReal cf);
    void setCollFactorG(LBMReal cf);
    void setCollFactorPh(LBMReal cf);
    void setDensityRatio(LBMReal dr);
    void setPhiBound(LBMReal phiL, LBMReal phiH);

    char getType();
    bool isPreCollision();
    virtual SPtr<BCAlgorithm> clone() = 0;
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

    LBMReal collFactor;
    LBMReal collFactorL, collFactorG, collFactorPh;
    LBMReal densityRatio;
    LBMReal phiL, phiH;
    int x1, x2, x3;

    LBMReal compressibleFactor;

    using CalcMacrosFct    = void (*)(const LBMReal *const &, LBMReal &, LBMReal &, LBMReal &, LBMReal &);
    using CalcFeqForDirFct = LBMReal (*)(const int &, const LBMReal &, const LBMReal &, const LBMReal &,
                                         const LBMReal &);
    using CalcFeqFct = void (*)(LBMReal *const &, const LBMReal &, const LBMReal &, const LBMReal &, const LBMReal &);

    CalcFeqForDirFct calcFeqsForDirFct;
    CalcMacrosFct calcMacrosFct;
    CalcFeqFct calcFeqFct;
};

#endif

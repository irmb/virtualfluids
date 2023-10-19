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

#ifndef MultiphaseBCStrategyType_H
#define MultiphaseBCStrategyType_H

struct MultiphaseBCStrategyType
{
    //static const char VelocityBCStrategy = 0;
    //static const char EqDensityBCStrategy = 1;
    //static const char NonEqDensityBCStrategy = 2;
    //static const char NoSlipBCStrategy = 3;
    //static const char SlipBCStrategy = 4;
    //static const char HighViscosityNoSlipBCStrategy = 5;
    //static const char ThinWallNoSlipBCStrategy = 6;
    //static const char VelocityWithDensityBCStrategy = 7;
    //static const char NonReflectingOutflowBCStrategy = 8;
    //static const char ThixotropyVelocityBCStrategy = 9;
    //static const char ThixotropyDensityBCStrategy = 10;
    //static const char ThixotropyNoSlipBCStrategy = 11;
    //static const char ThixotropyNonReflectingOutflowBCStrategy = 12;
    //static const char ThixotropyVelocityWithDensityBCStrategy = 13;
    //static const char RheologyBinghamModelNoSlipBCStrategy = 14;
    //static const char RheologyHerschelBulkleyModelNoSlipBCStrategy = 15;
    //static const char SimpleVelocityBCStrategy = 16;
    //static const char SimpleSlipBCStrategy = 17;
    //static const char RheologyPowellEyringModelNoSlipBCStrategy = 18;
    //static const char RheologyBinghamModelVelocityBCStrategy = 19;
    //static const char MultiphaseNoSlipBCStrategy = 20;
    //static const char MultiphaseVelocityBCStrategy = 21;
    //static const char NonReflectingInflowBCStrategy = 22;
    //static const char NonReflectingOutflowWithRelaxationBCStrategy = 23;
    //static const char MultiphasePressureBCStrategy = 24;


    static const char MultiphaseNoSlipBCStrategy = 0;
    static const char MultiphaseVelocityBCStrategy = 1;
    static const char MultiphasePressureBCStrategy = 2;
    static const char MultiphaseNonReflectingOutflowBCStrategy = 3;
    static const char MultiphaseSlipBCStrategy = 4;
};

#endif
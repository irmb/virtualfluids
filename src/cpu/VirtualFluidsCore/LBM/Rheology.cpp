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
//! \file Rheology.cpp
//! \ingroup LBM
//! \author Konstantin Kutscher, Martin Geier
//=======================================================================================
#include "Rheology.h"

SPtr<Rheology> Rheology::instance = SPtr<Rheology>();
LBMReal Rheology::tau0 = 0;
LBMReal Rheology::k = 0;
LBMReal Rheology::n = 1;
LBMReal Rheology::omegaMin = 0;
LBMReal Rheology::beta = 0;
LBMReal Rheology::c = 0;
LBMReal Rheology::mu0 = 0;

//////////////////////////////////////////////////////////////////////////
SPtr<Rheology> Rheology::getInstance()
{
   if (!instance)
      instance = SPtr<Rheology>(new Rheology());
   return instance;
}

void Rheology::setYieldStress(LBMReal yieldStress)
{
	tau0 = yieldStress;
}
LBMReal Rheology::getYieldStress() const
{
	return tau0;
}
void Rheology::setViscosityParameter(LBMReal kParameter)
{
	k = kParameter;
}
LBMReal Rheology::getViscosityParameter() const
{
	return k;
}
void Rheology::setPowerIndex(LBMReal index)
{
	n = index;
}
LBMReal Rheology::getPowerIndex() const
{
	return n;
}

void Rheology::setOmegaMin(LBMReal omega)
{
	omegaMin = omega;
}
LBMReal Rheology::getOmegaMin() const
{
	return omegaMin;
}

void Rheology::setBeta(LBMReal PowellEyringBeta)
{
	beta = PowellEyringBeta;
}

LBMReal Rheology::getBeta() const
{
	return beta;
}

void Rheology::setC(LBMReal PowellEyringC)
{
	c = PowellEyringC;
}

LBMReal Rheology::getC() const
{
	return c;
}

void Rheology::setMu0(LBMReal mu)
{
	mu0 = mu;
}

LBMReal Rheology::getMu0() const
{
	return mu0;
}

Rheology::Rheology()
{
}
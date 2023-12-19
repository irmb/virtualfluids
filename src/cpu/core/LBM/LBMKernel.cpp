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
//! \addtogroup cpu_LBM LBM
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#include "LBMKernel.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "DataSet3D.h"

LBMKernel::LBMKernel()

{
    this->setForcingX1(0.0);
    this->setForcingX2(0.0);
    this->setForcingX3(0.0);
    dataSet     = std::make_shared<DataSet3D>();
    this->nx[0] = 0;
    this->nx[1] = 0;
    this->nx[2] = 0;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setBCSet(SPtr<BCSet> bcp) { bcSet = bcp; }
//////////////////////////////////////////////////////////////////////////
SPtr<BCSet> LBMKernel::getBCSet() const { return bcSet; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setCollisionFactor(real collFactor) { this->collFactor = collFactor; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getCollisionFactor() const { return collFactor; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX1(real forcingX1)
{
    this->muForcingX1.SetExpr(UbSystem::toString(forcingX1, realLim::digits10));
    this->checkFunction(muForcingX1);
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX2(real forcingX2)
{
    this->muForcingX2.SetExpr(UbSystem::toString(forcingX2, realLim::digits10));
    this->checkFunction(muForcingX2);
}
void LBMKernel::setForcingX3(real forcingX3)
{
    this->muForcingX3.SetExpr(UbSystem::toString(forcingX3, realLim::digits10));
    this->checkFunction(muForcingX3);
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX1(const mu::Parser &parser)
{
    this->checkFunction(parser);
    this->muForcingX1 = parser;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX2(const mu::Parser &parser)
{
    this->checkFunction(parser);
    this->muForcingX2 = parser;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX3(const mu::Parser &parser)
{
    this->checkFunction(parser);
    this->muForcingX3 = parser;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX1(const std::string &muParserString)
{
    this->muForcingX1.SetExpr(muParserString);
    this->checkFunction(muForcingX1);
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX2(const std::string &muParserString)
{
    this->muForcingX2.SetExpr(muParserString);
    this->checkFunction(muForcingX2);
}
void LBMKernel::setForcingX3(const std::string &muParserString)
{
    this->muForcingX3.SetExpr(muParserString);
    this->checkFunction(muForcingX3);
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::checkFunction(mu::Parser fct)
{
    using namespace vf::basics::constant;

    double x1 = c1o1, x2 = c1o1, x3 = c1o1, dt = c1o1, nue = c1o1, rho = c1o1;
    fct.DefineVar("x1", &x1);
    fct.DefineVar("x2", &x2);
    fct.DefineVar("x3", &x3);
    fct.DefineVar("dt", &dt);
    fct.DefineVar("nue", &nue);
    fct.DefineVar("rho", &rho);

    try {
        fct.Eval();
        fct.ClearVar();
    } catch (mu::ParserError &e) {
        throw UbException(UB_EXARGS, "function: " + e.GetExpr() + (std::string) "error: " + e.GetMsg() +
                                         (std::string) ", only x1,x2,x3,dx are allowed as variables");
    }
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setGhostLayerWidth(int witdh) { ghostLayerWidth = witdh; }
//////////////////////////////////////////////////////////////////////////
int LBMKernel::getGhostLayerWidth() const { return ghostLayerWidth; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setIndex(int x1, int x2, int x3)
{
    this->ix1 = x1;
    this->ix2 = x2;
    this->ix3 = x3;
}
//////////////////////////////////////////////////////////////////////////
SPtr<DataSet3D> LBMKernel::getDataSet() const { return this->dataSet; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getDeltaT() const { return this->deltaT; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setDeltaT(real dt) { deltaT = dt; }
//////////////////////////////////////////////////////////////////////////
bool LBMKernel::getCompressible() const { return compressible; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setCompressible(bool val) { compressible = val; }
//////////////////////////////////////////////////////////////////////////
bool LBMKernel::getWithForcing() const { return withForcing; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setWithForcing(bool val) { withForcing = val; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setBlock(SPtr<Block3D> block) { this->block = block; }
//////////////////////////////////////////////////////////////////////////
SPtr<Block3D> LBMKernel::getBlock() const { return block.lock(); }
//////////////////////////////////////////////////////////////////////////
bool LBMKernel::getWithSpongeLayer() const { return withSpongeLayer; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setWithSpongeLayer(bool val) { withSpongeLayer = val; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setSpongeLayer(const mu::Parser &parser)
{
    this->checkFunction(parser);
    this->muSpongeLayer = parser;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setSpongeLayer(const std::string &muParserString)
{
    this->muSpongeLayer.SetExpr(muParserString);
    this->checkFunction(muSpongeLayer);
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setDataSet(SPtr<DataSet3D> dataSet) { this->dataSet = dataSet; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::swapDistributions() { dataSet->getFdistributions()->swap(); }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setNX(std::array<int, 3> nx) { this->nx = nx; }
//////////////////////////////////////////////////////////////////////////
std::array<int, 3> LBMKernel::getNX() { return nx; }
//////////////////////////////////////////////////////////////////////////
bool LBMKernel::isInsideOfDomain(const int &x1, const int &x2, const int &x3) const
{
    const SPtr<BCArray3D> bcArray = this->bcSet->getBCArray();
    return bcArray->isInsideOfDomain(x1, x2, x3, ghostLayerWidth);
}
//////////////////////////////////////////////////////////////////////////

void LBMKernel::setCollisionFactorMultiphase(real collFactorL, real collFactorG)
{
    this->collFactorL = collFactorL;
    this->collFactorG = collFactorG;
}
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getCollisionFactorL() const { return collFactorL; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getCollisionFactorG() const { return collFactorG; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setDensityRatio(real densityRatio) { this->densityRatio = densityRatio; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getDensityRatio() const { return densityRatio; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setMultiphaseModelParameters(real beta, real kappa)
{
    this->beta  = beta;
    this->kappa = kappa;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::getMultiphaseModelParameters(real &beta, real &kappa)
{
    beta  = this->beta;
    kappa = this->kappa;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setContactAngle(real contactAngle) { this->contactAngle = contactAngle; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getContactAngle() const { return contactAngle; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setPhiL(real phiL) { this->phiL = phiL; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setPhiH(real phiH) { this->phiH = phiH; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getPhiL() const { return phiL; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getPhiH() const { return phiH; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setPhaseFieldRelaxation(real tauH) { this->tauH = tauH; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getPhaseFieldRelaxation() const { return tauH; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setMobility(real mob) { this->mob = mob; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getMobility() const { return mob; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setInterfaceWidth(real w) { this->interfaceWidth = w; }
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getInterfaceWidth() const { return interfaceWidth; }
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setSigma(real sigma){ this->sigma = sigma;}
//////////////////////////////////////////////////////////////////////////
real LBMKernel::getSigma() const { return sigma;}

//! \}

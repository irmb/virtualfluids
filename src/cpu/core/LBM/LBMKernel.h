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

#ifndef LBMKERNEL_H
#define LBMKERNEL_H

#include "ILBMKernel.h"
#include "LBMSystem.h"
#include <PointerDefinitions.h>
#include <array>
#include <limits>
#include <muParser.h>

class BCSet;
class DataSet3D;
class Block3D;

//! \brief A base class provides basic functionality for LBM kernel
class LBMKernel : public ILBMKernel, public enableSharedFromThis<LBMKernel>
{
public:
    LBMKernel();

    virtual SPtr<LBMKernel> clone() = 0;

    void calculate(int step) override    = 0;
    real getCalculationTime() override = 0;

    void setBCSet(SPtr<BCSet> bcp) override;
    SPtr<BCSet> getBCSet() const override;

    void setCollisionFactor(real collFactor) override;
    real getCollisionFactor() const override;

    void setGhostLayerWidth(int witdh);
    int getGhostLayerWidth() const override;

    void setDataSet(SPtr<DataSet3D> dataSet);
    SPtr<DataSet3D> getDataSet() const override;

    void setForcingX1(real forcingX1);
    void setForcingX2(real forcingX2);
    void setForcingX3(real forcingX3);

    void setForcingX1(const mu::Parser &parser);
    void setForcingX2(const mu::Parser &parser);
    void setForcingX3(const mu::Parser &parser);

    void setForcingX1(const std::string &muParserString);
    void setForcingX2(const std::string &muParserString);
    void setForcingX3(const std::string &muParserString);

    void setIndex(int x1, int x2, int x3);

    real getDeltaT() const override;
    void setDeltaT(real dt);

    bool getCompressible() const override;
    void setCompressible(bool val);

    bool getWithForcing() const override;
    void setWithForcing(bool val);

    bool getWithSpongeLayer() const;
    void setWithSpongeLayer(bool val);

    void setSpongeLayer(const mu::Parser &parser);
    void setSpongeLayer(const std::string &muParserString);

    void setBlock(SPtr<Block3D> block);
    SPtr<Block3D> getBlock() const;

    bool isInsideOfDomain(const int &x1, const int &x2, const int &x3) const override;

    void swapDistributions() override;

    void setNX(std::array<int, 3> nx);
    std::array<int, 3> getNX();

    ///////// Extra methods for the multiphase kernel ////////////

    void setCollisionFactorMultiphase(real collFactorL, real collFactorG);
    real getCollisionFactorL() const;
    real getCollisionFactorG() const;
    void setDensityRatio(real densityRatio);
    real getDensityRatio() const;
    void setMultiphaseModelParameters(real beta, real kappa);
    void getMultiphaseModelParameters(real &beta, real &kappa);
    void setContactAngle(real contactAngle);
    real getContactAngle() const;
    void setPhiL(real phiL);
    void setPhiH(real phiH);
    real getPhiL() const;
    real getPhiH() const;
    void setPhaseFieldRelaxation(real tauH);
    real getPhaseFieldRelaxation() const;
    void setMobility(real mob);
    real getMobility() const;
    void setInterfaceWidth(real w);
    real getInterfaceWidth() const;
    void setSigma(real sigma);
    real getSigma() const;

protected:
    SPtr<DataSet3D> dataSet;
    SPtr<BCSet> bcSet;
    real collFactor;
    int ghostLayerWidth{ 1 };
    bool compressible{ false };

    // forcing
    bool withForcing{ false };
    mu::Parser muForcingX1;
    mu::Parser muForcingX2;
    mu::Parser muForcingX3;
    int ix1, ix2, ix3;
    real deltaT{ 1.0 };

    // sponge layer
    bool withSpongeLayer{ false };
    mu::Parser muSpongeLayer;

    WPtr<Block3D> block;

    std::array<int, 3> nx;

    // Multiphase model
    real collFactorL;
    real collFactorG;
    real densityRatio;
    real beta;
    real kappa;
    real sigma;
    real contactAngle;
    real phiL;
    real phiH;
    real tauH;
    real mob;
    real interfaceWidth { 4.0 };

private:
    void checkFunction(mu::Parser fct);
};

#endif

//! \}

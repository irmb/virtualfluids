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
//! \file IBsharpInterfaceLBMKernel.h
//! \ingroup LBMKernel
//! \author M. Geier, K. Kutscher, Hesameddin Safari
//=======================================================================================

#ifndef IBsharpInterfaceLBMKernel_H
#define IBsharpInterfaceLBMKernel_H

#include "BCSet.h"
#include "D3Q27System.h"
#include "LiggghtsCouplingLBMKernel.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"
#include "basics/utilities/UbTiming.h"
#include "IBdynamicsParticleData.h"

//! \brief  Multiphase Cascaded Cumulant LBM kernel.
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  M. Geier, K. Kutscher, Hesameddin Safari
class IBsharpInterfaceLBMKernel : public LiggghtsCouplingLBMKernel
{
public:
    IBsharpInterfaceLBMKernel();
    virtual ~IBsharpInterfaceLBMKernel(void) = default;
    void calculate(int step) override;
    SPtr<LBMKernel> clone() override;

    /// refactor
    // CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr pressure;

    double getCalculationTime() override
    {
        return .0;
    }

protected:
    virtual void initDataSet();
    void swapDistributions() override;

    void initForcing();

    void forwardInverseChimeraWithKincompressible(real &mfa, real &mfb, real &mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho);
    void backwardInverseChimeraWithKincompressible(real &mfa, real &mfb, real &mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho);
    void forwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2);
    void backwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2);

    real f1[D3Q27System::ENDF + 1];

    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr restDistributionsF;

    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsH1;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsH1;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr restDistributionsH1;

    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr pressureOld;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr p1Old;

    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr phaseField;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr phaseFieldOld;
    // CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr divU;

    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr rhoNode;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr vxNode;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr vyNode;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr vzNode;

    real h[D3Q27System::ENDF + 1];
    // real h2[D3Q27System::ENDF + 1];
    // real g  [D3Q27System::ENDF+1];
    real phi[D3Q27System::ENDF + 1];
    real phi2[D3Q27System::ENDF + 1];
    // real pr1[D3Q27System::ENDF+1];
    real phi_cutoff[D3Q27System::ENDF + 1];

    real gradX1_phi();
    real gradX2_phi();
    real gradX3_phi();
    real gradX1_rhoInv(real rhoL, real rhoDIV);
    real gradX2_rhoInv(real rhoL, real rhoDIV);
    real gradX3_rhoInv(real rhoL, real rhoDIV);
    real gradX1_phi2();
    real gradX2_phi2();
    real gradX3_phi2();
    void computePhasefield();
    void findNeighbors(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr ph /*Phase-Field*/, int x1, int x2, int x3);
    void findNeighbors2(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2, int x3);

    real nabla2_phi();

    real computeCurvature_phi();

    mu::value_type muX1, muX2, muX3;
    mu::value_type muDeltaT;
    mu::value_type muNu;
    mu::value_type muRho;
    mu::value_type muPhi;
    real forcingX1;
    real forcingX2;
    real forcingX3;

};

#endif

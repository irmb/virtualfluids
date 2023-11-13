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
//! \file MultiphaseScaleDistributionLBMKernel.h
//! \ingroup LBMKernel
//! \author M. Geier, K. Kutscher, Hesameddin Safari
//=======================================================================================

#ifndef MultiphaseScaleDistributionLBMKernel_H
#define MultiphaseScaleDistributionLBMKernel_H

#include "LBMKernel.h"
#include "BCSet.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

//! \brief  Multiphase Cascaded Cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  M. Geier, K. Kutscher, Hesameddin Safari
class MultiphaseScaleDistributionLBMKernel : public LBMKernel
{
public:
    MultiphaseScaleDistributionLBMKernel();
    virtual ~MultiphaseScaleDistributionLBMKernel(void) = default;
    void calculate(int step) override;
    SPtr<LBMKernel> clone() override;


    ///refactor
    //CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr pressure;


    double getCalculationTime() override { return .0; }
protected:
    virtual void initDataSet();
    void swapDistributions() override;

    void initForcing();

    void forwardInverseChimeraWithKincompressible(real& mfa, real& mfb, real& mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho);
    void backwardInverseChimeraWithKincompressible(real& mfa, real& mfb, real& mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho);
    void forwardChimera(real& mfa, real& mfb, real& mfc, real vv, real v2);
    void backwardChimera(real& mfa, real& mfb, real& mfc, real vv, real v2);

    real f1[D3Q27System::ENDF+1];

    CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
    CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
    CbArray3D<real,IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsF;

    CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr localDistributionsH1;
    CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsH1;
    CbArray3D<real,IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsH1;

    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsH2;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsH2;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsH2;

    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr pressureOld;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr p1Old;

    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr phaseField;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr phaseFieldOld;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr divU; 

    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr rhoNode;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr vxNode;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr vyNode;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr vzNode;

    real h  [D3Q27System::ENDF+1];
    real h2[D3Q27System::ENDF + 1];
    real g  [D3Q27System::ENDF+1];
    real phi[D3Q27System::ENDF+1];
    real phi2[D3Q27System::ENDF + 1];
    real pr1[D3Q27System::ENDF+1];
    real phi_cutoff[D3Q27System::ENDF+1];

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
    void findNeighbors(CbArray3D<real,IndexerX3X2X1>::CbArray3DPtr ph /*Phase-Field*/, int x1, int x2, int x3);
    void findNeighbors2(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2, int x3);
    bool isGas(real phiLim, real* phi, real* phi2);
    bool isGasBoundaryNow(real phiLim, real *phi);

    real nabla2_phi();

    real computeCurvature_phi();

    mu::value_type muX1,muX2,muX3;
    mu::value_type muDeltaT;
    mu::value_type muNu;
    mu::value_type muRho;
    real forcingX1;
    real forcingX2;
    real forcingX3;
};

/// @brief The function computes a fancy expression
/// @param phiLim 
/// @param phi 
/// @param phi2 
/// @return 
inline bool MultiphaseScaleDistributionLBMKernel::isGas(real phiLim, real* phi, real* phi2)
{
    using namespace vf::lbm::dir;
    return (phi2[d000] <= phiLim) || ((phi[d000] <= phiLim) &&
                        (
						(phi[dP00] > phiLim) ||
						(phi[dM00] > phiLim) ||
						(phi[d00P] > phiLim) ||
						(phi[d00M] > phiLim) ||
						(phi[d0M0] > phiLim) ||
						(phi[d0P0] > phiLim) ||
						(phi[dPP0] > phiLim) ||
						(phi[dPM0] > phiLim) ||
						(phi[dP0P] > phiLim) ||
						(phi[dP0M] > phiLim) ||
						(phi[dMP0] > phiLim) ||
						(phi[dMM0] > phiLim) ||
						(phi[dM0P] > phiLim) ||
						(phi[dM0M] > phiLim) ||
						(phi[d0PM] > phiLim) ||
						(phi[d0MM] > phiLim) ||
						(phi[d0PP] > phiLim) ||
						(phi[d0MP] > phiLim) ||
						(phi[dPPP] > phiLim) ||
						(phi[dPMP] > phiLim) ||
						(phi[dMPP] > phiLim) ||
						(phi[dMMP] > phiLim) ||
						(phi[dPPM] > phiLim) ||
						(phi[dPMM] > phiLim) ||
						(phi[dMPM] > phiLim) ||
						(phi[dMMM] > phiLim)
						));
}

inline bool MultiphaseScaleDistributionLBMKernel::isGasBoundaryNow(real phiLim, real *phi)
{
        using namespace vf::lbm::dir;
    return 
           ((phi[d000] <= phiLim) &&
            ((phi[dP00] > phiLim) || (phi[dM00] > phiLim) || (phi[d00P] > phiLim) || (phi[d00M] > phiLim) ||
             (phi[d0M0] > phiLim) || (phi[d0P0] > phiLim) || (phi[dPP0] > phiLim) || (phi[dPM0] > phiLim) ||
             (phi[dP0P] > phiLim) || (phi[dP0M] > phiLim) || (phi[dMP0] > phiLim) || (phi[dMM0] > phiLim) ||
             (phi[dM0P] > phiLim) || (phi[dM0M] > phiLim) || (phi[d0PM] > phiLim) || (phi[d0MM] > phiLim) ||
             (phi[d0PP] > phiLim) || (phi[d0MP] > phiLim) || (phi[dPPP] > phiLim) || (phi[dPMP] > phiLim) ||
             (phi[dMPP] > phiLim) || (phi[dMMP] > phiLim) || (phi[dPPM] > phiLim) || (phi[dPMM] > phiLim) ||
             (phi[dMPM] > phiLim) || (phi[dMMM] > phiLim)));
}

#endif

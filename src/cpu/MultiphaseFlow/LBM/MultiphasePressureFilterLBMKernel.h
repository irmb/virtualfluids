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
//! \file MultiphasePressureFilterLBMKernel.h
//! \ingroup LBMKernel
//! \author M. Geier, K. Kutscher, Hesameddin Safari
//=======================================================================================

#ifndef MultiphasePressureFilterLBMKernel_H
#define MultiphasePressureFilterLBMKernel_H

#include "LBMKernel.h"
#include "BCSet.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

//! \brief  Multiphase Cascaded Cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  M. Geier, K. Kutscher, Hesameddin Safari
class MultiphasePressureFilterLBMKernel : public LBMKernel
{
public:
    MultiphasePressureFilterLBMKernel();
    virtual ~MultiphasePressureFilterLBMKernel(void) = default;
    void calculate(int step) override;
    SPtr<LBMKernel> clone() override;
    real getCalculationTime() override { return .0; }

    void setPhaseFieldBC(real bc)
    {
        phaseFieldBC = bc;
    }
    real getPhaseFieldBC()
    {
        return phaseFieldBC;
    }

protected:
    virtual void initDataSet();
    void swapDistributions() override;

    void initForcing();

    void forwardInverseChimeraWithKincompressible(real& mfa, real& mfb, real& mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho);
    void backwardInverseChimeraWithKincompressible(real& mfa, real& mfb, real& mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho);
    void forwardChimera(real& mfa, real& mfb, real& mfc, real vv, real v2);
    void backwardChimera(real& mfa, real& mfb, real& mfc, real vv, real v2);

    CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
    CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
    CbArray3D<real,IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsF;

    CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr localDistributionsH1;
    CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsH1;
    CbArray3D<real,IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsH1;

    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr pressureOld;

    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr phaseField;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr divU; 

    real h  [D3Q27System::ENDF+1];
    real phi[D3Q27System::ENDF+1];

    real gradX1_phi();
    real gradX2_phi();
    real gradX3_phi();
    void computePhasefield();
    void findNeighbors(CbArray3D<real,IndexerX3X2X1>::CbArray3DPtr ph /*Phase-Field*/, int x1, int x2, int x3);

    real nabla2_phi();

    mu::value_type muX1,muX2,muX3;
    mu::value_type muDeltaT;
    mu::value_type muNu;
    mu::value_type muRho;
    real forcingX1;
    real forcingX2;
    real forcingX3;

    real phaseFieldBC { 0.0 }; // if 0.0 then light fluid on the wall, else if 1.0 havy fluid
};

#endif

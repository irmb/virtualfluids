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
//! \file MultiphaseTwoPhaseFieldsPressureFilterLBMKernel.h
//! \ingroup LBMKernel
//! \author Hesameddin Safari
//=======================================================================================

#ifndef MultiphaseTwoPhaseFieldsPressureFilterLBMKernel_H
#define MultiphaseTwoPhaseFieldsPressureFilterLBMKernel_H

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

//! \brief  Multiphase Cascaded Cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  H. Safari, K. Kutscher, M. Geier
class MultiphaseTwoPhaseFieldsPressureFilterLBMKernel : public LBMKernel
{
public:
   MultiphaseTwoPhaseFieldsPressureFilterLBMKernel();
   virtual ~MultiphaseTwoPhaseFieldsPressureFilterLBMKernel(void) = default;
   void calculate(int step) override;
   SPtr<LBMKernel> clone() override;
   void forwardInverseChimeraWithKincompressible(LBMReal& mfa, LBMReal& mfb, LBMReal& mfc, LBMReal vv, LBMReal v2, LBMReal Kinverse, LBMReal K, LBMReal oneMinusRho);
   void backwardInverseChimeraWithKincompressible(LBMReal& mfa, LBMReal& mfb, LBMReal& mfc, LBMReal vv, LBMReal v2, LBMReal Kinverse, LBMReal K, LBMReal oneMinusRho);
   void forwardChimera(LBMReal& mfa, LBMReal& mfb, LBMReal& mfc, LBMReal vv, LBMReal v2);
   void backwardChimera(LBMReal& mfa, LBMReal& mfb, LBMReal& mfc, LBMReal vv, LBMReal v2);

   ///refactor
   //CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr pressure;
   

   double getCalculationTime() override { return .0; }
protected:
   virtual void initDataSet();
   void swapDistributions() override;
   LBMReal f1[D3Q27System::ENDF+1];

   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsF;

   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributionsH1;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsH1;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsH1;

   CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsH2;
   CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsH2;
   CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsH2;

   //CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr   phaseField;

   CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr pressureOld;

   LBMReal h  [D3Q27System::ENDF+1];
   LBMReal h2[D3Q27System::ENDF + 1];
   LBMReal g  [D3Q27System::ENDF+1];
   LBMReal phi[D3Q27System::ENDF+1];
   LBMReal phi2[D3Q27System::ENDF + 1];
   LBMReal pr1[D3Q27System::ENDF+1];
   LBMReal phi_cutoff[D3Q27System::ENDF+1];

   LBMReal gradX1_phi();
   LBMReal gradX2_phi();
   LBMReal gradX3_phi();
   LBMReal gradX1_phi2();
   LBMReal gradX2_phi2();
   LBMReal gradX3_phi2();
   //LBMReal gradX1_pr1();
   //LBMReal gradX2_pr1();
   //LBMReal gradX3_pr1();
   //LBMReal dirgradC_phi(int n, int k);
   void computePhasefield();
   void findNeighbors(CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr ph /*Phase-Field*/, int x1, int x2, int x3);
   void findNeighbors2(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2, int x3);
   //void findNeighbors(CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr ph /*Phase-Field*/, CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr pf /*Pressure-Field*/, int x1, int x2, int x3);
   //void pressureFiltering(CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr pf /*Pressure-Field*/, CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr pf_filtered /*Pressure-Field*/);

   LBMReal nabla2_phi();


   mu::value_type muX1,muX2,muX3;
   mu::value_type muDeltaT;
   mu::value_type muNu;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;
};

#endif

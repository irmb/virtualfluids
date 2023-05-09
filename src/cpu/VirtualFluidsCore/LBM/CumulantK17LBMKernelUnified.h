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
//! \file CumulantK17LBMKernel.h
//! \ingroup LBM
//! \author Konstantin Kutscher, Martin Geier
//=======================================================================================

#ifndef CumulantK17LBMKernelUnified_h__
#define CumulantK17LBMKernelUnified_h__

#include "LBMKernel.h"
#include "BCSet.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

//! \brief   Compressible cumulant LBM kernel.
//! \details  LBM implementation that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//!
//! The model is publisched in
//! <a href="http://dx.doi.org/10.1016/j.jcp.2017.05.040"><b>[ Geier et al., (2017), 10.1016/j.jcp.2017.05.040]</b></a>,
//! <a href="http://dx.doi.org/10.1016/j.jcp.2017.07.004"><b>[ Geier et al., (2017), 10.1016/j.jcp.2017.07.004]</b></a>
//!
class CumulantK17LBMKernelUnified : public LBMKernel
{
public:
    CumulantK17LBMKernelUnified();
    ~CumulantK17LBMKernelUnified() = default;
    void calculate(int step) override;
    SPtr<LBMKernel> clone() override;
    real getCalculationTime() override { return .0; }

protected:
    virtual void initDataSet();
    real f[D3Q27System::ENDF + 1];

    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr restDistributions;

    mu::value_type muX1, muX2, muX3;
    mu::value_type muDeltaT;
    mu::value_type muNu;
    real forcingX1;
    real forcingX2;
    real forcingX3;
};


#endif // CumulantK17LBMKernel_h__
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
//! \author Konstantin Kutscher, Martin Geier
//=======================================================================================

#ifndef CPU_K17COMPRESSIBLENAVIERSTOKES_H
#define CPU_K17COMPRESSIBLENAVIERSTOKES_H

#include <memory>

#include <basics/DataTypes.h>

#include "LBMKernel.h"

//! \brief   Compressible cumulant LBM kernel.
//! \details  LBM implementation that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//!
//! The model is publisched in
//! <a href="http://dx.doi.org/10.1016/j.jcp.2017.05.040"><b>[ Geier et al., (2017), 10.1016/j.jcp.2017.05.040]</b></a>,
//! <a href="http://dx.doi.org/10.1016/j.jcp.2017.07.004"><b>[ Geier et al., (2017), 10.1016/j.jcp.2017.07.004]</b></a>
//!
class K17CompressibleNavierStokes : public LBMKernel
{
public:
    K17CompressibleNavierStokes();
    void calculate(int step) override;
    std::shared_ptr<LBMKernel> clone() override;
    real getCalculationTime() override
    {
        return .0;
    }

private:
    void initDataSet();
};

#endif

//! \}

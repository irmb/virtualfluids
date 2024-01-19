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
//! \addtogroup NumericalTestPostProcessing
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "TaylorGreenVortexUzLogFileDataImp.h"

std::shared_ptr<TaylorGreenVortexUzLogFileDataImp> TaylorGreenVortexUzLogFileDataImp::getNewInstance()
{
    return std::shared_ptr<TaylorGreenVortexUzLogFileDataImp>(new TaylorGreenVortexUzLogFileDataImp());
}

std::vector<int> TaylorGreenVortexUzLogFileDataImp::getL0()
{
    return l0;
}

std::vector<double> TaylorGreenVortexUzLogFileDataImp::getUz()
{
    return ux;
}

std::vector<double> TaylorGreenVortexUzLogFileDataImp::getAmplitude()
{
    return amp;
}

void TaylorGreenVortexUzLogFileDataImp::setL0(std::vector<int> l0)
{
    this->l0 = l0;
}

void TaylorGreenVortexUzLogFileDataImp::setUz(std::vector<double> ux)
{
    this->ux = ux;
}

void TaylorGreenVortexUzLogFileDataImp::setAmplitude(std::vector<double> amp)
{
    this->amp = amp;
}

TaylorGreenVortexUzLogFileDataImp::TaylorGreenVortexUzLogFileDataImp()
{
}

TaylorGreenVortexUzLogFileDataImp::~TaylorGreenVortexUzLogFileDataImp()
{
}

//! \}

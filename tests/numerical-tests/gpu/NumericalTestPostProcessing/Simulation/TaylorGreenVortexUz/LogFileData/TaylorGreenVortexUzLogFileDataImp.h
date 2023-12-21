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
#ifndef TGV_UZ_LOG_FILE_DATA_IMP_H
#define TGV_UZ_LOG_FILE_DATA_IMP_H

#include "TaylorGreenVortexUzLogFileData.h"

#include <memory>

class TaylorGreenVortexUzLogFileDataImp : public TaylorGreenVortexUzLogFileData
{
public:
    static std::shared_ptr<TaylorGreenVortexUzLogFileDataImp> getNewInstance();

    std::vector<int> getL0();
    std::vector<double> getUz();
    std::vector<double> getAmplitude();

    void setL0(std::vector<int> l0);
    void setUz(std::vector<double> ux);
    void setAmplitude(std::vector<double> amp);

    ~TaylorGreenVortexUzLogFileDataImp();

private:
    TaylorGreenVortexUzLogFileDataImp();

    std::vector<int> l0;
    std::vector<double> ux;
    std::vector<double> amp;
};
#endif
//! \}

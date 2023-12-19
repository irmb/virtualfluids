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
#ifndef VIRTUALFLUIDSPYTHONBINDINGS_D3Q27LBMSYSTEM_H
#define VIRTUALFLUIDSPYTHONBINDINGS_D3Q27LBMSYSTEM_H


#include "AbstractLBMSystem.h"

class D3Q27LBMSystem : public AbstractLBMSystem {
public:
    int getNumberOfDirections() override;

    std::shared_ptr<Interactor3D> makeInteractor() override;

    std::shared_ptr<Interactor3D>
    makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid, int type) override;

    std::shared_ptr<Interactor3D> makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                                                 std::shared_ptr<BC> bcAdapter, int type) override;

    std::shared_ptr<Interactor3D> makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                                                 std::shared_ptr<BC> bcAdapter, int type,
                                                 Interactor3D::Accuracy accuracy) override;
};


#endif //VIRTUALFLUIDSPYTHONBINDINGS_D3Q27LBMSYSTEM_H

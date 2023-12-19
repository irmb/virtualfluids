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
//! \addtogroup cpu_Visitors Visitors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher, Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef RefineCrossAndInsideGbObjectHelper_H
#define RefineCrossAndInsideGbObjectHelper_H

#include <PointerDefinitions.h>
#include <vector>

namespace vf::parallel {class Communicator;}
class Grid3D;
class GbObject3D;

//! \brief Refine blocks on base of bounding boxes.
//! \details You need to use <i>addGbObject()</i> to add corresponding bounding boxes. Then call <i>refine()</i>.
//! \author K. Kutscher
class RefineCrossAndInsideGbObjectHelper
{
public:
    //! Constructor
    //! \param grid a smart pointer to the grid object
    //! \param maxRefineLevel an integer for maximal refinement level
    RefineCrossAndInsideGbObjectHelper(SPtr<Grid3D> grid, int maxRefineLevel, std::shared_ptr<vf::parallel::Communicator> comm);
    virtual ~RefineCrossAndInsideGbObjectHelper();
    //! add geometric object
    //! \param object a smart pointer to bounding box
    //! \param refineLevel a value of refinement level for corresponding bounding box
    void addGbObject(SPtr<GbObject3D> object, int refineLevel);
    //! start refinement
    void refine();

private:
    SPtr<Grid3D> grid;
    std::vector<SPtr<GbObject3D>> objects;
    std::vector<int> levels;
    int maxRefineLevel;
    std::shared_ptr<vf::parallel::Communicator> comm;
};

#endif

//! \}

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
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef CoarsenCrossAndInsideGbObjectBlockVisitor_H
#define CoarsenCrossAndInsideGbObjectBlockVisitor_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class GbObject3D;
class Block3D;
class Grid3D;

//! \brief Refine blocks on base of bounding box which is defined with <i>geoObject</i>
//! \details The class uses a geometry object for define a bounding box. Inside and across this bounding box will be
//! grid on block basis refinement. \author K. Kutscher
class CoarsenCrossAndInsideGbObjectBlockVisitor : public Block3DVisitor
{
public:
    //! A default constructor
    CoarsenCrossAndInsideGbObjectBlockVisitor();
    //! A constructor
    //! \param geoObject a smart pointer to bounding box
    //! \param refineLevel an integer for refine on this level
    CoarsenCrossAndInsideGbObjectBlockVisitor(SPtr<GbObject3D> geoObject, int fineLevel, int coarseLevel);
    ~CoarsenCrossAndInsideGbObjectBlockVisitor() override;
    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;
    //////////////////////////////////////////////////////////////////////////
protected:
    SPtr<GbObject3D> geoObject;
    bool notActive{ true };
    int coarseLevel;
};

#endif

//! \}

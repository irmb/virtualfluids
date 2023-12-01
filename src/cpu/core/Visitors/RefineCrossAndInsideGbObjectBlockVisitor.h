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
//! \file RefineCrossAndInsideGbObjectBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher, Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef RefineCrossAndInsideGbObjectBlockVisitor_H
#define RefineCrossAndInsideGbObjectBlockVisitor_H

#include <PointerDefinitions.h>
#include <vector>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class GbObject3D;

//! \brief Refine blocks on base of bounding box which is defined with <i>geoObject</i>
//! \details The class uses a geometry object for define a bounding box. Inside and across this bounding box will be
//! grid on block basis refinement. \author K. Kucher
class RefineCrossAndInsideGbObjectBlockVisitor : public Block3DVisitor
{
public:
    //! A default constructor
    RefineCrossAndInsideGbObjectBlockVisitor();
    //! A constructor
    //! \param geoObject a smart pointer to bounding box
    //! \param refineLevel an integer for refine on this level
    RefineCrossAndInsideGbObjectBlockVisitor(SPtr<GbObject3D> geoObject, int refineLevel);
    ~RefineCrossAndInsideGbObjectBlockVisitor() override;

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

protected:
    SPtr<GbObject3D> geoObject;
    bool notActive{ true };
};

#endif

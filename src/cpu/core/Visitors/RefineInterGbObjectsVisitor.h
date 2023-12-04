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
//! \file RefineInterGbObjectsBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher, Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef RefineInterGbObjectsVisirtor_H
#define RefineInterGbObjectsVisirtor_H

#include <PointerDefinitions.h>
#include <vector>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class GbObject3D;

//////////////////////////////////////////////////////////////////////////
class RefineInterGbObjectsBlockVisitor : public Block3DVisitor
{
public:
    RefineInterGbObjectsBlockVisitor();
    RefineInterGbObjectsBlockVisitor(SPtr<GbObject3D> includeGbObject3D, SPtr<GbObject3D> excludeGbObject3D,
                                     int startlevel, int stoplevel);
    RefineInterGbObjectsBlockVisitor(std::vector<SPtr<GbObject3D>> includeGbObjects3D,
                                     std::vector<SPtr<GbObject3D>> excludeGbObjects3D, int startlevel, int stoplevel);
    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    std::vector<SPtr<GbObject3D>> includeGbObjects3D;
    std::vector<SPtr<GbObject3D>> excludeGbObjects3D;
};

#endif // RefineInterGbObjectsVisirtor_H

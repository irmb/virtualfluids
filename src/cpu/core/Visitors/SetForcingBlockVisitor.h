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
//! \file SetForcingBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef SetForcingBlockVisitor_h
#define SetForcingBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel.h"

class Block3D;
class Grid3D;

//! \brief Set forcing for all kernels of grid
//! \details This visitor is useful if you need to set or reset forcing in kernels (e.g. after restart because forcing
//! is not serializable). \author K. Kucher
class SetForcingBlockVisitor : public Block3DVisitor
{
public:
    SetForcingBlockVisitor(real forcingX1, real forcingX2, real forcingX3);

    SetForcingBlockVisitor(const mu::Parser &muForcingX1, const mu::Parser &muForcingX2, const mu::Parser &muForcingX3);

    SetForcingBlockVisitor(const std::string &sForcingX1, const std::string &sForcingX2, const std::string &sForcingX3);

    ~SetForcingBlockVisitor() override = default;

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    int ftype;
    real forcingX1;
    real forcingX2;
    real forcingX3;
    mu::Parser muForcingX1;
    mu::Parser muForcingX2;
    mu::Parser muForcingX3;
    std::string sForcingX1;
    std::string sForcingX2;
    std::string sForcingX3;
};

#endif

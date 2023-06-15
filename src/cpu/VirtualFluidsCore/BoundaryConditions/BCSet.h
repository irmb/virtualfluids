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
//! \file BCSet.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef BCSet_H
#define BCSet_H

#include <PointerDefinitions.h>
#include <vector>

class BCArray3D;
class BCStrategy;
class ILBMKernel;

//! A class provides an interface for boundary conditions in the calculation loop.
class BCSet
{
public:
    BCSet();
    BCSet(SPtr<ILBMKernel> kernel);
    virtual ~BCSet();
    virtual SPtr<BCArray3D> getBCArray();
    virtual void setBCArray(SPtr<BCArray3D> bcarray);
    virtual SPtr<BCSet> clone(SPtr<ILBMKernel> kernel);

    void addBC(SPtr<BCStrategy> bc);
    void applyPreCollisionBC();
    void applyPostCollisionBC();
    void clearBC();

protected:
    std::vector<SPtr<BCStrategy>> preBC;
    std::vector<SPtr<BCStrategy>> postBC;
    SPtr<BCArray3D> bcArray;

private:
};

#endif

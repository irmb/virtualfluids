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
//! \file KdSplitAlgorithm.h
//! \ingroup splitalgorithms
//! \author Soeren Textor, Sebastian Bindick
//=======================================================================================
#ifndef KDSPLITALGORITHM_H
#define KDSPLITALGORITHM_H

#include <geometry3d/GbTriFaceMesh3D.h>
//#include <geometry3d/KdTree/Node.h>
#include <geometry3d/KdTree/KdSplitCandidate.h>

#include <vector>

namespace Kd
{
template <typename T>
class Node;

template <typename T>
class SplitAlgorithm
{
public:
    virtual SplitCandidate<T> findBestSplitCandidate(const int &level, const int &maxLevel, Node<T> &node) const   = 0;
    virtual void distributeTriFaces(const SplitCandidate<T> &candidate,
                                    std::vector<GbTriFaceMesh3D::TriFace> &triFacesForChild1,
                                    std::vector<GbTriFaceMesh3D::TriFace> &triFacesForChild2, Node<T> &node) const = 0;
    virtual ~SplitAlgorithm() = default;
};
} // namespace Kd

#endif // KDSPLITALGORITHM_H

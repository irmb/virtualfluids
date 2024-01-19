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
//! \addtogroup geometry3d
//! \ingroup basics
//! \{
//! \author Soeren Textor, Sebastian Bindick
//=======================================================================================
#ifndef KDLINEINTERSECTIONHANDLER_H
#define KDLINEINTERSECTIONHANDLER_H

#include <basics/utilities/UbTuple.h>
//#include <geometry3d/GbTriFaceMesh3D.h>

#include <set>

// #ifdef CAB_RCF
// #  include <3rdParty/rcf/RcfSerializationIncludes.h>
// #end
namespace Kd
{
template <typename T>
class Node;

template <typename T>
class LineIntersectionHandler
{
public:
    virtual bool intersectLine(const UbTuple<T, T, T> &n1, const UbTuple<T, T, T> &n2, Node<T> &parent,
                               Node<T> *&child1, Node<T> *&child2) const = 0;
    virtual ~LineIntersectionHandler()                                   = default;
};
} // namespace Kd

// #if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
//    SF_NO_CTOR(Kd::LineIntersectionHandler<float>);
//    SF_NO_CTOR(Kd::LineIntersectionHandler<double>);
// #endif //RCF_USE_SF_SERIALIZATI
#endif // KDLINEINTERSECTIONHANDLER_H

//! \}

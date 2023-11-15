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
//! \file ArrowImp.cpp
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "ArrowImp.h"

#include "../Vertex/Vertex.h"

 std::shared_ptr<Arrow> ArrowImp::make(const Vertex &start, const Vertex &end)
{
    return std::shared_ptr<ArrowImp>(new ArrowImp(start, end));
}

ArrowImp::ArrowImp(const Vertex &start, const Vertex &end) : start(std::make_shared<Vertex>(start)), end(std::make_shared<Vertex>(end))
{

}

ArrowImp::~ArrowImp()
{

}

std::shared_ptr<Vertex> ArrowImp::getStart() const
{
    return this->start;
}

std::shared_ptr<Vertex> ArrowImp::getEnd() const
{
    return this->end;
}

void ArrowImp::print() const
{
    printf("v1: ");
    start->print();
    printf("v2: ");
    end->print();
}


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
//! \file TriangleException.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz, Martin Schoenherr
//=======================================================================================
#ifndef TriangleException_h
#define TriangleException_h

#include <exception>
#include <iostream>
#include <string>
#include <sstream>

class meshGenExcpetion : public std::exception {
public:
    virtual const char* what() const throw() = 0;
};

class nullVectorImpossibleToCalculateAngle : public meshGenExcpetion
{
    const char* what() const throw() {
        std::ostringstream getNr;
        getNr << "nullVectorImpossibleToCalculateAngle.";
        return getNr.str().c_str();
    }
};

class calculateAngleWhenTrianglesHaveNoCommonEdge : public meshGenExcpetion
{
    const char* what() const throw() {
        std::ostringstream getNr;
        getNr << "Triangles have no common Edge.";
        return getNr.str().c_str();
    }
};

class invalidTriangles : public meshGenExcpetion
{
    const char* what() const throw() {
        std::ostringstream getNr;
        getNr << "Triangles not valid.";
        return getNr.str().c_str();
    }
};

class invalidDelta : public meshGenExcpetion
{
	const char* what() const throw() {
		std::ostringstream getNr;
		getNr << "Delta cant be < Null. To enable no changes change delta to 1.0.";
		return getNr.str().c_str();
	}
};

class compareSameTriangleToFindNeighbor : public meshGenExcpetion
{
	const char* what() const throw() {
		std::ostringstream getNr;
		getNr << "Triangle Container function problem.";
		return getNr.str().c_str();
	}
};

class normalFromTwoAdjacentTrianglesShowInOppositeDirection : public meshGenExcpetion
{
	const char* what() const throw() {
		std::ostringstream getNr;
		getNr << "STL broken, it is not allowed that two adjacent Triangles have a normal that shows in the opposite direction.";
		return getNr.str().c_str();
	}
};


#endif

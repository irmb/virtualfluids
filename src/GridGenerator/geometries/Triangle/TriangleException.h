#ifndef meshGenExcpetion_h
#define meshGenExcpetion_h

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
